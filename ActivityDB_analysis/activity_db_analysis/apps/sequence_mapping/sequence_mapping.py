import pathlib
import re
from typing import Dict, List, Tuple

import pandas as pd
from loguru import logger

from activity_db_analysis.activity_database_interaction.database_names import SequenceTypes
from activity_db_analysis.db.connector import UnifiedPostgresConnector


# ======================================================================================
# Sequence type helpers
# ======================================================================================

def _pdb_sequence_type_order() -> List[SequenceTypes]:
    return [SequenceTypes.atomseq, SequenceTypes.seqres, SequenceTypes.seqresseqadv]


def _pdb_sequence_type_pretty_label(seq_type: SequenceTypes) -> str:
    mapping = {
        "atomseq": "ATOMSEQ",
        "seqres": "SEQRES",
        "seqresseqadv": "SEQADV",
    }
    return mapping.get(seq_type.name.lower(), seq_type.name.upper())


def _value_to_sequence_type() -> Dict[int, SequenceTypes]:
    return {int(member): member for member in SequenceTypes}


# ======================================================================================
# Coverage mode → MV tier mapping
# ======================================================================================

def _map_coverage_mode_to_tier(coverage_mode: str) -> str:
    cm = coverage_mode.lower()
    if cm == "full":
        return "full"
    if cm == "no-hit":
        return "no-hit"
    if cm == "identity-only":
        return "identity-only"
    raise ValueError(f"Unsupported coverage_mode '{coverage_mode}'")


# ======================================================================================
# MV loaders and transformations
# ======================================================================================

def _load_alignment_metrics_mv(
    db_handler: UnifiedPostgresConnector,
    coverage_mode: str,
    chain_level: bool,
) -> pd.DataFrame:
    """
    Load rows for the requested tier from the appropriate alignment metrics MV.

    Sequence-level (default):
        mv_alignment_sequence_metrics
    Chain-level (intersection of chains having ATOMSEQ, SEQRES, SEQADV):
        mv_alignment_sequence_metrics_chain

    Both MVs share the same column schema so downstream logic is identical.
    """
    tier = _map_coverage_mode_to_tier(coverage_mode)
    table_name = "mv_alignment_sequence_metrics_chain" if chain_level else "mv_alignment_sequence_metrics"

    sql = (
        "SELECT tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval, "
        "raw_pairs_interval, distinct_sequences_interval, total_sequences, mapped_sequences_any, "
        'interval_density, interval_occupancy, interval_cond_occupancy, '
        "raw_ge_100, raw_ge_95, raw_ge_90, raw_ge_85, raw_ge_80, "
        "distinct_ge_100, distinct_ge_95, distinct_ge_90, distinct_ge_85, distinct_ge_80 "
        f"FROM {table_name} "
        f"WHERE tier = '{tier}'"
    )
    rows = db_handler.execute(sql) or []
    df = pd.DataFrame(
        rows,
        columns=[
            "tier",
            "chembl_type",
            "pdb_stid",
            "pdb_sequence_type",
            "identity_interval",
            "raw_pairs_interval",
            "distinct_sequences_interval",
            "total_sequences",
            "mapped_sequences_any",
            "interval_density",
            "interval_occupancy",
            "interval_cond_occupancy",
            "raw_ge_100",
            "raw_ge_95",
            "raw_ge_90",
            "raw_ge_85",
            "raw_ge_80",
            "distinct_ge_100",
            "distinct_ge_95",
            "distinct_ge_90",
            "distinct_ge_85",
            "distinct_ge_80",
        ],
    )
    if df.empty:
        logger.warning("No rows found in %s for tier '%s'", table_name, tier)
    return df


def _build_right_bound_tables_from_mv(df_mv: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build right-bound (>= threshold) DataFrames for COMPONENT and VARIANT from MV snapshot.
    We derive RAW / DENSITY / OCC / COND_OCC from the cumulative columns (raw_ge_XX, distinct_ge_XX).
    """
    thresholds = [100, 95, 90, 85, 80]
    col_types = _pdb_sequence_type_order()
    col_labels = [_pdb_sequence_type_pretty_label(t) for t in col_types]

    if df_mv.empty:
        empty_cols = []
        for c in col_labels:
            empty_cols.extend([f"{c}_RAW", f"{c}_DENSITY", f"{c}_OCC", f"{c}_COND_OCC"])
        return pd.DataFrame(columns=empty_cols), pd.DataFrame(columns=empty_cols)

    # Keep one representative row per (chembl_type, pdb_stid) since cumulative metrics repeat per interval.
    collapsed = (
        df_mv.sort_values("identity_interval")
        .groupby(["chembl_type", "pdb_stid", "pdb_sequence_type"], as_index=False)
        .first()
    )

    stid_to_type = _value_to_sequence_type()

    def _side(chembl_type: str) -> pd.DataFrame:
        side_rows = collapsed[collapsed["chembl_type"] == chembl_type]
        if side_rows.empty:
            empty_cols = []
            for c in col_labels:
                empty_cols.extend([f"{c}_RAW", f"{c}_DENSITY", f"{c}_OCC", f"{c}_COND_OCC"])
            return pd.DataFrame(columns=empty_cols)

        records = []
        # Prepare dictionaries for totals only once from >=80 row
        totals_by_type: Dict[SequenceTypes, int] = {}
        mapped_by_type: Dict[SequenceTypes, int] = {}

        for thr in thresholds:
            label = f">={thr}"
            raw_col = f"raw_ge_{thr}"
            distinct_col = f"distinct_ge_{thr}"
            row_values_raw, row_values_density, row_values_occ, row_values_cond = [], [], [], []
            for st in col_types:
                stid = int(st)
                r_match = side_rows[side_rows["pdb_stid"] == stid]
                if r_match.empty:
                    row_values_raw.append(0)
                    row_values_density.append(0.0)
                    row_values_occ.append(0.0)
                    row_values_cond.append(0.0)
                    continue
                r = r_match.iloc[0]
                raw_val = int(r.get(raw_col, 0) or 0)
                distinct_val = float(r.get(distinct_col, 0.0) or 0.0)
                total_sequences = int(r.get("total_sequences", 0) or 0)
                mapped_sequences_any = int(r.get("mapped_sequences_any", 0) or 0)
                density = (raw_val / total_sequences) if total_sequences else 0.0
                occ = (distinct_val / total_sequences) if total_sequences else 0.0
                cond_occ = (distinct_val / mapped_sequences_any) if mapped_sequences_any else 0.0
                if thr == 80:  # capture totals at base floor
                    totals_by_type[stid_to_type[stid]] = total_sequences
                    mapped_by_type[stid_to_type[stid]] = mapped_sequences_any
                row_values_raw.append(raw_val)
                row_values_density.append(density)
                row_values_occ.append(occ)
                row_values_cond.append(cond_occ)
            records.append((label, row_values_raw, row_values_density, row_values_occ, row_values_cond))

        # Build DataFrame
        raw_df = pd.DataFrame(
            [r[1] for r in records],
            index=[r[0] for r in records],
            columns=[f"{c}_RAW" for c in col_labels],
        )
        dens_df = pd.DataFrame(
            [r[2] for r in records],
            index=[r[0] for r in records],
            columns=[f"{c}_DENSITY" for c in col_labels],
        )
        occ_df = pd.DataFrame(
            [r[3] for r in records],
            index=[r[0] for r in records],
            columns=[f"{c}_OCC" for c in col_labels],
        )
        cond_df = pd.DataFrame(
            [r[4] for r in records],
            index=[r[0] for r in records],
            columns=[f"{c}_COND_OCC" for c in col_labels],
        )
        out = pd.concat([raw_df, dens_df, occ_df, cond_df], axis=1)
        out.index.name = "Sequence identity threshold"
        out.columns.name = "Right-bound unique PDB–ChEMBL sequence pair metrics by PDB sequence type"
        # Totals
        out.loc["TOTAL_SEQUENCES", [f"{c}_RAW" for c in col_labels]] = [
            totals_by_type.get(t, 0) for t in col_types
        ]
        out.loc["MAPPED_SEQUENCES", [f"{c}_RAW" for c in col_labels]] = [
            mapped_by_type.get(t, 0) for t in col_types
        ]
        return out

    return _side("COMPONENT"), _side("VARIANT")


def _build_interval_tables_from_mv(df_mv: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build interval bucket tables for COMPONENT and VARIANT from MV snapshot.
    """
    if df_mv.empty:
        return pd.DataFrame(), pd.DataFrame()

    interval_order = ["[100]", "[95,100)", "[90,95)", "[85,90)", "[80,85)"]
    col_types = _pdb_sequence_type_order()
    col_labels = [_pdb_sequence_type_pretty_label(t) for t in col_types]
    stid_to_type = _value_to_sequence_type()

    def _side(chembl_type: str) -> pd.DataFrame:
        side_rows = df_mv[df_mv["chembl_type"] == chembl_type]
        if side_rows.empty:
            return pd.DataFrame()

        # Initialize collectors
        raw_by_type = {t: [0] * len(interval_order) for t in col_types}
        dens_by_type = {t: [0.0] * len(interval_order) for t in col_types}
        occ_by_type = {t: [0.0] * len(interval_order) for t in col_types}
        cond_by_type = {t: [0.0] * len(interval_order) for t in col_types}
        totals_by_type = {t: 0 for t in col_types}
        mapped_by_type = {t: 0 for t in col_types}

        for _, r in side_rows.iterrows():
            interval_label = r["identity_interval"]
            if interval_label not in interval_order:
                continue
            stid = int(r["pdb_stid"])
            seq_type = stid_to_type.get(stid)
            if seq_type not in col_types:
                continue
            idx = interval_order.index(interval_label)
            raw_by_type[seq_type][idx] = int(r["raw_pairs_interval"] or 0)
            dens_by_type[seq_type][idx] = float(r["interval_density"] or 0.0)
            occ_by_type[seq_type][idx] = float(r["interval_occupancy"] or 0.0)
            cond_by_type[seq_type][idx] = float(r["interval_cond_occupancy"] or 0.0)
            totals_by_type[seq_type] = int(r["total_sequences"] or 0)
            mapped_by_type[seq_type] = int(r["mapped_sequences_any"] or 0)

        # Build dataframes
        def _collect(grid):
            return [
                [grid[t][i] for t in col_types]
                for i in range(len(interval_order))
            ]

        df_raw = pd.DataFrame(
            _collect(raw_by_type),
            index=interval_order,
            columns=[f"{c}_RAW" for c in col_labels],
        )
        df_dens = pd.DataFrame(
            _collect(dens_by_type),
            index=interval_order,
            columns=[f"{c}_DENSITY" for c in col_labels],
        )
        df_occ = pd.DataFrame(
            _collect(occ_by_type),
            index=interval_order,
            columns=[f"{c}_OCC" for c in col_labels],
        )
        df_cond = pd.DataFrame(
            _collect(cond_by_type),
            index=interval_order,
            columns=[f"{c}_COND_OCC" for c in col_labels],
        )
        out = pd.concat([df_raw, df_dens, df_occ, df_cond], axis=1)
        out.index.name = "Sequence identity interval"
        out.columns.name = "Interval-bucketed unique PDB–ChEMBL sequence pair metrics by PDB sequence type"
        out.loc["TOTAL_SEQUENCES", [f"{c}_RAW" for c in col_labels]] = [
            totals_by_type[t] for t in col_types
        ]
        out.loc["MAPPED_SEQUENCES", [f"{c}_RAW" for c in col_labels]] = [
            mapped_by_type[t] for t in col_types
        ]
        return out

    return _side("COMPONENT"), _side("VARIANT")


# ======================================================================================
# CSV + generic TeX writer (full metrics)
# ======================================================================================

def _write_csv_and_tex(df: pd.DataFrame, csv_path: pathlib.Path, caption: str, label: str) -> None:
    df.to_csv(csv_path)
    logger.info(f"Wrote {csv_path}")
    try:
        latex_table_body = df.to_latex(index=True, escape=True)
        latex_table_body = re.sub(r'(?m)^(\s*)\[', r'\1\\relax [', latex_table_body)
        latex_table = (
            "\\begin{table}[!htbp]\n"
            "\\tiny\n"
            "\\centering\n"
            f"{latex_table_body}\n"
            f"\\caption{{{caption} This table lists multiple metrics (RAW / DENSITY / OCC / COND\\_OCC).}}\n"
            f"\\label{{tab:{label}}}\n"
            "\\end{table}\n"
        )
        tex_path = csv_path.with_suffix(".tex")
        tex_path.write_text(latex_table, encoding="utf-8")
        logger.info(f"Wrote {tex_path}")
    except Exception as e:
        logger.error(f"Failed to write TeX for {csv_path}: {e}")


# ======================================================================================
# Simplified publication LaTeX tables (interval counts + proportions)
# (Re-used from earlier implementation, adapted to MV-based DataFrames)
# ======================================================================================

def _simplified_interval_tables(
    interval_combined_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    col_map = {
        "COMPONENT_ATOMSEQ_RAW": "COMPONENT_ATOMSEQ_RAW",
        "COMPONENT_SEQRES_RAW": "COMPONENT_SEQRES_RAW",
        "COMPONENT_SEQADV_RAW": "COMPONENT_SEQADV_RAW",
        "VARIANT_ATOMSEQ_RAW": "VARIANT_ATOMSEQ_RAW",
        "VARIANT_SEQRES_RAW": "VARIANT_SEQRES_RAW",
        "VARIANT_SEQADV_RAW": "VARIANT_SEQADV_RAW",
    }
    expected_cols = list(col_map.values())
    missing = [c for c in expected_cols if c not in interval_combined_df.columns]
    if missing:
        raise ValueError(f"Interval combined DataFrame missing expected RAW columns: {missing}")

    interval_rows = [r for r in interval_combined_df.index if r.startswith('[')]
    counts_df = interval_combined_df.loc[interval_rows, expected_cols].copy()

    counts_df.loc["Total"] = counts_df.sum(axis=0)

    totals_sequences: Dict[str, int] = {}
    if "TOTAL_SEQUENCES" in interval_combined_df.index:
        total_row = interval_combined_df.loc["TOTAL_SEQUENCES"]
        for st in ("ATOMSEQ", "SEQRES", "SEQADV"):
            for side in ("COMPONENT", "VARIANT"):
                col = f"{side}_{st}_RAW"
                if col in total_row and pd.notna(total_row[col]):
                    totals_sequences[st] = int(total_row[col])
                    break

    prop_df = counts_df.drop(index="Total").copy()
    for col in prop_df.columns:
        s = prop_df[col].sum()
        prop_df[col] = (prop_df[col] / s * 100.0) if s else 0.0
    prop_df.loc["Sum"] = {c: 100.0 for c in prop_df.columns}

    return counts_df, prop_df, totals_sequences


def _write_simplified_counts_latex(
    counts_df: pd.DataFrame,
    totals_sequences: Dict[str, int],
    out_path: pathlib.Path,
    coverage_mode: str,
    chain_level: bool,
) -> None:
    def _fmt_int(v) -> str:
        return f"{int(v)}"

    lines: List[str] = []
    lines.append(r"\begin{table}[!htbp]")
    lines.append(r"  \centering")
    cap = (
        "Percent identity interval-binned raw counts of unique \\pdb{}-\\chembl sequence pairs"
        f"{' (chain-level)' if chain_level else ''} (coverage mode: {coverage_mode}). "
        "Counts shown separately for Component and Variant by \\pdb sequence type. "
        "Bottom ``Total'' row sums interval counts. "
        f"Unique sequence counts: \\atom = \\num{{{totals_sequences.get('ATOMSEQ',0)}}}; "
        f"\\seqres = \\num{{{totals_sequences.get('SEQRES',0)}}}; "
        f"\\seqadv = \\num{{{totals_sequences.get('SEQADV',0)}}}."
    )
    lines.append(f"  \\caption{{{cap}}}")
    label_suffix = coverage_mode.replace('-', '_') + ("_chain_level" if chain_level else "")
    lines.append(f"  \\label{{tab:target_mapping_interval_counts_{label_suffix}}}")
    lines.append("  \\begin{tabular}{")
    lines.append("    l S[table-format=6] S[table-format=5] S[table-format=5] S[table-format=6] S[table-format=5] S[table-format=4]")
    lines.append("  }")
    lines.append("    \\toprule")
    lines.append("    \\multirow{2}{*}{\\% identity} & \\multicolumn{3}{c}{Component} & \\multicolumn{3}{c}{Variant} \\\\")
    lines.append("    \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}")
    lines.append("    & {ATOM} & {SEQRES} & {SEQADV} & {ATOM} & {SEQRES} & {SEQADV} \\\\")
    lines.append("    \\midrule")

    interval_rows = [r for r in counts_df.index if r.startswith('[')]
    for r in interval_rows:
        lines.append(
            f"    {r:10} & {_fmt_int(counts_df.loc[r,'COMPONENT_ATOMSEQ_RAW']):>6} "
            f"& {_fmt_int(counts_df.loc[r,'COMPONENT_SEQRES_RAW']):>5} "
            f"& {_fmt_int(counts_df.loc[r,'COMPONENT_SEQADV_RAW']):>5} "
            f"& {_fmt_int(counts_df.loc[r,'VARIANT_ATOMSEQ_RAW']):>6} "
            f"& {_fmt_int(counts_df.loc[r,'VARIANT_SEQRES_RAW']):>5} "
            f"& {_fmt_int(counts_df.loc[r,'VARIANT_SEQADV_RAW']):>4} \\\\"
        )

    tot = counts_df.loc["Total"]
    lines.append("    \\midrule")
    lines.append(
        "    Total       &"
        f"{_fmt_int(tot['COMPONENT_ATOMSEQ_RAW']):>6} & {_fmt_int(tot['COMPONENT_SEQRES_RAW']):>5} & {_fmt_int(tot['COMPONENT_SEQADV_RAW']):>5} &"
        f"{_fmt_int(tot['VARIANT_ATOMSEQ_RAW']):>6} & {_fmt_int(tot['VARIANT_SEQRES_RAW']):>5} & {_fmt_int(tot['VARIANT_SEQADV_RAW']):>4} \\\\"
    )
    lines.append("    \\bottomrule")
    lines.append("  \\end{tabular}")
    lines.append("\\end{table}")
    out_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"Wrote simplified interval counts LaTeX table: {out_path}")


def _write_simplified_proportions_latex(
    prop_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    out_path: pathlib.Path,
    coverage_mode: str,
    chain_level: bool,
) -> None:
    example_interval = "[95,100)"
    example_numer = example_denom = example_percent = None
    if example_interval in counts_df.index:
        try:
            numer = counts_df.loc[example_interval, "COMPONENT_ATOMSEQ_RAW"]
            denom = counts_df.loc["Total", "COMPONENT_ATOMSEQ_RAW"]
            if denom:
                example_numer = int(numer)
                example_denom = int(denom)
                example_percent = numer / denom * 100.0
        except Exception as e:
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)

    lines: List[str] = []
    lines.append(r"\begin{table}[!htbp]")
    lines.append(r"  \centering")
    cap_prefix = (
        "Percent identity interval composition (percent of RAW unique \\pdb{}-\\chembl sequence pairs "
        "within each sequence type) for identity $\\geq \\percent{80}$"
        f"{' (chain-level)' if chain_level else ''} (coverage mode: {coverage_mode}). "
        "Normalized stacked bar inputs. "
    )
    if example_numer is not None:
        cap_prefix += (
            f"Example: Component-\\atom{{}} [95,100) = \\num{{{example_numer}}} / "
            f"\\num{{{example_denom}}} = \\percent{{{example_percent:.2f}}}. "
        )
    cap_prefix += "Columns sum to \\percent{100.00}."
    lines.append(f"  \\caption{{{cap_prefix}}}")
    label_suffix = coverage_mode.replace('-', '_') + ("_chain_level" if chain_level else "")
    lines.append(f"  \\label{{tab:target_mapping_interval_proportions_{label_suffix}}}")
    lines.append("  \\begin{tabular}{ l S[table-format=2.2] S[table-format=2.2] S[table-format=2.2] S[table-format=2.2] S[table-format=2.2] S[table-format=2.2] }")
    lines.append("    \\toprule")
    lines.append("    \\multirow{2}{*}{\\% Identity} & \\multicolumn{3}{c}{Component (\\%)} & \\multicolumn{3}{c}{Variant (\\%)} \\\\")
    lines.append("    \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}")
    lines.append("    & {ATOM} & {SEQRES} & {SEQADV} & {ATOM} & {SEQRES} & {SEQADV} \\\\")
    lines.append("    \\midrule")

    interval_rows = [r for r in prop_df.index if r.startswith('[')]
    for r in interval_rows:
        lines.append(
            f"    {r:10} & "
            f"{prop_df.loc[r,'COMPONENT_ATOMSEQ_RAW']:6.2f} & {prop_df.loc[r,'COMPONENT_SEQRES_RAW']:6.2f} & {prop_df.loc[r,'COMPONENT_SEQADV_RAW']:6.2f} & "
            f"{prop_df.loc[r,'VARIANT_ATOMSEQ_RAW']:6.2f} & {prop_df.loc[r,'VARIANT_SEQRES_RAW']:6.2f} & {prop_df.loc[r,'VARIANT_SEQADV_RAW']:6.2f} \\\\"
        )

    lines.append("    \\midrule")
    lines.append("    Sum        &100.00 &100.00 &100.00 &100.00 &100.00 &100.00 \\\\")
    lines.append("    \\bottomrule")
    lines.append("  \\end{tabular}")
    lines.append("\\end{table}")
    out_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"Wrote simplified interval proportions LaTeX table: {out_path}")


def _emit_simplified_interval_publication_tables(
    interval_combined_df: pd.DataFrame,
    output_dir: pathlib.Path,
    coverage_mode: str,
    chain_level: bool,
    filename_prefix: str,
) -> None:
    try:
        counts_df, prop_df, totals_sequences = _simplified_interval_tables(interval_combined_df)
        coverage_mode_filename = coverage_mode.replace("-", "_")
        if chain_level:
            coverage_mode_filename += "_chain_level"
        counts_path = output_dir / f"{filename_prefix}target_matching_interval_counts_{coverage_mode_filename}.tex"
        props_path = output_dir / f"{filename_prefix}target_matching_interval_proportions_{coverage_mode_filename}.tex"
        _write_simplified_counts_latex(
            counts_df=counts_df,
            totals_sequences=totals_sequences,
            out_path=counts_path,
            coverage_mode=coverage_mode,
            chain_level=chain_level,
        )
        _write_simplified_proportions_latex(
            prop_df=prop_df,
            counts_df=counts_df,
            out_path=props_path,
            coverage_mode=coverage_mode,
            chain_level=chain_level,
        )
    except Exception as e:
        logger.error(
            f"Failed to generate simplified publication interval tables for coverage mode '{coverage_mode}' "
            f"(chain_level={chain_level}): {e}"
        )


# ======================================================================================
# Public entry point
# ======================================================================================

def generate_target_matching_reports(
    db_handler: UnifiedPostgresConnector,
    output_file_directory: pathlib.Path,
    chain_level: bool = False,
    coverage_mode: str = "full",
    min_identity_floor: int = 80,  # retained for signature compatibility (MV fixed at 80 floor)
) -> None:
    """
    Generate right-bound (>= thresholds) and interval (disjoint buckets) reports for unique
    PDB–ChEMBL sequence pairs stratified by:
      - ChEMBL sequence classification: COMPONENT vs VARIANT
      - PDB sequence type: ATOMSEQ, SEQRES, SEQADV
      - Identity thresholds: >=100, >=95, >=90, >=85, >=80
      - Disjoint identity intervals: [100], [95,100), [90,95), [85,90), [80,85)

    Data source (pre-aggregated materialized views):
      - Sequence-level (default): mv_alignment_sequence_metrics
      - Chain-level (intersection of chains containing all three PDB sequence types): mv_alignment_sequence_metrics_chain

    Both views share an identical column schema:
        tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval,
        raw_pairs_interval, distinct_sequences_interval, total_sequences, mapped_sequences_any,
        interval_density, interval_occupancy, interval_cond_occupancy,
        raw_ge_100, raw_ge_95, raw_ge_90, raw_ge_85, raw_ge_80,
        distinct_ge_100, distinct_ge_95, distinct_ge_90, distinct_ge_85, distinct_ge_80

    Notes:
      - coverage_mode selects the tier (full | no-hit | identity-only) mapped to MV tier labels.
      - min_identity_floor is retained for API compatibility but the MVs are already fixed at >=80.
      - Chain-level output now truly restricts the cohort; previously it was a placeholder warning.

    Outputs (per invocation):
      1. Right-bound CSV + TeX
      2. Interval CSV + TeX
      3. Simplified publication LaTeX (interval counts + proportions)

    The simplified tables use only RAW counts (plus derived percentages) and feed downstream plotting.
    """
    output_dir = (output_file_directory or pathlib.Path(".")).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Writing target mapping reports to: {output_dir}")
    prefix = "chain_level_" if chain_level else ""
    caption_suffix = " (chain-level)" if chain_level else ""
    label_suffix = "_chain_level" if chain_level else ""
    logger.info(f"Coverage mode: {coverage_mode}  | Chain-level: {'ENABLED' if chain_level else 'DISABLED'}")

    # Single MV load
    mv_df = _load_alignment_metrics_mv(
        db_handler=db_handler,
        coverage_mode=coverage_mode,
        chain_level=chain_level,
    )

    # Right-bound
    logger.info("Deriving right-bound metrics from mv_alignment_sequence_metrics")
    comp_right_df, var_right_df = _build_right_bound_tables_from_mv(mv_df)
    comp_right_df = comp_right_df.add_prefix("COMPONENT_")
    var_right_df = var_right_df.add_prefix("VARIANT_")
    right_combined_df = pd.concat([comp_right_df, var_right_df], axis=1)
    coverage_mode_filename = coverage_mode.replace("-", "_")
    right_path = output_dir / f"{prefix}target_matches_by_sequences_right_bound_{coverage_mode_filename}.csv"
    _write_csv_and_tex(
        df=right_combined_df,
        csv_path=right_path,
        caption=f"Raw, density, occupancy and conditional occupancy metrics for unique PDB–ChEMBL sequence pairs by PDB sequence type (right-bound thresholds){caption_suffix}. Coverage mode: {coverage_mode}.",
        label=f"target_matches_by_sequences_right_bound_{coverage_mode_filename}{label_suffix}",
    )

    # Intervals
    logger.info("Deriving interval metrics from mv_alignment_sequence_metrics")
    comp_int_df, var_int_df = _build_interval_tables_from_mv(mv_df)
    comp_int_df = comp_int_df.add_prefix("COMPONENT_")
    var_int_df = var_int_df.add_prefix("VARIANT_")
    interval_combined_df = pd.concat([comp_int_df, var_int_df], axis=1)
    interval_path = output_dir / f"{prefix}target_matches_by_sequences_sequence_identity_ranges_{coverage_mode_filename}.csv"
    _write_csv_and_tex(
        df=interval_combined_df,
        csv_path=interval_path,
        caption=f"Raw, density, occupancy and conditional occupancy metrics for unique PDB–ChEMBL sequence pairs by PDB sequence type (sequence identity intervals){caption_suffix}. Coverage mode: {coverage_mode}.",
        label=f"target_matches_by_sequences_sequence_identity_ranges_{coverage_mode_filename}{label_suffix}",
    )

    # Simplified publication tables
    logger.info("Generating simplified publication interval tables (counts + proportions).")
    _emit_simplified_interval_publication_tables(
        interval_combined_df=interval_combined_df,
        output_dir=output_dir,
        coverage_mode=coverage_mode,
        chain_level=chain_level,
        filename_prefix=prefix,
    )

    # Summary logging
    try:
        right_tex_path = right_path.with_suffix(".tex")
        interval_tex_path = interval_path.with_suffix(".tex")
        logger.info("Summary of generated target mapping report artifacts:")
        logger.info(f"- Right-bound CSV: {right_path}")
        if right_tex_path.exists():
            logger.info(f"- Right-bound TeX: {right_tex_path}")
        logger.info(f"- Interval CSV:    {interval_path}")
        if interval_tex_path.exists():
            logger.info(f"- Interval TeX:    {interval_tex_path}")

        coverage_mode_filename = coverage_mode.replace("-", "_")
        if chain_level:
            coverage_mode_filename += "_chain_level"
        simp_counts = output_dir / f"{prefix}target_matching_interval_counts_{coverage_mode_filename}.tex"
        simp_props = output_dir / f"{prefix}target_matching_interval_proportions_{coverage_mode_filename}.tex"
        if simp_counts.exists():
            logger.info(f"- Simplified counts TeX: {simp_counts}")
        if simp_props.exists():
            logger.info(f"- Simplified proportions TeX: {simp_props}")
    except Exception as e:
        logger.warning(f"Failed to emit summary of generated report files: {e}")
