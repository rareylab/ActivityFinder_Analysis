"""
LigandMatchLevelCounts
======================

Generates a compact, standalone table summarizing counts of PDB–ChEMBL small‑molecule
matches by matching confidence level (smmclid). The mapping table stores only the
strongest (highest) confidence level per (pdblid, clid) pair, so each row represents
a unique ligand pair.

Column abbreviations (kept minimal for page fit):
  level  ........ matching confidence level (higher = more stringent identifier / structure concordance)
  map (n,% ) .... total ligand pairs at the level (column % across levels)

  P(n,% ) ....... distinct PDB ligand ids at that level (column %)
  C(n,% ) ....... distinct ChEMBL ligand ids at that level (column %)
  P_ex .......... PDB ligand ids exclusive to exactly this level
  C_ex .......... ChEMBL ligand ids exclusive to exactly this level

Global distinct counts for the full table (not repeated in rows):
  P_all (global PDB ligands)  = reported in LaTeX caption
  C_all (global ChEMBL ligands) = reported in LaTeX caption

Percentages are per column across levels; exclusives are absolute.
"""

from __future__ import annotations

import dataclasses
import pathlib
from typing import Optional

import pandas as pd
from loguru import logger

from activity_db_analysis.siteconfig.app_config import AppConfig
from activity_db_analysis.db.connector import UnifiedPostgresConnector


# ---------------------------------------------------------------------------
# SQL QUERIES
# ---------------------------------------------------------------------------

RAW_COUNTS_SQL = """
/*
Per-level aggregation ONLY (decoupled from global / universe counts).

This query yields:
  - mapping_rows: total (pdblid, clid) pairs assigned to the level
  - distinct_pdb_ligands_level / distinct_chembl_ligands_level
  - distinct_pdb_ligands_level_exclusive / distinct_chembl_ligands_level_exclusive
    (ligands whose id appears in exactly ONE confidence level)
  - duplicate_pair_entry_count: sanity check (should be zero)
Global mapping-ligand counts and full ligand universe sizes are now fetched
separately in Python to avoid cross-join induced NULL -> 0 propagation.
*/
WITH base AS (
  SELECT pdblid, clid, smmclid
  FROM pdb_ligand_chembl_ligand
),
pair_level AS (
  SELECT pdblid, clid, smmclid, COUNT(*) AS occurrences
  FROM base
  GROUP BY pdblid, clid, smmclid
),
pdb_level_counts AS (
  SELECT
    smmclid,
    COUNT(DISTINCT pdblid) AS distinct_pdb_ligands_level,
    COUNT(DISTINCT pdblid) FILTER (
      WHERE pdblid IN (
        SELECT pdblid FROM base GROUP BY pdblid HAVING COUNT(DISTINCT smmclid)=1
      )
    ) AS distinct_pdb_ligands_level_exclusive
  FROM base
  GROUP BY smmclid
),
chembl_level_counts AS (
  SELECT
    smmclid,
    COUNT(DISTINCT clid) AS distinct_chembl_ligands_level,
    COUNT(DISTINCT clid) FILTER (
      WHERE clid IN (
        SELECT clid FROM base GROUP BY clid HAVING COUNT(DISTINCT smmclid)=1
      )
    ) AS distinct_chembl_ligands_level_exclusive
  FROM base
  GROUP BY smmclid
)
SELECT
  smmcl.smmclid,
  smmcl.confidence_level,
  smmcl.comment,
  COALESCE(SUM(pl.occurrences), 0) AS mapping_rows,
  COALESCE(plc.distinct_pdb_ligands_level, 0) AS distinct_pdb_ligands_level,
  COALESCE(clc.distinct_chembl_ligands_level, 0) AS distinct_chembl_ligands_level,
  COALESCE(plc.distinct_pdb_ligands_level_exclusive, 0) AS distinct_pdb_ligands_level_exclusive,
  COALESCE(clc.distinct_chembl_ligands_level_exclusive, 0) AS distinct_chembl_ligands_level_exclusive,
  COALESCE(SUM(CASE WHEN pl.occurrences > 1 THEN 1 ELSE 0 END), 0) AS duplicate_pair_entry_count
FROM small_molecule_matching_confidence_level smmcl
LEFT JOIN pair_level pl        ON pl.smmclid = smmcl.smmclid
LEFT JOIN pdb_level_counts plc ON plc.smmclid = smmcl.smmclid
LEFT JOIN chembl_level_counts clc ON clc.smmclid = smmcl.smmclid
GROUP BY
  smmcl.smmclid, smmcl.confidence_level, smmcl.comment,
  plc.distinct_pdb_ligands_level, clc.distinct_chembl_ligands_level,
  plc.distinct_pdb_ligands_level_exclusive, clc.distinct_chembl_ligands_level_exclusive
ORDER BY smmcl.smmclid DESC;
"""

# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class OutputPaths:
    raw_csv: pathlib.Path
    raw_tex: Optional[pathlib.Path]


# ---------------------------------------------------------------------------
# CORE EXTRACTION FUNCTIONS
# ---------------------------------------------------------------------------

def _execute_query(db: UnifiedPostgresConnector, sql: str) -> pd.DataFrame:
    """
    Execute a SELECT-style statement and return a DataFrame.
    Uses the unified connector's raw_connection to obtain a cursor.
    Adds lightweight defensive checks and ignores static type noise
    related to pandas / psycopg2 dynamic attributes.
    """
    rows: list[tuple] | None = None
    columns: list[str] = []
    try:
        raw = getattr(db, "raw_connection", None)  # type: ignore[attr-defined]
        if raw is None:
            logger.error("UnifiedPostgresConnector has no active raw_connection (not connected).")
            return pd.DataFrame()
        with raw.cursor() as cur:  # type: ignore[no-untyped-call]
            cur.execute(sql)  # type: ignore[no-untyped-call]
            rows = cur.fetchall()  # type: ignore[no-untyped-call]
            columns = [d[0] for d in (cur.description or [])]  # type: ignore[index]
    except Exception as e:
        logger.error(f"Query execution failed: {e}")
        return pd.DataFrame()
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns)


def _add_proportion(df: pd.DataFrame, count_column: str, pct_column: str) -> pd.DataFrame:
    """
    Safely add a percentage column (0–100) for count_column relative to its column total.
    Converts potential Decimal values to float to avoid TypeError during arithmetic.
    """
    if count_column not in df.columns:
        df[pct_column] = 0.0
        return df

    # Convert to numeric (float) to avoid Decimal * float issues
    df[count_column] = pd.to_numeric(df[count_column], errors="coerce").fillna(0).astype(float)
    total = float(df[count_column].sum())
    if total == 0:
        df[pct_column] = 0.0
    else:
        df[pct_column] = (df[count_column] / total) * 100.0
    return df


def fetch_raw_counts_with_proportions(db: UnifiedPostgresConnector) -> pd.DataFrame:
    df = _execute_query(db, RAW_COUNTS_SQL)
    if df.empty:
        logger.warning("Raw counts query returned no rows.")
        return df
    # Add proportion columns (percent)
    df = _add_proportion(df, "mapping_rows", "mapping_rows_%")
    # Percentages are column-wise relative; for uniqueness across all levels we only show % for level-local distinct counts, not globals
    if "distinct_pdb_ligands_level" in df.columns:
        df = _add_proportion(df, "distinct_pdb_ligands_level", "distinct_pdb_ligands_level_%")
    if "distinct_chembl_ligands_level" in df.columns:
        df = _add_proportion(df, "distinct_chembl_ligands_level", "distinct_chembl_ligands_level_%")
    # duplicate_pair_entry_count processing removed (column no longer emitted)
    # Fetch global + universe counts independently (decoupled from per-level aggregation)
    global_p = global_c = total_pdblids = total_clids = 0
    try:
        with db.raw_connection.cursor() as cur:  # type: ignore[attr-defined]
            # Global distinct ligands in the mapping table
            cur.execute("SELECT COUNT(DISTINCT pdblid), COUNT(DISTINCT clid) FROM pdb_ligand_chembl_ligand;")
            row = cur.fetchone()
            if row:
                global_p, global_c = int(row[0]), int(row[1])
            # Full ligand universes
            cur.execute("SELECT COUNT(DISTINCT pdblid) FROM pdb_ligand;")
            total_pdblids = int((cur.fetchone() or [0])[0] or 0)
            cur.execute("SELECT COUNT(DISTINCT clid) FROM chembl_ligand;")
            total_clids = int((cur.fetchone() or [0])[0] or 0)
    except Exception as e:  # pragma: no cover
        logger.warning(f"Independent global/universe count queries failed: {e}")

    df.attrs["global_p"] = global_p
    df.attrs["global_c"] = global_c
    df.attrs["total_pdblids"] = total_pdblids
    df.attrs["total_clids"] = total_clids
    # Round selected percentage columns
    for c in ["mapping_rows_%", "distinct_pdb_ligands_level_%", "distinct_chembl_ligands_level_%"]:
        if c in df.columns:
            df[c] = df[c].round(2)

    # Compact presentation: merge count + percent into a single string column to save width
    # We now split count + percentage into separate columns (no parentheses, no % sign in cells).
    def _format_count(series: pd.Series) -> pd.Series:
        return series.astype(int).map(lambda x: f"\\num{{{x}}}")
    def _format_pct(series: pd.Series) -> pd.Series:
        return series.map(lambda x: f"\\num{{{x:.2f}}}")

    total_pdblids_univ = max(df.attrs.get("total_pdblids", 1), 1)
    total_clids_univ = max(df.attrs.get("total_clids", 1), 1)
    cols_dict = {
        "level": df.get("smmclid", pd.Series(range(1, len(df)+1))),
        "map_n": _format_count(df["mapping_rows"]),
        "map_%": _format_pct(df["mapping_rows_%"]),

        "P_n": _format_count(df["distinct_pdb_ligands_level"]),
        "P_%": _format_pct(df["distinct_pdb_ligands_level_%"]),
        "C_n": _format_count(df["distinct_chembl_ligands_level"]),
        "C_%": _format_pct(df["distinct_chembl_ligands_level_%"]),
        # Compute exclusive percentages relative to the total ligand universes
        "P_ex_n": df.get("distinct_pdb_ligands_level_exclusive", pd.Series([0]*len(df))).astype(int).map(lambda x: f"\\num{{{x}}}"),
        "P_ex_%univ": df.get("distinct_pdb_ligands_level_exclusive", pd.Series([0]*len(df))).astype(float).map(lambda x: f"\\num{{{(x/total_pdblids_univ)*100:.2f}}}"),
        "C_ex_n": df.get("distinct_chembl_ligands_level_exclusive", pd.Series([0]*len(df))).astype(int).map(lambda x: f"\\num{{{x}}}"),
        "C_ex_%univ": df.get("distinct_chembl_ligands_level_exclusive", pd.Series([0]*len(df))).astype(float).map(lambda x: f"\\num{{{(x/total_clids_univ)*100:.2f}}}"),
    }
    df_compact = pd.DataFrame(cols_dict)
    # Propagate metadata attributes (global/universe counts) so the caption
    # generation logic downstream can access them. Without this, constructing
    # a new DataFrame dropped attrs, leading to zeros in the LaTeX caption.
    for k in ["global_p", "global_c", "total_pdblids", "total_clids"]:
        if k in df.attrs:
            df_compact.attrs[k] = df.attrs[k]

    return df_compact


# (Removed fetch_best_counts_with_proportions: best-specific aggregation no longer required.)


# ---------------------------------------------------------------------------
# LATEX GENERATION
# ---------------------------------------------------------------------------

def _escape_latex(s: str, preserve_math_subscripts: bool = False) -> str:
    if s is None:
        return ""
    # When preserve_math_subscripts is True we avoid escaping '_' that occur inside math
    # environments delimited by $ ... $, so that expressions like $P_{\\mathrm{ex}}$
    # are not rendered as $P\_{\\mathrm{ex}}$ (which would show the underscore literally).
    if not preserve_math_subscripts:
        repl = {
            "&": r"\&",
            "_": r"\_",
            "%": r"\%",
            "#": r"\#",
            "$": r"\$",
            "~": r"\textasciitilde{}",
            "^": r"\textasciicircum{}",
        }
        for k, v in repl.items():
            s = s.replace(k, v)
        return s
    # Math‑aware escaping:
    repl_non_math = {
        "&": r"\&",
        "%": r"\%",
        "#": r"\#",
        "$": r"\$",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
        # underscore intentionally omitted here (handled contextually)
    }
    out = []
    in_math = False
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "$":
            # Toggle math mode
            in_math = not in_math
            out.append(ch)
            i += 1
            continue
        if in_math:
            # Inside math: leave '_' alone; still escape other specials except braces (macros use them)
            if ch in repl_non_math:
                out.append(repl_non_math[ch])
            else:
                out.append(ch)
        else:
            # Outside math: escape underscore and other specials
            if ch == "_":
                out.append(r"\_")
            elif ch in repl_non_math:
                out.append(repl_non_math[ch])
            else:
                out.append(ch)
        i += 1
    return "".join(out)


def dataframe_to_latex_table(
    df: pd.DataFrame,
    caption: str,
    label: str,
    column_formats: Optional[str] = None,
    float_format: str = "{:.2f}",
) -> str:
    """
    Build a LaTeX table using booktabs + row coloring.
    """
    if df.empty:
        return (
            "\\begin{table}[!h]\n"
            "\\footnotesize\n"
            "\\centering\n"
            f"\\caption{{{_escape_latex(caption)} (No data returned.)}}\n"
            f"\\label{{tab:{_escape_latex(label)}}}\n"
            "\\begin{tabular}{l}\n\\toprule\nNo data \\\\\n\\bottomrule\n\\end{tabular}\n"
            "\\end{table}\n"
        )

    # Determine grouped headers (pairs *_n / *_% or *_%univ)
    cols = list(df.columns)
    groups = []
    i = 0
    while i < len(cols):
        c = cols[i]
        if c.endswith("_n") and i + 1 < len(cols):
            nxt = cols[i+1]
            if nxt.startswith(c[:-2]) and (nxt.endswith("_%") or nxt.endswith("_%univ")):
                groups.append((c[:-2], [c, nxt]))
                i += 2
                continue
        if c == "level":
            groups.append(("level", [c]))
        else:
            groups.append((c, [c]))
        i += 1
    if column_formats is None:
        fmts = ["r" if g[0] == "level" else "r" * len(g[1]) for g in groups]
        column_formats = "".join(fmts)
    # First header row (group names)
    header_row1_cells = []
    for name, members in groups:
        if name == "level":
            header_row1_cells.append("\\textbf{level}")
        else:
            # Render group headers with math subscripts for *_ex groups
            if name in {"P_ex", "C_ex"}:
                base = name.split('_')[0]
                macro = "$\\mathrm{{{}}}_{{\\mathrm{{ex}}}}$".format(base)
                header_row1_cells.append("\\multicolumn{{{}}}{{c}}{{\\textbf{{{}}}}}".format(len(members), macro))
            elif name in {"P", "C"}:
                macro = "$\\mathrm{{{}}}$".format(name)
                header_row1_cells.append("\\multicolumn{{{}}}{{c}}{{\\textbf{{{}}}}}".format(len(members), macro))
            else:
                header_row1_cells.append("\\multicolumn{{{}}}{{c}}{{\\textbf{{{}}}}}".format(len(members), _escape_latex(name)))
    header_row1 = " & ".join(header_row1_cells) + " \\\\"
    # Second header row (sub-columns)
    sub_header_cells = []
    for name, members in groups:
        if name == "level":
            sub_header_cells.append("")  # blank under 'level'
        else:
            for m in members:
                if m.endswith("_n"):
                    sub_header_cells.append("\\textbf{n}")
                elif m.endswith("_%univ"):
                    sub_header_cells.append("\\textbf{\\%univ}")
                elif m.endswith("_%"):
                    sub_header_cells.append("\\textbf{\\%}")
                else:
                    sub_header_cells.append(f"\\textbf{{{_escape_latex(m)}}}")
    header_row2 = " & ".join(sub_header_cells) + " \\\\"

    # Build header cells explicitly to avoid any nested f-string escaping side-effects
    # (Old single-row header removed; replaced by two-level header above)

    body_lines = []
    for _, row in df.iterrows():
        row_cells = []
        for c in cols:
            val = row[c]
            if isinstance(val, float):
                if pd.isna(val):
                    row_cells.append("")
                else:
                    row_cells.append(float_format.format(val))
            else:
                row_cells.append(_escape_latex(str(val)))
        body_lines.append(" & ".join(row_cells) + " \\\\")

    table = [
        "\\begin{table}[!h]",
        "\\footnotesize",
        "\\centering",
        f"\\caption{{{_escape_latex(caption, preserve_math_subscripts=True)}}}",
        f"\\label{{tab:{_escape_latex(label)}}}",
        "\\rowcolors{2}{gray!15}{white}",
        f"\\begin{{tabular}}{{{column_formats}}}",
        "\\toprule",
        header_row1,
        header_row2,
        "\\midrule",
        *body_lines,
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}",
        "",
    ]
    return "\n".join(table)


def write_latex_tables(
    raw_df: pd.DataFrame,
    raw_tex_path: pathlib.Path,
) -> None:
    raw_global_p = raw_df.attrs.get("global_p", 0)
    raw_global_c = raw_df.attrs.get("global_c", 0)
    total_pdblids = raw_df.attrs.get("total_pdblids", 0)
    total_clids = raw_df.attrs.get("total_clids", 0)
    raw_caption = (
        "PDB–ChEMBL small‑molecule mapping summary by confidence level (smmclid; higher = stricter structural / identifier agreement). "
        "Each row lists a unique best (pdblid, clid) pair placed at its highest confidence level after canonical normalization. "
        "Columns are grouped: map, $\\mathrm{P}$, $\\mathrm{C}$, $\\mathrm{P}_{\\mathrm{ex}}$, $\\mathrm{C}_{\\mathrm{ex}}$ each split into a count subcolumn (n) and a percentage subcolumn (% or %univ). "
        "All map counts are unique (pdblid, clid) pairs with no duplication across or within levels. "
        "$\\mathrm{P}_{\\mathrm{ex}}$ / $\\mathrm{C}_{\\mathrm{ex}}$ percentages are relative to the ligand universes (PDB = \\num{" + f"{total_pdblids}" + "}, ChEMBL = \\num{" + f"{total_clids}" + "}). "
        "Global unique ligand totals: P_all = \\num{" + f"{raw_global_p}" + "}, C_all = \\num{" + f"{raw_global_c}" + "}. "
        "Column percentages (except *_ex) are computed across levels within the column; exclusives use total ligand universes. "
        "Data source: pdb_ligand_chembl_ligand (post canonicalization retaining only highest-confidence matches)."
    )
    raw_table = dataframe_to_latex_table(
        raw_df,
        caption=raw_caption,
        label="ligand_match_levels",
    )
    raw_tex_path.write_text(raw_table)
    logger.info(f"Wrote LaTeX table: {raw_tex_path}")


# ---------------------------------------------------------------------------
# OUTPUT HANDLING
# ---------------------------------------------------------------------------

def _resolve_output_dir(output_dir: Optional[pathlib.Path]) -> pathlib.Path:
    if output_dir is None:
        # Default to tables subdirectory in current working directory
        output_dir = pathlib.Path.cwd() / "tables"
    output_dir = output_dir.expanduser().absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def build_output_paths(
    output_dir: pathlib.Path,
    raw_basename: str = "ligand_match_levels.csv",
    write_tex: bool = True,
) -> OutputPaths:
    raw_csv = output_dir / raw_basename
    raw_tex = raw_csv.with_suffix(".tex") if write_tex else None
    return OutputPaths(raw_csv=raw_csv, raw_tex=raw_tex)


# ---------------------------------------------------------------------------
# PUBLIC RUN FUNCTION
# ---------------------------------------------------------------------------

def run(
    output_dir: Optional[pathlib.Path] = None,
    raw_basename: str = "ligand_match_levels.csv",
    write_tex: bool = True,
    app_config: Optional[pathlib.Path] = None,
) -> None:
    """
    Extract ligand mapping confidence level counts (single aggregation; data already best-per-pair)
    and write CSV (and optional LaTeX).

    Parameters
    ----------
    output_dir : Optional[pathlib.Path]
        Destination directory for CSV/TeX. Defaults to thesis tables directory.
    raw_basename : str
        Filename for raw mapping counts CSV.
    write_tex : bool
        Whether to also emit LaTeX table files.
    app_config : Optional[pathlib.Path]
        Explicit path to appconfig.toml (otherwise assume already loaded).
    """
    # Load explicit app config if provided
    if app_config is not None:
        from activity_db_analysis.siteconfig.app_config import parse_app_config_toml
        cfg = parse_app_config_toml(app_config)
        if isinstance(cfg, dict):
            AppConfig.set_config_file_parameters(cfg)
        else:
            logger.warning(f"Could not parse configuration at {app_config}")

    out_dir = _resolve_output_dir(output_dir)
    paths = build_output_paths(
        output_dir=out_dir,
        raw_basename=raw_basename,
        write_tex=write_tex,
    )

    from activity_db_analysis.siteconfig.app_config import AppConfig as _AC
    with _AC.build_activity_sql_handler_from_config() as db:
        logger.info("Fetching ligand matching level counts...")
        raw_df = fetch_raw_counts_with_proportions(db)
        raw_df.to_csv(paths.raw_csv, index=False)
        logger.info(f"Wrote CSV file: {paths.raw_csv}")
        if write_tex and paths.raw_tex:
            write_latex_tables(
                raw_df=raw_df,
                raw_tex_path=paths.raw_tex,
            )


# ---------------------------------------------------------------------------
# OPTIONAL: SIMPLE MAIN (manual execution)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Minimal manual run (relies on AppConfig already loaded if needed)
    run()
