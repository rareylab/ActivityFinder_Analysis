#!/usr/bin/env python3
"""
PairIdentityDistributionActivity
================================

Activity-qualified, coverage-agnostic pair identity distribution.

(Updated) Implementation now sources data from the staged ligand‑first
pipeline culminating in the materialized view
`mv_ligand_supported_pair_best_identity` (replacing the older
`mv_activity_supported_pair_best_identity`). The staged chain enforces:

  (pdblid → clid) → activities (aid) → assays (asid) → targets (tid / chembl_tid)
  → target sequences (sid) → best BLAST alignments (query_sid, hit_sid, identity≥80)
  → structural remap back to (pdblid, clid).

The final MV stores one row per qualifying (pdblid, clid) with:
  * max_percent_identity (>=80)
  * best_smmclid
  * contributing_alignment_rows
  * pre-bucketed identity_interval

Derived (post-query) field added here:
  * confidence_level (compressed 1..5 scale mapped from best_smmclid where smmclid 6 → confidence_level 5;
    historical smmclid 5 retired). This column is intended for plotting / legends so that
    structural confidence levels appear visually as 1,2,3,4,5 instead of 1,2,3,4,6.

The analysis now performs only lightweight ordering and percentage derivations.

Purpose
-------
Compute the distribution of distinct ligand mapping pairs (pdblid, clid) over
maximum BLAST percent identity intervals (mathematical interval notation) while
REQUIRING that each pair is supported by at least one ChEMBL activity row of
*any* kind (no pchembl_value filter). Structural confidence (smmclid) is
summarized as the highest (best) structural matching level observed for that
pair. Coverage filters (query_coverage / hit_coverage) are intentionally
ignored; the only sequence alignment constraint is percent_identity >= 80.

Identity Intervals (mutually exclusive, descending):
    [100]
    [95,100)
    [90,95)
    [85,90)
    [80,85)

Outputs
-------
1. CSV:  pair_identity_distribution_activity.csv
2. LaTeX: pair_identity_distribution_activity.tex (unless --no-tex)
(In cumulative mode: pair_identity_distribution_activity_cumulative.csv / .tex)
3. Optional derived metrics (added in-memory, not persisted in DB):
   - interval_total
   - pct_within_interval
   - pct_global

SQL Source
----------
Uses the constant PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_SQL defined in
ligand_mapping_sql.py. (Refactored)
  * NOW: Reads directly from final MV `mv_ligand_supported_pair_best_identity`
    (refresh staged chain first) – one pre-aggregated row per (pdblid, clid)
  * Avoids raw BLAST / target / assay / activity joins at analysis runtime
  * Provides: max_percent_identity, best_smmclid, contributing_alignment_rows, identity_interval
  * Buckets are mutually exclusive; each pair appears exactly once

Assumptions
-----------
- Database connection parameters are supplied via CLI or loaded configuration.
- Required base tables and indexes exist (see ActivityDB schema docs).
- The SQL constant reflects the corrected join path (assay.tid = target.tid).
- Counts may be lower than coverage-filtered universes (because we enforce
  activity presence).

CLI Integration
---------------
Add to your main CLI (example snippet):

    @app.command()
    def pair_identity_distribution_activity(...):
        from activity_db_analysis.apps.ligand_mapping.pair_identity_distribution_activity import run
        run(...)

Example Invocation
------------------
    uv run activity_db_analysis/cli.py pair-identity-distribution-activity \
        --output-dir ./tables

"""

from __future__ import annotations

import pathlib
from typing import Optional

import pandas as pd
from loguru import logger

from activity_db_analysis.db.connector import UnifiedPostgresConnector
from activity_db_analysis.siteconfig.app_config import AppConfig

from .latex_utils import build_latex_table
from .ligand_mapping_sql import (
    PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_CUMULATIVE_SQL,
    PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_SQL,
)

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

DEFAULT_OUTPUT_DIR = pathlib.Path.cwd() / "tables"

INTERVAL_ORDER = ["[100]", "[95,100)", "[90,95)", "[85,90)", "[80,85)"]


# --------------------------------------------------------------------------------------
# Internal Helpers
# --------------------------------------------------------------------------------------


def _resolve_output_dir(out: Optional[pathlib.Path]) -> pathlib.Path:
    out = (out or DEFAULT_OUTPUT_DIR).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    return out


# Connection helper removed; unified configuration + environment variables are enforced.


def _run_distribution_query(
    db: UnifiedPostgresConnector, cumulative: bool
) -> pd.DataFrame:
    sql = (
        PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_CUMULATIVE_SQL
        if cumulative
        else PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_SQL
    )
    logger.info(
        "Executing {}activity-qualified pair identity distribution SQL...",
        "cumulative " if cumulative else "",
    )
    rows = db.execute(sql) or []
    if not rows:
        logger.warning("Pair identity distribution query returned 0 rows")
        return pd.DataFrame()
    # Derive columns via a lightweight probe
    raw = getattr(db, "raw_connection", None)
    cols: list[str] = []
    if raw is not None:
        try:
            with raw.cursor() as cur:  # type: ignore[no-untyped-call]
                cur.execute(sql)  # type: ignore[no-untyped-call]
                cols = [d[0] for d in (cur.description or [])]  # type: ignore[index]
        except Exception:
            cols = [f"col_{i}" for i in range(len(rows[0]))]
    else:
        cols = [f"col_{i}" for i in range(len(rows[0]))]
    df = pd.DataFrame(rows, columns=cols)
    logger.info("Query returned {} rows", len(df))
    return df


def _augment(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add interval totals, global totals, percentages; enforce interval ordering.
    """
    if df.empty:
        return df

    # Coerce numeric
    df["pair_count"] = (
        pd.to_numeric(df["pair_count"], errors="coerce").fillna(0).astype(int)
    )
    df["smmclid"] = pd.to_numeric(df["smmclid"], errors="coerce").astype(int)

    # Categorical ordering
    df["identity_interval"] = pd.Categorical(
        df["identity_interval"], categories=INTERVAL_ORDER, ordered=True
    )

    # Interval totals
    interval_totals = (
        df.groupby("identity_interval")["pair_count"].sum().rename("interval_total")
    )
    df = df.merge(interval_totals, on="identity_interval", how="left")

    # Global total
    global_total = float(df["pair_count"].sum())
    if global_total == 0:
        df["pct_within_interval"] = 0.0
        df["pct_global"] = 0.0
    else:
        df["pct_within_interval"] = (
            df["pair_count"] / df["interval_total"] * 100.0
        ).round(2)
        df["pct_global"] = (df["pair_count"] / global_total * 100.0).round(3)

    # Prune retired structural confidence level 5 (Standard InChI) if present
    if "smmclid" in df.columns and (df["smmclid"] == 5).any():
        before = len(df)
        df = df[df["smmclid"] != 5]
        removed = before - len(df)
        logger.info("Pruned {} rows with retired smmclid=5 (Standard InChI).", removed)

    # Derive compressed confidence_level 1..5 (retiring gap at smmclid 5):
    #   smmclid: 1 2 3 4 6
    #   conf  : 1 2 3 4 5
    CONFIDENCE_LEVEL_MAP = {1: 1, 2: 2, 3: 3, 4: 4, 6: 5}
    df["confidence_level"] = df["smmclid"].map(CONFIDENCE_LEVEL_MAP).astype(int)

    # Sort final using the compressed confidence_level for presentation
    df = df.sort_values(["identity_interval", "confidence_level", "smmclid"])
    return df


def _build_latex(df: pd.DataFrame) -> str:
    """
    Construct LaTeX table (auto-format numeric columns) using existing generic utility.
    """
    caption = (
        "Activity-qualified distribution of distinct ligand mapping pairs "
        "$(\\mathrm{pdblid},\\mathrm{clid})$ across maximum BLAST sequence identity "
        "intervals (coverage ignored). Each pair appears exactly once, assigned via its "
        "maximum observed percent identity (>=80) among activity-supported sequence "
        "mappings; structural confidence (smmclid) is the highest structural matching "
        "level observed. Percent columns: within-interval and global."
    )
    label = "tab:pair_identity_distribution_activity"
    header_aliases = {
        "identity_interval": "Identity Interval",
        "confidence_level": "Confidence Level",
        "smmclid": "smmclid",
        "pair_count": "Pairs",
        "interval_total": "Interval Total",
        "pct_within_interval": "\\% Interval",
        "pct_global": "\\% Global",
    }
    # Present confidence_level first (compressed 1..5), keep raw smmclid for traceability.
    ordered_cols = [
        "identity_interval",
        "confidence_level",
        "smmclid",
        "pair_count",
        "interval_total",
        "pct_within_interval",
        "pct_global",
    ]
    present = [c for c in ordered_cols if c in df.columns]
    # Wrap identity_interval values in math mode parentheses for LaTeX
    # Example: [100] -> \([100]\) so that brackets render inside math mode.
    df_for_latex = df[present].copy()
    if "identity_interval" in df_for_latex.columns:

        def _wrap_interval(val):
            sval = str(val).strip()
            # Already properly math-wrapped
            if sval.startswith(r"\(") and sval.endswith(r"\)"):
                return sval
            # Wrap in math mode so LaTeX renders brackets correctly: [100] -> \([100]\)
            return r"\(" + sval + r"\)"

        df_for_latex["identity_interval"] = df_for_latex["identity_interval"].map(
            _wrap_interval
        )

    return build_latex_table(
        df_for_latex,
        caption=caption,
        label=label,
        header_aliases=header_aliases,
        auto_format=True,
        squeeze=True,
    )


# --------------------------------------------------------------------------------------
# Public Entry Point
# --------------------------------------------------------------------------------------


def run(
    output_dir: Optional[pathlib.Path] = None,
    no_tex: bool = False,
    cumulative: bool = False,
    app_config: Optional[pathlib.Path] = None,
) -> None:
    """
    Public entry point: generate pair identity distribution (activity-qualified) CSV (+ optional LaTeX).
    Connection overrides removed; unified AppConfig + env password policy enforced.
    Execute the activity-qualified pair identity distribution analysis.

    Parameters
    ----------
    output_dir : Optional[pathlib.Path]
        Where to write outputs (CSV + optional LaTeX).
    no_tex : bool
        If True, skip LaTeX generation.
    cumulative : bool
        If True, use multi-bucket (cumulative) MV so a pair contributes to every
        identity interval for which it has at least one qualifying alignment. Output
        filenames gain suffix '_cumulative'.
    app_config : Optional[pathlib.Path]
        Path to application config file. If not provided, uses default location.

    Notes
    -----
    SSH tunnel management is handled at the CLI layer. Database connections
    use the unified AppConfig mechanism with environment variable password sourcing.
    """
    logger.info("=" * 90)
    logger.info("Starting Activity-Qualified Pair Identity Distribution")
    logger.info("=" * 90)

    out_dir = _resolve_output_dir(output_dir)
    # Unified connection context (no per-command overrides)
    from activity_db_analysis.siteconfig.app_config import AppConfig as _AC

    with _AC.build_activity_sql_handler_from_config() as db:
        try:
            df = _run_distribution_query(db, cumulative=cumulative)
            aug_df = _augment(df)

            suffix = "_cumulative" if cumulative else ""
            csv_path = out_dir / f"pair_identity_distribution_activity{suffix}.csv"
            aug_df.to_csv(csv_path, index=False)
            logger.info("Wrote CSV: {}", csv_path)

            if not no_tex:
                tex_path = out_dir / f"pair_identity_distribution_activity{suffix}.tex"
                tex_content = _build_latex(aug_df)
                tex_path.write_text(tex_content)
                logger.info("Wrote LaTeX: {}", tex_path)
            else:
                logger.info("Skipping LaTeX output (--no-tex).")

            logger.info(
                "Completed {}activity-qualified pair identity distribution.",
                "cumulative " if cumulative else "",
            )

        except Exception as exc:
            logger.error("Pair identity distribution failed: {}", exc)
            raise
        finally:
            try:
                if hasattr(db, "connection") and db.connection:
                    db.connection.close()
                    logger.info("Database connection closed.")
            except Exception as e:
                logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)


if __name__ == "__main__":
    run()
