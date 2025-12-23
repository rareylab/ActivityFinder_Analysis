#!/usr/bin/env python3
"""
sequence_mapping_parity_check.py
================================

Purpose
-------
Runtime parity / regression safety net that compares:

1. Legacy (inline) best-hit + coverage-filtered logic
2. Unified tier-view based logic (v_best_alignment_id*, mv_best_alignment_identity_floor)

for both:
  - Right-bound thresholds (>= 100,95,90,85,80)
  - Interval buckets ([100], [95,100), [90,95), [85,90), [80,85))

Dimensions validated:
  - chembl_type (COMPONENT / VARIANT)
  - pdb_stid (sequence type id)
  - pair_count_raw (multi-mapping raw pair rows)
  - distinct_sequences (distinct pdb_sid)

Coverage modes validated (mapped to old/new):
  identity-only  -> v_best_alignment_id
  no-hit         -> v_best_alignment_id_query
  full           -> v_best_alignment_id_query_hit

If ANY discrepancy is detected the script prints a detailed diff report
and exits with non-zero status (unless --allow-differences is supplied).

Usage Examples
--------------
# Basic (Activity DB / ChEMBL config already loaded elsewhere)
python -m activity_db_analysis.apps.sequence_mapping.sequence_mapping_parity_check --coverage-modes full no-hit identity-only

# With explicit connection overrides (falls back to defaults if omitted)
python -m activity_db_analysis.apps.sequence_mapping.sequence_mapping_parity_check \
   --host localhost --port 5432 --user myuser --database-name activitydb \
   --coverage-modes full identity-only

# Permit differences (informational run only)
python -m activity_db_analysis.apps.sequence_mapping.sequence_mapping_parity_check --allow-differences

Exit Codes
----------
0 = Perfect parity (or differences allowed)
1 = Differences detected (and not allowed)
2 = Runtime / connection / SQL error

Notes
-----
- This script re-implements the *legacy* logic purely for comparison; it does
  NOT rely on the refactored PDBChEMBLTargetMatchesExpressions builder.
- Chain-level parity can be added later if required (the unified MV currently
  does not include chain-level intersection logic; enabling it here without
  a legacy counterpart would bias results).

"""

from __future__ import annotations

import argparse
import sys
import textwrap
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

from loguru import logger

# Reuse existing connector abstraction
try:
    from activity_db_analysis.db.connector import (
        UnifiedPostgresConnector,
        build_unified_connector,
    )
except Exception as e:  # pragma: no cover
    logger.error(
        "Failed to import database connector. Ensure project environment is set up. Error: {}",
        e,
    )
    sys.exit(2)


RIGHT_THRESHOLDS: List[int] = [100, 95, 90, 85, 80]
INTERVALS: List[Tuple[int, int]] = [
    (100, 101),  # Represent [100] using left=100, right exclusive=101
    (95, 100),
    (90, 95),
    (85, 90),
    (80, 85),
]

COVERAGE_MODE_TO_TIER_VIEW = {
    "identity-only": "v_best_alignment_id",
    "no-hit": "v_best_alignment_id_query",
    "full": "v_best_alignment_id_query_hit",
}

VALID_COVERAGE_MODES = list(COVERAGE_MODE_TO_TIER_VIEW.keys())

# --------------------------------------------------------------------------------------
# Data Structures
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class MetricKey:
    coverage_mode: str
    mode_type: str  # 'right' or 'interval'
    threshold_or_interval: str  # e.g. '>=95' or '[90,95)'
    chembl_type: str
    pdb_stid: int


@dataclass
class MetricRecord:
    pair_count_raw: int
    distinct_sequences: int


# --------------------------------------------------------------------------------------
# Legacy SQL builders (inline replication of pre-refactor logic)
# --------------------------------------------------------------------------------------


def _legacy_coverage_predicate(mode: str) -> str:
    m = mode.lower()
    if m == "full":
        return (
            "bm.query_coverage IS NOT NULL "
            "AND bm.hit_coverage IS NOT NULL "
            "AND bm.query_coverage >= 80.0 "
            "AND bm.hit_coverage >= 80.0"
        )
    if m == "no-hit":
        return (
            "bm.query_coverage IS NOT NULL "
            "AND bm.hit_coverage IS NOT NULL "
            "AND bm.query_coverage >= 80.0"
        )
    if m == "identity-only":
        return "bm.query_coverage IS NOT NULL AND bm.hit_coverage IS NOT NULL"
    raise ValueError(f"Unsupported coverage mode '{mode}'")


def build_legacy_right_bound_sql(threshold: int, coverage_mode: str) -> str:
    coverage_sql = _legacy_coverage_predicate(coverage_mode)
    return f"""
    WITH chembl_hit AS (
        SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
        UNION ALL
        SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
    ),
    ranked AS (
        SELECT
            ps.sid AS pdb_sid,
            st.stid AS pdb_stid,
            ch.chembl_type,
            ch.hit_sid,
            bm.percent_identity,
            bm.query_coverage,
            bm.hit_coverage,
            ROW_NUMBER() OVER (
                PARTITION BY ps.sid, st.stid, ch.chembl_type, ch.hit_sid
                ORDER BY bm.percent_identity DESC, bm.query_coverage DESC, bm.hit_coverage DESC
            ) AS rn
        FROM pdb_sequence ps
        JOIN sequence_type st ON st.stid = ps.stid
        JOIN blast_match bm   ON bm.query_sid = ps.sid
        JOIN chembl_hit ch ON ch.hit_sid = bm.hit_sid
        WHERE {coverage_sql}
    ),
    best AS (
        SELECT * FROM ranked WHERE rn = 1
    ),
    threshold_pairs AS (
        SELECT *
        FROM best
        WHERE percent_identity >= {threshold}
    )
    SELECT chembl_type, pdb_stid,
           COUNT(*) AS pair_count_raw,
           COUNT(DISTINCT pdb_sid) AS distinct_sequences
    FROM threshold_pairs
    GROUP BY chembl_type, pdb_stid;
    """


def build_legacy_interval_sql(min_floor: int, coverage_mode: str) -> str:
    coverage_sql = _legacy_coverage_predicate(coverage_mode)
    return f"""
    WITH chembl_hit AS (
        SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
        UNION ALL
        SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
    ),
    ranked AS (
        SELECT
            ps.sid AS pdb_sid,
            st.stid AS pdb_stid,
            ch.chembl_type,
            ch.hit_sid,
            bm.percent_identity,
            bm.query_coverage,
            bm.hit_coverage,
            ROW_NUMBER() OVER (
                PARTITION BY ps.sid, st.stid, ch.chembl_type, ch.hit_sid
                ORDER BY bm.percent_identity DESC, bm.query_coverage DESC, bm.hit_coverage DESC
            ) AS rn
        FROM pdb_sequence ps
        JOIN sequence_type st ON st.stid = ps.stid
        JOIN blast_match bm   ON bm.query_sid = ps.sid
        JOIN chembl_hit ch ON ch.hit_sid = bm.hit_sid
        WHERE {coverage_sql}
    ),
    best AS (
        SELECT * FROM ranked WHERE rn = 1 AND percent_identity >= {min_floor}
    ),
    binned AS (
        SELECT
            chembl_type,
            pdb_stid,
            pdb_sid,
            percent_identity,
            CASE
              WHEN percent_identity = 100 THEN '[100]'
              WHEN percent_identity >= 95 THEN '[95,100)'
              WHEN percent_identity >= 90 THEN '[90,95)'
              WHEN percent_identity >= 85 THEN '[85,90)'
              ELSE '[80,85)'
            END AS interval_label
        FROM best
    )
    SELECT interval_label,
           chembl_type,
           pdb_stid,
           COUNT(*) AS pair_count_raw,
           COUNT(DISTINCT pdb_sid) AS distinct_sequences
    FROM binned
    GROUP BY interval_label, chembl_type, pdb_stid;
    """


# --------------------------------------------------------------------------------------
# Unified (tier view) SQL builders
# --------------------------------------------------------------------------------------


def build_unified_right_bound_sql(threshold: int, coverage_mode: str) -> str:
    tier_view = COVERAGE_MODE_TO_TIER_VIEW[coverage_mode]
    return f"""
    WITH tier AS (
        SELECT
            ba.query_sid AS pdb_sid,
            ps.stid AS pdb_stid,
            ba.hit_sid
        FROM {tier_view} ba
        JOIN pdb_sequence ps ON ps.sid = ba.query_sid
    ),
    chembl_hit AS (
        SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
        UNION ALL
        SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
    ),
    joined AS (
        SELECT t.pdb_sid, t.pdb_stid, ch.chembl_type, ba.percent_identity
        FROM tier t
        JOIN {tier_view} ba ON ba.query_sid = t.pdb_sid AND ba.hit_sid = t.hit_sid
        JOIN chembl_hit ch ON ch.hit_sid = t.hit_sid
    ),
    threshold_pairs AS (
        SELECT * FROM joined WHERE percent_identity >= {threshold}
    )
    SELECT chembl_type, pdb_stid,
           COUNT(*) AS pair_count_raw,
           COUNT(DISTINCT pdb_sid) AS distinct_sequences
    FROM threshold_pairs
    GROUP BY chembl_type, pdb_stid;
    """


def build_unified_interval_sql(min_floor: int, coverage_mode: str) -> str:
    tier_view = COVERAGE_MODE_TO_TIER_VIEW[coverage_mode]
    return f"""
    WITH base AS (
        SELECT
            ba.query_sid AS pdb_sid,
            ba.hit_sid,
            ps.stid AS pdb_stid,
            ba.percent_identity
        FROM {tier_view} ba
        JOIN pdb_sequence ps ON ps.sid = ba.query_sid
        WHERE ba.percent_identity >= {min_floor}
    ),
    chembl_hit AS (
        SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
        UNION ALL
        SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
    ),
    joined AS (
        SELECT b.pdb_sid,
               b.pdb_stid,
               ch.chembl_type,
               b.percent_identity,
               CASE
                  WHEN b.percent_identity = 100 THEN '[100]'
                  WHEN b.percent_identity >= 95 THEN '[95,100)'
                  WHEN b.percent_identity >= 90 THEN '[90,95)'
                  WHEN b.percent_identity >= 85 THEN '[85,90)'
                  ELSE '[80,85)'
               END AS interval_label
        FROM base b
        JOIN chembl_hit ch ON ch.hit_sid = b.hit_sid
    )
    SELECT interval_label,
           chembl_type,
           pdb_stid,
           COUNT(*) AS pair_count_raw,
           COUNT(DISTINCT pdb_sid) AS distinct_sequences
    FROM joined
    GROUP BY interval_label, chembl_type, pdb_stid;
    """


# --------------------------------------------------------------------------------------
# Execution Helpers
# --------------------------------------------------------------------------------------


def fetch_as_dict(db: UnifiedPostgresConnector, sql: str) -> List[Dict]:
    rows = db.execute(sql) or []
    if not rows:
        return []
    # Attempt to fetch column names if cursor metadata available (UnifiedPostgresConnector expected to handle)
    # Fallback: map positional indexes generically.
    try:
        # Ugly but we assume internal connector returns list of tuples only.
        # Without metadata, caller should know column order.
        # We'll just return positional map inside calling code instead of here.
        pass
    except Exception as e:
        logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)
    # Just return list[tuple]; caller will interpret.
    return rows  # type: ignore


def to_metric_map_right(
    rows: List[Tuple], coverage_mode: str, threshold: int
) -> Dict[MetricKey, MetricRecord]:
    out: Dict[MetricKey, MetricRecord] = {}
    for chembl_type, pdb_stid, pair_count_raw, distinct_sequences in rows:
        mk = MetricKey(
            coverage_mode=coverage_mode,
            mode_type="right",
            threshold_or_interval=f">={threshold}",
            chembl_type=str(chembl_type),
            pdb_stid=int(pdb_stid),
        )
        out[mk] = MetricRecord(
            pair_count_raw=int(pair_count_raw or 0),
            distinct_sequences=int(distinct_sequences or 0),
        )
    return out


def to_metric_map_interval(
    rows: List[Tuple], coverage_mode: str
) -> Dict[MetricKey, MetricRecord]:
    out: Dict[MetricKey, MetricRecord] = {}
    for (
        interval_label,
        chembl_type,
        pdb_stid,
        pair_count_raw,
        distinct_sequences,
    ) in rows:
        mk = MetricKey(
            coverage_mode=coverage_mode,
            mode_type="interval",
            threshold_or_interval=str(interval_label),
            chembl_type=str(chembl_type),
            pdb_stid=int(pdb_stid),
        )
        out[mk] = MetricRecord(
            pair_count_raw=int(pair_count_raw or 0),
            distinct_sequences=int(distinct_sequences or 0),
        )
    return out


# --------------------------------------------------------------------------------------
# Parity Logic
# --------------------------------------------------------------------------------------


@dataclass
class Difference:
    key: MetricKey
    field: str
    legacy: int
    unified: int
    delta: int


def compare_metric_maps(
    legacy: Dict[MetricKey, MetricRecord], unified: Dict[MetricKey, MetricRecord]
) -> List[Difference]:
    diffs: List[Difference] = []
    all_keys = set(legacy.keys()) | set(unified.keys())
    for k in sorted(
        all_keys,
        key=lambda x: (
            x.coverage_mode,
            x.mode_type,
            x.threshold_or_interval,
            x.chembl_type,
            x.pdb_stid,
        ),
    ):
        l = legacy.get(k)
        u = unified.get(k)
        if l is None or u is None:
            diffs.append(
                Difference(
                    k,
                    "record_presence",
                    1 if l else 0,
                    1 if u else 0,
                    (1 if u else 0) - (1 if l else 0),
                )
            )
            continue
        if l.pair_count_raw != u.pair_count_raw:
            diffs.append(
                Difference(
                    k,
                    "pair_count_raw",
                    l.pair_count_raw,
                    u.pair_count_raw,
                    u.pair_count_raw - l.pair_count_raw,
                )
            )
        if l.distinct_sequences != u.distinct_sequences:
            diffs.append(
                Difference(
                    k,
                    "distinct_sequences",
                    l.distinct_sequences,
                    u.distinct_sequences,
                    u.distinct_sequences - l.distinct_sequences,
                )
            )
    return diffs


def format_differences(diffs: List[Difference]) -> str:
    if not diffs:
        return "No differences detected."
    lines = []
    lines.append("Differences Detected (legacy vs unified):")
    lines.append("-" * 110)
    lines.append(
        "coverage | mode | thr/interval | chembl | stid | field               | legacy | unified | delta"
    )
    lines.append("-" * 110)
    for d in diffs:
        lines.append(
            f"{d.key.coverage_mode:11} {d.key.mode_type:7} {d.key.threshold_or_interval:13} "
            f"{d.key.chembl_type:8} {d.key.pdb_stid:5} {d.field:19} "
            f"{d.legacy:7} {d.unified:8} {d.delta:+6}"
        )
    return "\n".join(lines)


# --------------------------------------------------------------------------------------
# Main Execution
# --------------------------------------------------------------------------------------


def run_parity_check(
    db: UnifiedPostgresConnector,
    coverage_modes: List[str],
    allow_differences: bool,
    min_identity_floor: int = 80,
) -> int:
    all_diffs: List[Difference] = []

    for coverage_mode in coverage_modes:
        mode_norm = coverage_mode.lower()
        if mode_norm not in VALID_COVERAGE_MODES:
            logger.error("Unsupported coverage mode '{}'", coverage_mode)
            return 2

        # Right-bound thresholds
        legacy_right_map: Dict[MetricKey, MetricRecord] = {}
        unified_right_map: Dict[MetricKey, MetricRecord] = {}
        for thr in RIGHT_THRESHOLDS:
            legacy_sql = build_legacy_right_bound_sql(thr, mode_norm)
            unified_sql = build_unified_right_bound_sql(thr, mode_norm)
            legacy_rows = db.execute(legacy_sql) or []
            unified_rows = db.execute(unified_sql) or []
            legacy_right_map.update(to_metric_map_right(legacy_rows, mode_norm, thr))
            unified_right_map.update(to_metric_map_right(unified_rows, mode_norm, thr))

        diffs_right = compare_metric_maps(legacy_right_map, unified_right_map)
        if diffs_right:
            logger.warning(
                "Right-bound differences for coverage_mode='{}' found: {}",
                coverage_mode,
                len(diffs_right),
            )
        all_diffs.extend(diffs_right)

        # Interval buckets
        legacy_interval_sql = build_legacy_interval_sql(min_identity_floor, mode_norm)
        unified_interval_sql = build_unified_interval_sql(min_identity_floor, mode_norm)
        legacy_interval_rows = db.execute(legacy_interval_sql) or []
        unified_interval_rows = db.execute(unified_interval_sql) or []
        legacy_interval_map = to_metric_map_interval(legacy_interval_rows, mode_norm)
        unified_interval_map = to_metric_map_interval(unified_interval_rows, mode_norm)
        diffs_interval = compare_metric_maps(legacy_interval_map, unified_interval_map)
        if diffs_interval:
            logger.warning(
                "Interval differences for coverage_mode='{}' found: {}",
                coverage_mode,
                len(diffs_interval),
            )
        all_diffs.extend(diffs_interval)

    # Reporting
    if all_diffs:
        diff_report = format_differences(all_diffs)
        logger.info("\n{}\n", diff_report)
        if not allow_differences:
            logger.error(
                "Parity check FAILED: differences detected (rerun with --allow-differences to ignore)."
            )
            return 1
        else:
            logger.info("Differences detected but allowed by flag.")
            return 0
    else:
        logger.success("Parity check PASSED: no differences detected.")
        return 0


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="sequence_mapping_parity_check",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__),
    )
    p.add_argument(
        "--coverage-modes",
        nargs="+",
        default=["identity-only", "no-hit", "full"],
        help=f"Coverage modes to test (subset of {VALID_COVERAGE_MODES}).",
    )
    p.add_argument("--host", help="DB host override")
    p.add_argument("--port", type=int, help="DB port override")
    p.add_argument("--user", help="DB user override")
    p.add_argument(
        "--password", help="DB password override (use env / config preferred)"
    )
    p.add_argument("--database-name", help="DB name override")
    p.add_argument(
        "--allow-differences",
        action="store_true",
        help="Do not exit non-zero on differences.",
    )
    p.add_argument(
        "--min-identity-floor",
        type=int,
        default=80,
        help="Minimum identity floor for interval bin inclusion.",
    )
    return p.parse_args(list(argv) if argv is not None else None)


def build_connector_from_args(ns: argparse.Namespace) -> UnifiedPostgresConnector:
    """
    Build connector using AppConfig-based configuration (CLI tunnel management).

    CLI overrides (host, port, user, password, db_name) are passed to AppConfig
    which merges them with appconfig.toml settings. SSH tunneling is managed at
    the CLI layer (--manage-tunnel, --tunnel-mode).
    """
    from activity_db_analysis.siteconfig.app_config import AppConfig

    cli_overrides = {
        k: v
        for k, v in {
            "host": ns.host,
            "port": ns.port,
            "user": ns.user,
            "password": ns.password,
            "name": ns.database_name,
        }.items()
        if v is not None
    }

    if cli_overrides:
        logger.info(
            "Building DB connector with CLI overrides: {}",
            {k: v for k, v in cli_overrides.items() if k != "password"},
        )

    # Use AppConfig which handles config file + CLI overrides + env var password
    return AppConfig.build_activity_sql_handler_from_config(cli_overrides=cli_overrides)


def main(argv: Optional[Iterable[str]] = None) -> int:
    ns = parse_args(argv)
    try:
        db = build_connector_from_args(ns)
    except Exception as e:
        logger.exception("Failed to establish database connection: {}", e)
        return 2

    try:
        rc = run_parity_check(
            db=db,
            coverage_modes=ns.coverage_modes,
            allow_differences=ns.allow_differences,
            min_identity_floor=ns.min_identity_floor,
        )
        return rc
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 2
    except Exception as e:
        logger.exception("Unexpected error during parity check: {}", e)
        return 2
    finally:
        try:
            db.close()
        except Exception as e:
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)


if __name__ == "__main__":
    sys.exit(main())
