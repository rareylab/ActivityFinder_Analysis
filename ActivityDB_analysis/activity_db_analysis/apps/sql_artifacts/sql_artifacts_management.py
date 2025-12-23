#!/usr/bin/env python3
"""
SQL Artifacts Management for Ligand Mapping Analysis
====================================================

A comprehensive tool for managing all SQL artifacts defined in the ligand mapping
DDL file, including functions, materialized views, views, and indexes.

This application provides fine-grained control over:
- Creating individual or groups of artifacts
- Dropping artifacts with dependency checking
- Refreshing materialized views (regular and concurrent)
- Status monitoring and reporting
- Batch operations with proper dependency ordering

All SQL definitions are embedded in Python for full control and type safety.
Database connections follow the same patterns as other apps in the repository,
including SSH tunnel support via AppConfig.

Usage Examples:
    # Individual operations
    python sql_artifacts_management.py create mv_ligand_sequence_target_fact
    python sql_artifacts_management.py refresh mv_ligand_sequence_target_fact --concurrent
    python sql_artifacts_management.py drop v_ligand_mapping_retention_summary --cascade

    # Batch operations
    python sql_artifacts_management.py create-all --type materialized_view
    python sql_artifacts_management.py status

    # SSH tunnels are managed by the CLI root - configure appconfig.toml [ssh] section
    uv run actfind sql-manager create-all
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple

import typer
from loguru import logger
from psycopg2 import sql

from activity_db_analysis.db.connector import (
    UnifiedPostgresConnector,
    build_unified_connector,
)
from activity_db_analysis.siteconfig.app_config import AppConfig

# ======================================================================================
# Data Models and Enums
# ======================================================================================


class ArtifactType(str, Enum):
    FUNCTION = "function"
    MATERIALIZED_VIEW = "materialized_view"
    VIEW = "view"
    INDEX = "index"


@dataclass
class SQLArtifact:
    name: str
    sql: str
    drop_sql: str
    dependencies: List[str]
    description: str
    artifact_type: ArtifactType

    # Optional specialized operations
    refresh_sql: Optional[str] = None
    refresh_concurrent_sql: Optional[str] = None
    status_sql: Optional[str] = None

    def can_refresh(self) -> bool:
        return self.artifact_type == ArtifactType.MATERIALIZED_VIEW

    def supports_concurrent_refresh(self) -> bool:
        return (
            self.artifact_type == ArtifactType.MATERIALIZED_VIEW
            and self.refresh_concurrent_sql is not None
        )


# ======================================================================================
# SQL Artifact Definitions (Embedded from DDL file)
# ======================================================================================

ARTIFACTS: Dict[str, SQLArtifact] = {
    # Functions
    "f_ligand_sequence_target_base": SQLArtifact(
        name="f_ligand_sequence_target_base",
        sql="""CREATE OR REPLACE FUNCTION f_ligand_sequence_target_base(
    identity_threshold NUMERIC DEFAULT 80,
    apply_query_coverage_filter BOOLEAN DEFAULT TRUE,
    query_coverage_threshold NUMERIC DEFAULT 80,
    apply_hit_coverage_filter BOOLEAN DEFAULT TRUE,
    hit_coverage_threshold NUMERIC DEFAULT 80
)
RETURNS TABLE (
    query_sid INTEGER,
    pdb_stid INTEGER,
    hit_sid INTEGER,
    chembl_type TEXT,
    chembl_tid INTEGER,
    pdblid INTEGER,
    clid INTEGER,
    smmclid INTEGER,
    percent_identity NUMERIC,
    query_coverage NUMERIC,
    hit_coverage NUMERIC,
    chain_id INTEGER,
    pdb_id INTEGER
) LANGUAGE sql STABLE PARALLEL SAFE AS
$$
WITH chembl_hit AS (
    SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type
    FROM chembl_component_sequence ccs
    UNION ALL
    SELECT cvs.sid AS hit_sid, 'VARIANT'::text
    FROM chembl_variant_sequence cvs
),
ranked AS (
    SELECT
        ps.sid AS query_sid,
        ps.cid AS chain_id,
        ps.pid AS pdb_id,
        st.stid AS pdb_stid,
        bm.hit_sid,
        ch.chembl_type,
        bm.percent_identity,
        bm.query_coverage,
        bm.hit_coverage,
        ROW_NUMBER() OVER (
            PARTITION BY ps.sid, bm.hit_sid, ch.chembl_type
            ORDER BY bm.percent_identity DESC,
                     bm.query_coverage DESC,
                     bm.hit_coverage DESC
        ) AS rn
    FROM pdb_sequence ps
    JOIN sequence_type st ON st.stid = ps.stid
    JOIN blast_match bm ON bm.query_sid = ps.sid
    JOIN chembl_hit ch ON ch.hit_sid = bm.hit_sid
),
best AS (
    SELECT
        r.query_sid,
        r.pdb_stid,
        r.hit_sid,
        r.chembl_type,
        r.percent_identity,
        r.query_coverage,
        r.hit_coverage,
        r.chain_id,
        r.pdb_id
    FROM ranked r
    WHERE r.rn = 1
      AND r.percent_identity IS NOT NULL
      AND (identity_threshold IS NULL OR r.percent_identity >= identity_threshold)
      AND r.query_coverage IS NOT NULL
      AND (
          NOT apply_query_coverage_filter
          OR (query_coverage_threshold IS NOT NULL AND r.query_coverage >= query_coverage_threshold)
      )
      AND r.hit_coverage IS NOT NULL
      AND (
          NOT apply_hit_coverage_filter
          OR (hit_coverage_threshold IS NOT NULL AND r.hit_coverage >= hit_coverage_threshold)
      )
),
seq_target AS (
    SELECT st.sid AS hit_sid, st.chembl_tid
    FROM sequence_target st
),
seq_to_ligand AS (
    SELECT plpi.pid, plpi.cid, plpi.pdblid
    FROM pdb_ligand_pdb_id plpi
),
fact AS (
    SELECT
        b.query_sid,
        b.pdb_stid,
        b.hit_sid,
        b.chembl_type,
        st.chembl_tid,
        sl.pdblid,
        plc.clid,
        plc.smmclid,
        b.percent_identity,
        b.query_coverage,
        b.hit_coverage,
        b.chain_id,
        b.pdb_id
    FROM best b
    JOIN seq_target st ON st.hit_sid = b.hit_sid
    JOIN pdb_sequence ps ON ps.sid = b.query_sid
    JOIN seq_to_ligand sl ON sl.pid = ps.pid AND sl.cid = ps.cid
    JOIN pdb_ligand_chembl_ligand plc ON plc.pdblid = sl.pdblid
)
SELECT * FROM fact;
$$;""",
        drop_sql="DROP FUNCTION IF EXISTS f_ligand_sequence_target_base CASCADE",
        dependencies=[],
        description="Generalized base function (identity/query/hit coverage filters) enforcing non-null metrics",
        artifact_type=ArtifactType.FUNCTION,
    ),
    # New unified-tier based base function leveraging mv_best_alignment_identity_floor and tier flags.
    # This removes internal ROW_NUMBER ranking & coverage threshold logic; instead it selects
    # pre-ranked best (query_sid, hit_sid) rows and filters by pass flags (query_cov_pass, hit_cov_pass).
    # Parameters:
    #   identity_threshold      : minimum percent identity (floor, default 80)
    #   require_query_coverage  : if TRUE enforce query_cov_pass (>=80)
    #   require_hit_coverage    : if TRUE enforce hit_cov_pass   (>=80)
    "f_ligand_sequence_target_base_v2": SQLArtifact(
        name="f_ligand_sequence_target_base_v2",
        sql="""CREATE OR REPLACE FUNCTION f_ligand_sequence_target_base_v2(
    identity_threshold NUMERIC DEFAULT 80,
    require_query_coverage BOOLEAN DEFAULT TRUE,
    require_hit_coverage BOOLEAN DEFAULT TRUE
)
RETURNS TABLE (
    query_sid INTEGER,
    pdb_stid INTEGER,
    hit_sid INTEGER,
    chembl_type TEXT,
    chembl_tid INTEGER,
    pdblid INTEGER,
    clid INTEGER,
    smmclid INTEGER,
    percent_identity NUMERIC,
    query_coverage NUMERIC,
    hit_coverage NUMERIC,
    chain_id INTEGER,
    pdb_id INTEGER
) LANGUAGE sql STABLE PARALLEL SAFE AS
$$
WITH base AS (
    SELECT
        ba.query_sid,
        st.stid AS pdb_stid,
        ba.hit_sid,
        ch.chembl_type,
        ba.percent_identity,
        ba.query_coverage,
        ba.hit_coverage,
        ps.cid AS chain_id,
        ps.pid AS pdb_id,
        ba.query_cov_pass,
        ba.hit_cov_pass
    FROM mv_best_alignment_identity_floor ba
    JOIN pdb_sequence ps ON ps.sid = ba.query_sid
    JOIN sequence_type st ON st.stid = ps.stid
    JOIN (
        SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
        UNION ALL
        SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
    ) ch ON ch.hit_sid = ba.hit_sid
    WHERE (identity_threshold IS NULL OR ba.percent_identity >= identity_threshold)
      AND (NOT require_query_coverage OR ba.query_cov_pass)
      AND (NOT require_hit_coverage OR ba.hit_cov_pass)
),
seq_target AS (
    SELECT st.sid AS hit_sid, st.chembl_tid
    FROM sequence_target st
),
seq_to_ligand AS (
    SELECT plpi.pid, plpi.cid, plpi.pdblid
    FROM pdb_ligand_pdb_id plpi
),
fact AS (
    SELECT
        b.query_sid,
        b.pdb_stid,
        b.hit_sid,
        b.chembl_type,
        st.chembl_tid,
        sl.pdblid,
        plc.clid,
        plc.smmclid,
        b.percent_identity,
        b.query_coverage,
        b.hit_coverage,
        b.chain_id,
        b.pdb_id
    FROM base b
    JOIN seq_target st ON st.hit_sid = b.hit_sid
    JOIN pdb_sequence ps ON ps.sid = b.query_sid
    JOIN seq_to_ligand sl ON sl.pid = ps.pid AND sl.cid = ps.cid
    JOIN pdb_ligand_chembl_ligand plc ON plc.pdblid = sl.pdblid
)
SELECT * FROM fact;
$$;""",
        drop_sql="DROP FUNCTION IF EXISTS f_ligand_sequence_target_base_v2(NUMERIC,BOOLEAN,BOOLEAN) CASCADE",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Unified-tier base function (v2) using pre-ranked best alignment MV + coverage pass flags.",
        artifact_type=ArtifactType.FUNCTION,
    ),
    "refresh_ligand_mapping_mvs": SQLArtifact(
        name="refresh_ligand_mapping_mvs",
        sql="""CREATE OR REPLACE FUNCTION refresh_ligand_mapping_mvs(
    refresh_concurrently BOOLEAN DEFAULT FALSE
)
RETURNS TABLE(
    mv_name TEXT,
    rows_before BIGINT,
    rows_after BIGINT,
    refresh_time INTERVAL
) AS $$
DECLARE
    start_time TIMESTAMP;
    end_time TIMESTAMP;
    strict_before BIGINT;
    strict_after BIGINT;
    relaxed_before BIGINT;
    relaxed_after BIGINT;
    identity_only_before BIGINT;
    identity_only_after BIGINT;
BEGIN
    -- Get initial row counts
    SELECT COUNT(*) INTO strict_before FROM mv_ligand_sequence_target_fact;
    SELECT COUNT(*) INTO relaxed_before FROM mv_ligand_sequence_target_fact_nohitcov;
    SELECT COUNT(*) INTO identity_only_before FROM mv_ligand_sequence_target_fact_identity_only;

    -- Refresh strict MV
    start_time := clock_timestamp();
    IF refresh_concurrently THEN
        REFRESH MATERIALIZED VIEW CONCURRENTLY mv_ligand_sequence_target_fact;
    ELSE
        REFRESH MATERIALIZED VIEW mv_ligand_sequence_target_fact;
    END IF;
    end_time := clock_timestamp();
    SELECT COUNT(*) INTO strict_after FROM mv_ligand_sequence_target_fact;

    -- Return strict result
    mv_name := 'mv_ligand_sequence_target_fact';
    rows_before := strict_before;
    rows_after := strict_after;
    refresh_time := end_time - start_time;
    RETURN NEXT;

    -- Refresh relaxed MV
    start_time := clock_timestamp();
    IF refresh_concurrently THEN
        REFRESH MATERIALIZED VIEW CONCURRENTLY mv_ligand_sequence_target_fact_nohitcov;
    ELSE
        REFRESH MATERIALIZED VIEW mv_ligand_sequence_target_fact_nohitcov;
    END IF;
    end_time := clock_timestamp();
    SELECT COUNT(*) INTO relaxed_after FROM mv_ligand_sequence_target_fact_nohitcov;

    -- Return relaxed result
    mv_name := 'mv_ligand_sequence_target_fact_nohitcov';
    rows_before := relaxed_before;
    rows_after := relaxed_after;
    refresh_time := end_time - start_time;
    RETURN NEXT;

    -- Refresh identity-only MV
    start_time := clock_timestamp();
    IF refresh_concurrently THEN
        REFRESH MATERIALIZED VIEW CONCURRENTLY mv_ligand_sequence_target_fact_identity_only;
    ELSE
        REFRESH MATERIALIZED VIEW mv_ligand_sequence_target_fact_identity_only;
    END IF;
    end_time := clock_timestamp();
    SELECT COUNT(*) INTO identity_only_after FROM mv_ligand_sequence_target_fact_identity_only;

    -- Return identity-only result
    mv_name := 'mv_ligand_sequence_target_fact_identity_only';
    rows_before := identity_only_before;
    rows_after := identity_only_after;
    refresh_time := end_time - start_time;
    RETURN NEXT;

    -- Update statistics
    ANALYZE mv_ligand_sequence_target_fact;
    ANALYZE mv_ligand_sequence_target_fact_nohitcov;
    ANALYZE mv_ligand_sequence_target_fact_identity_only;

    RETURN;
END;
$$ LANGUAGE plpgsql;""",
        drop_sql="DROP FUNCTION IF EXISTS refresh_ligand_mapping_mvs CASCADE",
        dependencies=[
            "mv_ligand_sequence_target_fact",
            "mv_ligand_sequence_target_fact_nohitcov",
            "mv_ligand_sequence_target_fact_identity_only",
        ],
        description="Refresh function for materialized views",
        artifact_type=ArtifactType.FUNCTION,
    ),
    "get_ligand_mapping_mv_status": SQLArtifact(
        name="get_ligand_mapping_mv_status",
        sql="""CREATE OR REPLACE FUNCTION get_ligand_mapping_mv_status()
RETURNS TABLE(
    matview_name TEXT,
    is_populated BOOLEAN,
    row_count BIGINT,
    distinct_quadruplets BIGINT,
    size_pretty TEXT,
    index_count INT,
    total_index_size TEXT
) AS $$
BEGIN
    RETURN QUERY
    WITH mv_stats AS (
        SELECT
            mv.matviewname::text AS name,
            mv.ispopulated,
            CASE
                WHEN mv.matviewname = 'mv_ligand_sequence_target_fact' THEN
                    (SELECT COUNT(*) FROM mv_ligand_sequence_target_fact)
                WHEN mv.matviewname = 'mv_ligand_sequence_target_fact_nohitcov' THEN
                    (SELECT COUNT(*) FROM mv_ligand_sequence_target_fact_nohitcov)
                WHEN mv.matviewname = 'mv_ligand_sequence_target_fact_identity_only' THEN
                    (SELECT COUNT(*) FROM mv_ligand_sequence_target_fact_identity_only)
                ELSE 0
            END AS rows,
            CASE
                WHEN mv.matviewname = 'mv_ligand_sequence_target_fact' THEN
                    (SELECT COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid))
                     FROM mv_ligand_sequence_target_fact)
                WHEN mv.matviewname = 'mv_ligand_sequence_target_fact_nohitcov' THEN
                    (SELECT COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid))
                     FROM mv_ligand_sequence_target_fact_nohitcov)
                WHEN mv.matviewname = 'mv_ligand_sequence_target_fact_identity_only' THEN
                    (SELECT COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid))
                     FROM mv_ligand_sequence_target_fact_identity_only)
                ELSE 0
            END AS quads,
            pg_size_pretty(pg_total_relation_size(mv.schemaname||'.'||mv.matviewname)) AS size,
            COUNT(i.indexname) AS idx_count,
            pg_size_pretty(COALESCE(SUM(pg_relation_size(i.schemaname||'.'||i.indexname)), 0)) AS idx_size
        FROM pg_matviews mv
        LEFT JOIN pg_indexes i ON i.tablename = mv.matviewname AND i.schemaname = mv.schemaname
        WHERE mv.schemaname = 'public'
            AND mv.matviewname LIKE '%ligand_sequence_target_fact%'
        GROUP BY mv.matviewname, mv.ispopulated, mv.schemaname
    )
    SELECT
        name,
        ispopulated,
        rows,
        quads,
        size,
        idx_count::int,
        idx_size
    FROM mv_stats
    ORDER BY name;
END;
$$ LANGUAGE plpgsql;""",
        drop_sql="DROP FUNCTION IF EXISTS get_ligand_mapping_mv_status CASCADE",
        dependencies=[],
        description="Status checking function for materialized views (includes identity-only tier)",
        artifact_type=ArtifactType.FUNCTION,
    ),
    "analyze_ligand_mapping_performance": SQLArtifact(
        name="analyze_ligand_mapping_performance",
        sql="""CREATE OR REPLACE FUNCTION analyze_ligand_mapping_performance()
RETURNS TABLE(
    check_name TEXT,
    status TEXT,
    details TEXT
) AS $$
BEGIN
    -- Check MV population (all three tiers)
    RETURN QUERY
    SELECT
        'MV Population'::text,
        CASE
            WHEN bool_and(ispopulated) THEN 'OK'::text
            ELSE 'ERROR'::text
        END,
        'All mapping MVs populated: ' || bool_and(ispopulated)::text
    FROM pg_matviews
    WHERE matviewname IN (
        'mv_ligand_sequence_target_fact',
        'mv_ligand_sequence_target_fact_nohitcov',
        'mv_ligand_sequence_target_fact_identity_only'
    );

    -- Per-MV row + quadruplet counts
    RETURN QUERY
    SELECT
        'MV Rows'::text AS check_name,
        mv_name::text   AS status,
        ('rows='||row_count||' quads='||quad_count)::text AS details
    FROM (
        SELECT 'mv_ligand_sequence_target_fact' AS mv_name,
               COUNT(*) AS row_count,
               COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid)) AS quad_count
          FROM mv_ligand_sequence_target_fact
        UNION ALL
        SELECT 'mv_ligand_sequence_target_fact_nohitcov',
               COUNT(*),
               COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid))
          FROM mv_ligand_sequence_target_fact_nohitcov
        UNION ALL
        SELECT 'mv_ligand_sequence_target_fact_identity_only',
               COUNT(*),
               COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid))
          FROM mv_ligand_sequence_target_fact_identity_only
    ) s;

    -- Check index usage (pattern still matches all three tiers)
    RETURN QUERY
    SELECT
        'Index Usage'::text,
        CASE
            WHEN AVG(idx_scan::numeric / NULLIF(seq_scan + idx_scan, 0)) > 0.9 THEN 'GOOD'::text
            WHEN AVG(idx_scan::numeric / NULLIF(seq_scan + idx_scan, 0)) > 0.5 THEN 'OK'::text
            ELSE 'POOR'::text
        END,
        'Average index usage: ' ||
        ROUND(AVG(idx_scan::numeric / NULLIF(seq_scan + idx_scan, 0)) * 100, 2)::text || '%'
    FROM pg_stat_user_tables
    WHERE schemaname = 'public'
        AND tablename LIKE '%ligand_sequence_target_fact%';

    -- Check table bloat
    RETURN QUERY
    SELECT
        'Table Bloat'::text,
        CASE
            WHEN AVG(n_dead_tup::numeric / NULLIF(n_live_tup + n_dead_tup, 0)) < 0.1 THEN 'OK'::text
            WHEN AVG(n_dead_tup::numeric / NULLIF(n_live_tup + n_dead_tup, 0)) < 0.2 THEN 'WARNING'::text
            ELSE 'HIGH'::text
        END,
        'Average dead tuple ratio: ' ||
        ROUND(AVG(n_dead_tup::numeric / NULLIF(n_live_tup + n_dead_tup, 0)) * 100, 2)::text || '%'
    FROM pg_stat_user_tables
    WHERE schemaname = 'public'
        AND tablename LIKE '%ligand_sequence_target_fact%';

    RETURN;
END;
$$ LANGUAGE plpgsql;""",
        drop_sql="DROP FUNCTION IF EXISTS analyze_ligand_mapping_performance CASCADE",
        dependencies=[],
        description="Performance analysis function (extended: includes per-MV row/quad counts)",
        artifact_type=ArtifactType.FUNCTION,
    ),
    # Materialized Views
    # ------------------------------------------------------------------
    # NEW FOUNDATIONAL ACTIVITY / ASSAY LIGAND MVs (Phase 1)
    # ------------------------------------------------------------------
    "mv_pdb_ligand_chembl_activity_assay_base": SQLArtifact(
        name="mv_pdb_ligand_chembl_activity_assay_base",
        sql="""CREATE MATERIALIZED VIEW mv_pdb_ligand_chembl_activity_assay_base AS
WITH mapped_pdb_ligand AS (
    SELECT
        plpi.pid,
        plpi.cid,
        plpi.pdblid,
        plc.clid,
        plc.smmclid
    FROM pdb_ligand_pdb_id plpi
    JOIN pdb_ligand_chembl_ligand plc ON plc.pdblid = plpi.pdblid
),
chembl_ligand_activity AS (
    SELECT
        cl.clid,
        cl.molregno,
        a.aid,
        a.asid,
        a.chembl_activity_id,
        a.chembl_assay_id,
        a.standard_type,
        a.standard_relation,
        a.standard_value,
        a.standard_units,
        a.pchembl_value,
        a.standard_flag,
        a.dvid,
        a.csid
    FROM chembl_ligand cl
    JOIN activity a ON a.molregno = cl.molregno
),
activity_assay_target AS (
    SELECT
        cla.clid,
        cla.molregno,
        cla.aid,
        cla.asid,
        cla.chembl_activity_id,
        cla.chembl_assay_id,
        cla.standard_type,
        cla.standard_relation,
        cla.standard_value,
        cla.standard_units,
        cla.pchembl_value,
        cla.standard_flag,
        cla.dvid,
        cla.csid,
        ass.tid,
        ass.variant_id,
        t.chembl_tid
    FROM chembl_ligand_activity cla
    JOIN assay ass ON ass.asid = cla.asid
    JOIN target t ON t.tid = ass.tid
)
SELECT
    p.pid,
    p.pdb_id              AS pdb_code,
    mpl.cid,
    mpl.pdblid,
    mpl.clid,
    mpl.smmclid,
    aat.molregno,
    aat.asid,
    aat.tid,
    aat.chembl_tid,
    aat.variant_id
FROM mapped_pdb_ligand mpl
JOIN activity_assay_target aat ON aat.clid = mpl.clid
JOIN pdb p ON p.pid = mpl.pid;
""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_pdb_ligand_chembl_activity_assay_base CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_pdb_ligand_chembl_activity_assay_base",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_pdb_ligand_chembl_activity_assay_base",
        dependencies=[],
        description="Full ligand→activity→assay convergence (no filters) collapsed to assay-level (no AID retention).",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_pdb_ligand_assay_link": SQLArtifact(
        name="mv_pdb_ligand_assay_link",
        sql="""CREATE MATERIALIZED VIEW mv_pdb_ligand_assay_link AS
SELECT DISTINCT
    pid,
    pdb_code,
    cid,
    pdblid,
    clid,
    smmclid,
    molregno,
    asid,
    tid,
    chembl_tid,
    variant_id
FROM mv_pdb_ligand_chembl_activity_assay_base;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_pdb_ligand_assay_link CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_pdb_ligand_assay_link",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_pdb_ligand_assay_link",
        dependencies=["mv_pdb_ligand_chembl_activity_assay_base"],
        description="Collapsed ligand↔assay link (distinct) used as convergence base for sequence mapping.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # NEW SEQUENCE PRESENCE FUNCTIONS (PARAMETERIZED)
    # ------------------------------------------------------------------
    "f_component_sequence_presence": SQLArtifact(
        name="f_component_sequence_presence",
        sql="""CREATE OR REPLACE FUNCTION f_component_sequence_presence(
    identity_threshold NUMERIC DEFAULT 80,
    apply_query_cov BOOLEAN DEFAULT FALSE,
    query_cov_threshold NUMERIC DEFAULT 0.80,
    apply_hit_cov BOOLEAN DEFAULT FALSE,
    hit_cov_threshold NUMERIC DEFAULT 0.80
) RETURNS TABLE(
    pid INTEGER,
    pdb_code TEXT,
    chain_id INTEGER,
    query_sid INTEGER,
    hit_sid INTEGER,
    chembl_tid INTEGER,
    tid INTEGER,
    variant_id BIGINT,
    pdblid INTEGER,
    clid INTEGER,
    smmclid INTEGER,
    percent_identity NUMERIC,
    query_coverage NUMERIC,
    hit_coverage NUMERIC,
    route_type TEXT
) LANGUAGE sql STABLE PARALLEL SAFE AS
$$
-- FIX (2025-01-14): Re-implemented presence logic as a 3-step pipeline:
--   1. best_blast_matches: DISTINCT ON (query_sid, hit_sid) picks the top BLAST
--      alignment per pair using deterministic ordering (identity, coverage, bmid).
--   2. sequence_mappings: attach all PDB (pid,cid) sequences whose sid matched.
--   3. Final SELECT: join ALL ligand mappings (no ROW_NUMBER ligand pruning).
-- Previous implementation applied ROW_NUMBER() over (pid,cid,stid,hit_sid) after
-- joining ligands and filtered rn=1, discarding valid ligand combinations.
WITH best_blast_matches AS (
    SELECT DISTINCT ON (query_sid, hit_sid)
        query_sid,
        hit_sid,
        percent_identity,
        query_coverage,
        hit_coverage,
        bmid
    FROM blast_match
    WHERE percent_identity >= identity_threshold
      AND (NOT apply_query_cov OR query_cov_threshold IS NULL OR query_coverage IS NULL OR query_coverage >= query_cov_threshold)
      AND (NOT apply_hit_cov OR hit_cov_threshold IS NULL OR hit_coverage IS NULL OR hit_coverage >= hit_cov_threshold)
    ORDER BY query_sid,
             hit_sid,
             percent_identity DESC,
             query_coverage DESC NULLS LAST,
             hit_coverage DESC NULLS LAST,
             bmid ASC
),
sequence_mappings AS (
    SELECT
        ps.pid,
        p.pdb_id AS pdb_code,
        ps.cid AS chain_id,
        ps.sid AS query_sid,
        bbm.hit_sid,
        bbm.percent_identity,
        bbm.query_coverage,
        bbm.hit_coverage,
        st.chembl_tid
    FROM pdb_sequence ps
    JOIN pdb p ON p.pid = ps.pid
    JOIN best_blast_matches bbm ON bbm.query_sid = ps.sid
    JOIN chembl_component_sequence ccs ON ccs.sid = bbm.hit_sid
    JOIN sequence_target st ON st.sid = ccs.sid
)
SELECT
    sm.pid,
    sm.pdb_code,
    sm.chain_id,
    sm.query_sid,
    sm.hit_sid,
    sm.chembl_tid,
    la.tid,
    la.variant_id,
    la.pdblid,
    la.clid,
    la.smmclid,
    sm.percent_identity,
    sm.query_coverage,
    sm.hit_coverage,
    'COMPONENT'::text AS route_type
FROM sequence_mappings sm
JOIN mv_pdb_ligand_assay_link la
  ON la.pid = sm.pid
 AND la.cid = sm.chain_id
 AND la.chembl_tid = sm.chembl_tid
 AND la.variant_id IS NULL;
$$;""",
        drop_sql="DROP FUNCTION IF EXISTS f_component_sequence_presence(NUMERIC,BOOLEAN,NUMERIC,BOOLEAN,NUMERIC) CASCADE",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Component route presence (best BLAST per (query_sid,hit_sid); ALL ligands retained).",
        artifact_type=ArtifactType.FUNCTION,
    ),
    "f_variant_sequence_presence": SQLArtifact(
        name="f_variant_sequence_presence",
        sql="""CREATE OR REPLACE FUNCTION f_variant_sequence_presence(
    identity_threshold NUMERIC DEFAULT 80,
    apply_query_cov BOOLEAN DEFAULT FALSE,
    query_cov_threshold NUMERIC DEFAULT 0.80,
    apply_hit_cov BOOLEAN DEFAULT FALSE,
    hit_cov_threshold NUMERIC DEFAULT 0.80
) RETURNS TABLE(
    pid INTEGER,
    pdb_code TEXT,
    chain_id INTEGER,
    query_sid INTEGER,
    hit_sid INTEGER,
    chembl_tid INTEGER,
    tid INTEGER,
    variant_id BIGINT,
    pdblid INTEGER,
    clid INTEGER,
    smmclid INTEGER,
    percent_identity NUMERIC,
    query_coverage NUMERIC,
    hit_coverage NUMERIC,
    route_type TEXT
) LANGUAGE sql STABLE PARALLEL SAFE AS
$$
-- FIX (2025-01-14): Align variant route with component route 3-step logic:
--   1. Best BLAST per (query_sid, hit_sid) via DISTINCT ON
--   2. Map to all (pid,cid) sequences referencing that query_sid
--   3. Join ALL ligand mappings constrained by variant_id; no ligand pruning.
WITH best_blast_matches AS (
    SELECT DISTINCT ON (query_sid, hit_sid)
        query_sid,
        hit_sid,
        percent_identity,
        query_coverage,
        hit_coverage,
        bmid
    FROM blast_match
    WHERE percent_identity >= identity_threshold
      AND (NOT apply_query_cov OR query_cov_threshold IS NULL OR query_coverage IS NULL OR query_coverage >= query_cov_threshold)
      AND (NOT apply_hit_cov OR hit_cov_threshold IS NULL OR hit_coverage IS NULL OR hit_coverage >= hit_cov_threshold)
    ORDER BY query_sid,
             hit_sid,
             percent_identity DESC,
             query_coverage DESC NULLS LAST,
             hit_coverage DESC NULLS LAST,
             bmid ASC
),
sequence_mappings AS (
    SELECT
        ps.pid,
        p.pdb_id AS pdb_code,
        ps.cid AS chain_id,
        ps.sid AS query_sid,
        bbm.hit_sid,
        bbm.percent_identity,
        bbm.query_coverage,
        bbm.hit_coverage,
        cvs.variant_id
    FROM pdb_sequence ps
    JOIN pdb p ON p.pid = ps.pid
    JOIN best_blast_matches bbm ON bbm.query_sid = ps.sid
    JOIN chembl_variant_sequence cvs ON cvs.sid = bbm.hit_sid
)
SELECT
    sm.pid,
    sm.pdb_code,
    sm.chain_id,
    sm.query_sid,
    sm.hit_sid,
    la.chembl_tid,
    la.tid,
    la.variant_id,
    la.pdblid,
    la.clid,
    la.smmclid,
    sm.percent_identity,
    sm.query_coverage,
    sm.hit_coverage,
    'VARIANT'::text AS route_type
FROM sequence_mappings sm
JOIN mv_pdb_ligand_assay_link la
  ON la.pid = sm.pid
 AND la.cid = sm.chain_id
 AND la.variant_id = sm.variant_id;
$$;""",
        drop_sql="DROP FUNCTION IF EXISTS f_variant_sequence_presence(NUMERIC,BOOLEAN,NUMERIC,BOOLEAN,NUMERIC) CASCADE",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Variant route presence (best BLAST per (query_sid,hit_sid); ALL ligands for matching variant).",
        artifact_type=ArtifactType.FUNCTION,
    ),
    # ------------------------------------------------------------------
    # PRESENCE MATERIALIZED VIEWS (IDENTITY-ONLY & COVERAGE)
    # ------------------------------------------------------------------
    "mv_component_sequence_presence_identity_only": SQLArtifact(
        name="mv_component_sequence_presence_identity_only",
        sql="CREATE MATERIALIZED VIEW mv_component_sequence_presence_identity_only AS SELECT * FROM f_component_sequence_presence(80, FALSE, 0.80, FALSE, 0.80);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_component_sequence_presence_identity_only CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_component_sequence_presence_identity_only",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_sequence_presence_identity_only",
        dependencies=["f_component_sequence_presence"],
        description="Component route presence (identity ≥80 only).",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_variant_sequence_presence_identity_only": SQLArtifact(
        name="mv_variant_sequence_presence_identity_only",
        sql="CREATE MATERIALIZED VIEW mv_variant_sequence_presence_identity_only AS SELECT * FROM f_variant_sequence_presence(80, FALSE, 0.80, FALSE, 0.80);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_variant_sequence_presence_identity_only CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_variant_sequence_presence_identity_only",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_variant_sequence_presence_identity_only",
        dependencies=["f_variant_sequence_presence"],
        description="Variant route presence (identity ≥80 only).",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_component_sequence_presence_cov": SQLArtifact(
        name="mv_component_sequence_presence_cov",
        sql="CREATE MATERIALIZED VIEW mv_component_sequence_presence_cov AS SELECT * FROM f_component_sequence_presence(80, TRUE, 0.80, TRUE, 0.80);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_component_sequence_presence_cov CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_component_sequence_presence_cov",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_sequence_presence_cov",
        dependencies=["f_component_sequence_presence"],
        description="Component route presence (identity + query & hit coverage ≥80%).",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_variant_sequence_presence_cov": SQLArtifact(
        name="mv_variant_sequence_presence_cov",
        sql="CREATE MATERIALIZED VIEW mv_variant_sequence_presence_cov AS SELECT * FROM f_variant_sequence_presence(80, TRUE, 0.80, TRUE, 0.80);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_variant_sequence_presence_cov CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_variant_sequence_presence_cov",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_variant_sequence_presence_cov",
        dependencies=["f_variant_sequence_presence"],
        description="Variant route presence (identity + query & hit coverage ≥80%).",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # OVERLAP & UNION COUNTS
    # ------------------------------------------------------------------
    "mv_component_variant_presence_overlap": SQLArtifact(
        name="mv_component_variant_presence_overlap",
        sql="""CREATE MATERIALIZED VIEW mv_component_variant_presence_overlap AS
SELECT
    c.pdb_code,
    c.chembl_tid,
    c.clid,
    c.pdblid,
    c.percent_identity AS component_identity,
    v.percent_identity AS variant_identity
FROM mv_component_sequence_presence_identity_only c
JOIN mv_variant_sequence_presence_identity_only v
  ON v.pdb_code = c.pdb_code
 AND v.chembl_tid = c.chembl_tid
 AND v.clid = c.clid
 AND v.pdblid = c.pdblid;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_component_variant_presence_overlap CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_component_variant_presence_overlap",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_variant_presence_overlap",
        dependencies=[
            "mv_component_sequence_presence_identity_only",
            "mv_variant_sequence_presence_identity_only",
        ],
        description="Overlap set of presence tuples reachable via both component and variant paths.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_component_presence_counts": SQLArtifact(
        name="mv_component_presence_counts",
        sql="""CREATE MATERIALIZED VIEW mv_component_presence_counts AS
SELECT
    COUNT(DISTINCT pdb_code)   AS pdbs,
    COUNT(DISTINCT chembl_tid) AS targets,
    COUNT(DISTINCT clid)       AS chembl_ligands,
    COUNT(DISTINCT pdblid)     AS pdb_ligands,
    COUNT(*)                   AS presence_rows,
    now()                      AS generated_at
FROM mv_component_sequence_presence_identity_only;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_component_presence_counts CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_component_presence_counts",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_presence_counts",
        dependencies=["mv_component_sequence_presence_identity_only"],
        description="Counts for component presence universe.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_variant_presence_counts": SQLArtifact(
        name="mv_variant_presence_counts",
        sql="""CREATE MATERIALIZED VIEW mv_variant_presence_counts AS
SELECT
    COUNT(DISTINCT pdb_code)   AS pdbs,
    COUNT(DISTINCT chembl_tid) AS targets,
    COUNT(DISTINCT clid)       AS chembl_ligands,
    COUNT(DISTINCT pdblid)     AS pdb_ligands,
    COUNT(*)                   AS presence_rows,
    now()                      AS generated_at
FROM mv_variant_sequence_presence_identity_only;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_variant_presence_counts CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_variant_presence_counts",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_variant_presence_counts",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Counts for variant presence universe.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_presence_union_counts": SQLArtifact(
        name="mv_presence_union_counts",
        sql="""CREATE MATERIALIZED VIEW mv_presence_union_counts AS
WITH unified AS (
    SELECT pdb_code, chembl_tid, clid, pdblid FROM mv_component_sequence_presence_identity_only
    UNION
    SELECT pdb_code, chembl_tid, clid, pdblid FROM mv_variant_sequence_presence_identity_only
)
SELECT
    COUNT(DISTINCT pdb_code)   AS pdbs,
    COUNT(DISTINCT chembl_tid) AS targets,
    COUNT(DISTINCT clid)       AS chembl_ligands,
    COUNT(DISTINCT pdblid)     AS pdb_ligands,
    COUNT(*)                   AS presence_rows,
    now()                      AS generated_at
FROM unified;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_presence_union_counts CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_presence_union_counts",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_union_counts",
        dependencies=[
            "mv_component_sequence_presence_identity_only",
            "mv_variant_sequence_presence_identity_only",
        ],
        description="Union counts across component and variant presence universes.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_presence_overlap_counts": SQLArtifact(
        name="mv_presence_overlap_counts",
        sql="""CREATE MATERIALIZED VIEW mv_presence_overlap_counts AS
SELECT
    COUNT(DISTINCT pdb_code)   AS pdbs,
    COUNT(DISTINCT chembl_tid) AS targets,
    COUNT(DISTINCT clid)       AS chembl_ligands,
    COUNT(DISTINCT pdblid)     AS pdb_ligands,
    COUNT(*)                   AS overlap_rows,
    now()                      AS generated_at
FROM mv_component_variant_presence_overlap;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_presence_overlap_counts CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_presence_overlap_counts",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_overlap_counts",
        dependencies=["mv_component_variant_presence_overlap"],
        description="Overlap counts between component and variant presence universes.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # QUADRUPLET MATRIX & COUNTS (NEW)
    # ------------------------------------------------------------------
    "mv_presence_quadruplet_matrix": SQLArtifact(
        name="mv_presence_quadruplet_matrix",
        sql="""CREATE MATERIALIZED VIEW mv_presence_quadruplet_matrix AS
WITH component AS (
    SELECT DISTINCT pdb_code, pdblid, chembl_tid, clid
    FROM mv_component_sequence_presence_identity_only
),
variant AS (
    SELECT DISTINCT pdb_code, pdblid, chembl_tid, clid
    FROM mv_variant_sequence_presence_identity_only
),
full_union AS (
    SELECT
        COALESCE(c.pdb_code, v.pdb_code)     AS pdb_code,
        COALESCE(c.pdblid,   v.pdblid)       AS pdblid,
        COALESCE(c.chembl_tid, v.chembl_tid) AS chembl_tid,
        COALESCE(c.clid,     v.clid)         AS clid,
        (c.pdb_code IS NOT NULL) AS in_component,
        (v.pdb_code IS NOT NULL) AS in_variant
    FROM component c
    FULL OUTER JOIN variant v
      USING (pdb_code, pdblid, chembl_tid, clid)
)
SELECT
    pdb_code,
    pdblid,
    chembl_tid,
    clid,
    in_component,
    in_variant,
    CASE
        WHEN in_component AND in_variant THEN 'both'
        WHEN in_component THEN 'component_only'
        WHEN in_variant THEN 'variant_only'
        ELSE 'error'
    END AS membership
FROM full_union;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_presence_quadruplet_matrix CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_presence_quadruplet_matrix",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_quadruplet_matrix",
        dependencies=[
            "mv_component_sequence_presence_identity_only",
            "mv_variant_sequence_presence_identity_only",
        ],
        description="Distinct quadruplets (pdb_code,pdblid,chembl_tid,clid) with component/variant membership flags.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_presence_quadruplet_counts": SQLArtifact(
        name="mv_presence_quadruplet_counts",
        sql="""CREATE MATERIALIZED VIEW mv_presence_quadruplet_counts AS
WITH base AS (
    SELECT * FROM mv_presence_quadruplet_matrix
),
component AS (
    SELECT
        COUNT(*)                       AS quadruplets,
        COUNT(DISTINCT pdb_code)       AS pdbs,
        COUNT(DISTINCT pdblid)         AS pdb_ligands,
        COUNT(DISTINCT clid)           AS chembl_ligands,
        COUNT(DISTINCT chembl_tid)     AS targets
    FROM base
    WHERE in_component
),
variant AS (
    SELECT
        COUNT(*)                       AS quadruplets,
        COUNT(DISTINCT pdb_code)       AS pdbs,
        COUNT(DISTINCT pdblid)         AS pdb_ligands,
        COUNT(DISTINCT clid)           AS chembl_ligands,
        COUNT(DISTINCT chembl_tid)     AS targets
    FROM base
    WHERE in_variant
),
union_all AS (
    SELECT
        COUNT(*)                       AS quadruplets,
        COUNT(DISTINCT pdb_code)       AS pdbs,
        COUNT(DISTINCT pdblid)         AS pdb_ligands,
        COUNT(DISTINCT clid)           AS chembl_ligands,
        COUNT(DISTINCT chembl_tid)     AS targets
    FROM base
),
overlap AS (
    SELECT
        COUNT(*)                       AS quadruplets,
        COUNT(DISTINCT pdb_code)       AS pdbs,
        COUNT(DISTINCT pdblid)         AS pdb_ligands,
        COUNT(DISTINCT clid)           AS chembl_ligands,
        COUNT(DISTINCT chembl_tid)     AS targets
    FROM base
    WHERE in_component AND in_variant
)
SELECT 'component' AS route, * , now() AS generated_at FROM component
UNION ALL
SELECT 'variant',  * , now() FROM variant
UNION ALL
SELECT 'union',    * , now() FROM union_all
UNION ALL
SELECT 'overlap',  * , now() FROM overlap;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_presence_quadruplet_counts CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_presence_quadruplet_counts",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_quadruplet_counts",
        dependencies=["mv_presence_quadruplet_matrix"],
        description="Pre-aggregated quadruplet, PDB, ligand, and target counts for component/variant/union/overlap.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # QUADRUPLET MATRIX INDEXES (NEW)
    # ------------------------------------------------------------------
    "mv_presence_quadruplet_matrix_quad_idx": SQLArtifact(
        name="mv_presence_quadruplet_matrix_quad_idx",
        sql="CREATE INDEX mv_presence_quadruplet_matrix_quad_idx ON mv_presence_quadruplet_matrix (pdb_code, pdblid, chembl_tid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_presence_quadruplet_matrix_quad_idx",
        dependencies=["mv_presence_quadruplet_matrix"],
        description="Composite quadruplet index to accelerate membership-filtered distinct scans.",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_presence_quadruplet_matrix_membership_idx": SQLArtifact(
        name="mv_presence_quadruplet_matrix_membership_idx",
        sql="CREATE INDEX mv_presence_quadruplet_matrix_membership_idx ON mv_presence_quadruplet_matrix (in_component, in_variant);",
        drop_sql="DROP INDEX IF EXISTS mv_presence_quadruplet_matrix_membership_idx",
        dependencies=["mv_presence_quadruplet_matrix"],
        description="Boolean membership flag index to speed component/variant/overlap filtering.",
        artifact_type=ArtifactType.INDEX,
    ),
    # ------------------------------------------------------------------
    # UNION VIEW (NON-MV) FOR ANALYSIS
    # ------------------------------------------------------------------
    "v_component_variant_presence_union": SQLArtifact(
        name="v_component_variant_presence_union",
        sql="""CREATE OR REPLACE VIEW v_component_variant_presence_union AS
SELECT *
FROM mv_component_sequence_presence_identity_only
UNION ALL
SELECT *
FROM mv_variant_sequence_presence_identity_only;""",
        drop_sql="DROP VIEW IF EXISTS v_component_variant_presence_union CASCADE",
        dependencies=[
            "mv_component_sequence_presence_identity_only",
            "mv_variant_sequence_presence_identity_only",
        ],
        description="Union (ALL) of component and variant presence rows with route_type flag preserved.",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_presence_quadruplet_enriched": SQLArtifact(
        name="v_presence_quadruplet_enriched",
        sql="""CREATE OR REPLACE VIEW v_presence_quadruplet_enriched AS
SELECT
    m.pdb_code,
    m.pdblid,
    pl.ligand_id      AS pdb_ligand_id,
    m.clid,
    cl.chembl_id      AS chembl_ligand_id,
    m.chembl_tid,
    t.chembl_id       AS chembl_target_id,
    m.in_component,
    m.in_variant,
    m.membership
FROM mv_presence_quadruplet_matrix m
JOIN pdb_ligand pl      ON pl.pdblid = m.pdblid
JOIN chembl_ligand cl   ON cl.clid   = m.clid
JOIN target t           ON t.tid     = m.chembl_tid;""",
        drop_sql="DROP VIEW IF EXISTS v_presence_quadruplet_enriched CASCADE",
        dependencies=[
            "mv_presence_quadruplet_matrix",
        ],
        description="Human-readable quadruplet view (adds textual PDB ligand, ChEMBL ligand, and target identifiers).",
        artifact_type=ArtifactType.VIEW,
    ),
    # (Existing materialized view entries follow)
    "mv_ligand_sequence_target_fact": SQLArtifact(
        name="mv_ligand_sequence_target_fact",
        sql="CREATE MATERIALIZED VIEW mv_ligand_sequence_target_fact AS SELECT * FROM f_ligand_sequence_target_base(80, TRUE, 80, TRUE, 80);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_sequence_target_fact CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_sequence_target_fact",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_ligand_sequence_target_fact",
        dependencies=["f_ligand_sequence_target_base"],
        description="Strict materialized view (identity + query + hit coverage) with complete non-null metrics",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_ligand_sequence_target_fact_nohitcov": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov",
        sql="CREATE MATERIALIZED VIEW mv_ligand_sequence_target_fact_nohitcov AS SELECT * FROM f_ligand_sequence_target_base(80, TRUE, 80, FALSE, NULL);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_sequence_target_fact_nohitcov CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_sequence_target_fact_nohitcov",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_ligand_sequence_target_fact_nohitcov",
        dependencies=["f_ligand_sequence_target_base"],
        description="Identity + query coverage filtered (no hit coverage filter) with complete non-null metrics",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_active_pair_ligands": SQLArtifact(
        name="mv_active_pair_ligands",
        sql="""CREATE MATERIALIZED VIEW mv_active_pair_ligands AS
SELECT DISTINCT
    plc.clid,
    cl.molregno
FROM pdb_ligand_chembl_ligand plc
JOIN chembl_ligand cl ON cl.clid = plc.clid;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_active_pair_ligands CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_active_pair_ligands",
        dependencies=[],
        description="Active pair ligands for activity qualification",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_ligand_activity_target": SQLArtifact(
        name="mv_ligand_activity_target",
        sql="""CREATE MATERIALIZED VIEW mv_ligand_activity_target AS
SELECT DISTINCT
    ap.clid,
    a.aid,
    ay.asid,
    t.tid,
    t.chembl_tid,
    a.molregno
FROM mv_active_pair_ligands ap
JOIN activity a ON a.molregno = ap.molregno
JOIN assay ay   ON ay.asid = a.asid
JOIN target t   ON t.tid = ay.tid
WHERE t.chembl_tid IS NOT NULL;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_activity_target CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_activity_target",
        dependencies=["mv_active_pair_ligands"],
        description="Ligand to activity to target mapping",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # Unified Best Alignment Layer (single MV + tiered views)
    # Provides a canonical best (query_sid, hit_sid) alignment at identity >=80
    # with non-null coverage metrics and boolean pass flags. Downstream layers
    # (sequence mapping reports, ligand mapping presence / retention) can
    # reference the tier views instead of duplicating ROW_NUMBER logic.
    # Tiers (mnemonic):
    #   v_best_alignment_id            : identity floor only
    #   v_best_alignment_id_query      : identity + query coverage >=80
    #   v_best_alignment_id_query_hit  : identity + query coverage >=80 + hit coverage >=80
    # ------------------------------------------------------------------
    "mv_best_alignment_identity_floor": SQLArtifact(
        name="mv_best_alignment_identity_floor",
        sql="""CREATE MATERIALIZED VIEW mv_best_alignment_identity_floor AS
WITH ranked AS (
    SELECT
        bm.query_sid,
        bm.hit_sid,
        bm.percent_identity,
        bm.query_coverage,
        bm.hit_coverage,
        ROW_NUMBER() OVER (
            PARTITION BY bm.query_sid, bm.hit_sid
            ORDER BY
                bm.percent_identity DESC,
                bm.query_coverage   DESC,
                bm.hit_coverage     DESC,
                bm.bmid             ASC
        ) AS rn
    FROM blast_match bm
    WHERE bm.percent_identity >= 80
      AND bm.query_coverage IS NOT NULL
      AND bm.hit_coverage   IS NOT NULL
)
SELECT
    query_sid,
    hit_sid,
    percent_identity,
    query_coverage,
    hit_coverage,
    (query_coverage >= 80.0) AS query_cov_pass,
    (hit_coverage   >= 80.0) AS hit_cov_pass
FROM ranked
WHERE rn = 1;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_best_alignment_identity_floor CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_best_alignment_identity_floor",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_best_alignment_identity_floor",
        dependencies=[],
        description="Canonical best alignment (identity >=80, non-null coverages) with pass flags.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "v_best_alignment_id": SQLArtifact(
        name="v_best_alignment_id",
        sql="""CREATE OR REPLACE VIEW v_best_alignment_id AS
SELECT * FROM mv_best_alignment_identity_floor;""",
        drop_sql="DROP VIEW IF EXISTS v_best_alignment_id CASCADE",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Tier view: identity floor only (no coverage thresholds).",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_best_alignment_id_query": SQLArtifact(
        name="v_best_alignment_id_query",
        sql="""CREATE OR REPLACE VIEW v_best_alignment_id_query AS
SELECT *
FROM mv_best_alignment_identity_floor
WHERE query_cov_pass;""",
        drop_sql="DROP VIEW IF EXISTS v_best_alignment_id_query CASCADE",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Tier view: identity floor + query coverage >=80%.",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_best_alignment_id_query_hit": SQLArtifact(
        name="v_best_alignment_id_query_hit",
        sql="""CREATE OR REPLACE VIEW v_best_alignment_id_query_hit AS
SELECT *
FROM mv_best_alignment_identity_floor
WHERE query_cov_pass
  AND hit_cov_pass;""",
        drop_sql="DROP VIEW IF EXISTS v_best_alignment_id_query_hit CASCADE",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Tier view: identity floor + query coverage >=80% + hit coverage >=80%.",
        artifact_type=ArtifactType.VIEW,
    ),
    # Indexes for unified best alignment MV (enable CONCURRENT refresh + filters)
    "mv_best_alignment_identity_floor_unique_idx": SQLArtifact(
        name="mv_best_alignment_identity_floor_unique_idx",
        sql="CREATE UNIQUE INDEX mv_best_alignment_identity_floor_unique_idx ON mv_best_alignment_identity_floor (query_sid, hit_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity_floor_unique_idx",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Unique index (required for CONCURRENT refresh, ensures 1 row per pair).",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_best_alignment_identity_floor_query_sid_idx": SQLArtifact(
        name="mv_best_alignment_identity_floor_query_sid_idx",
        sql="CREATE INDEX mv_best_alignment_identity_floor_query_sid_idx ON mv_best_alignment_identity_floor (query_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity_floor_query_sid_idx",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Lookup index on query_sid.",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_best_alignment_identity_floor_hit_sid_idx": SQLArtifact(
        name="mv_best_alignment_identity_floor_hit_sid_idx",
        sql="CREATE INDEX mv_best_alignment_identity_floor_hit_sid_idx ON mv_best_alignment_identity_floor (hit_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity_floor_hit_sid_idx",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Lookup index on hit_sid.",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_best_alignment_identity_floor_identity_desc_idx": SQLArtifact(
        name="mv_best_alignment_identity_floor_identity_desc_idx",
        sql="CREATE INDEX mv_best_alignment_identity_floor_identity_desc_idx ON mv_best_alignment_identity_floor (percent_identity DESC);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity_floor_identity_desc_idx",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Ordering / range filter index on percent_identity (DESC).",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_best_alignment_identity_floor_cov_flags_idx": SQLArtifact(
        name="mv_best_alignment_identity_floor_cov_flags_idx",
        sql="CREATE INDEX mv_best_alignment_identity_floor_cov_flags_idx ON mv_best_alignment_identity_floor (query_cov_pass, hit_cov_pass);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity_floor_cov_flags_idx",
        dependencies=["mv_best_alignment_identity_floor"],
        description="Boolean coverage pass flag index to accelerate tier view filtering.",
        artifact_type=ArtifactType.INDEX,
    ),
    # ------------------------------------------------------------------
    # NEW: Chain-level totals MV (precomputes total & mapped chain counts per tier and sequence type)
    # ------------------------------------------------------------------
    "mv_sequence_totals_chain": SQLArtifact(
        name="mv_sequence_totals_chain",
        sql="""CREATE MATERIALIZED VIEW mv_sequence_totals_chain AS
WITH intersection_chains AS (
    SELECT pid, cid
    FROM (
        SELECT ps.pid, ps.cid, ps.stid
        FROM pdb_sequence ps
        WHERE ps.stid IN (1,4,7)
        GROUP BY ps.pid, ps.cid, ps.stid
    ) s
    GROUP BY pid, cid
    HAVING COUNT(DISTINCT stid) = 3
),
all_chains AS (
    SELECT ps.pid, ps.cid, st.stid, st.type AS pdb_sequence_type
    FROM pdb_sequence ps
    JOIN sequence_type st ON st.stid = ps.stid
    JOIN intersection_chains ic ON ic.pid = ps.pid AND ic.cid = ps.cid
),
tier_union AS (
    SELECT 'identity-only'::text AS tier, query_sid FROM v_best_alignment_id
    UNION ALL
    SELECT 'no-hit'::text, query_sid FROM v_best_alignment_id_query
    UNION ALL
    SELECT 'full'::text, query_sid FROM v_best_alignment_id_query_hit
),
chain_sid AS (
    SELECT ps.sid AS pdb_sid, ps.pid, ps.cid, ps.stid
    FROM pdb_sequence ps
),
mapped AS (
    SELECT tu.tier, cs.pid, cs.cid, cs.stid
    FROM tier_union tu
    JOIN chain_sid cs ON cs.pdb_sid = tu.query_sid
)
SELECT
    ac.stid AS pdb_stid,
    ac.pdb_sequence_type,
    tu.tier,
    COUNT(DISTINCT (ac.pid, ac.cid)) AS total_chains,
    COUNT(DISTINCT (m.pid, m.cid))   AS mapped_chains_any
FROM all_chains ac
CROSS JOIN (VALUES ('identity-only'), ('no-hit'), ('full')) AS tu(tier)
LEFT JOIN mapped m
  ON m.pid = ac.pid AND m.cid = ac.cid AND m.stid = ac.stid AND m.tier = tu.tier
GROUP BY ac.stid, ac.pdb_sequence_type, tu.tier;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_sequence_totals_chain CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_sequence_totals_chain",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_sequence_totals_chain",
        dependencies=[
            "v_best_alignment_id",
            "v_best_alignment_id_query",
            "v_best_alignment_id_query_hit",
        ],
        description="Chain-level totals & mapped chain counts per tier and sequence type (precomputed).",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # NEW: Unified alignment sequence metrics MV (interval + right-bound in one structure)
    # ------------------------------------------------------------------
    "mv_alignment_sequence_metrics": SQLArtifact(
        name="mv_alignment_sequence_metrics",
        sql="""CREATE MATERIALIZED VIEW mv_alignment_sequence_metrics AS
WITH tier_union AS (
    SELECT 'identity-only'::text AS tier, * FROM v_best_alignment_id
    UNION ALL
    SELECT 'no-hit', * FROM v_best_alignment_id_query
    UNION ALL
    SELECT 'full', * FROM v_best_alignment_id_query_hit
),
pair_base AS (
    SELECT
        tu.tier,
        ps.sid AS pdb_sid,
        st.stid AS pdb_stid,
        st.type AS pdb_sequence_type,
        tu.hit_sid,
        tu.percent_identity
    FROM tier_union tu
    JOIN pdb_sequence ps ON ps.sid = tu.query_sid
    JOIN sequence_type st ON st.stid = ps.stid
),
chembl_type_map AS (
    SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
    UNION ALL
    SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
),
joined AS (
    SELECT pb.*, ctm.chembl_type
    FROM pair_base pb
    JOIN chembl_type_map ctm ON ctm.hit_sid = pb.hit_sid
),
interval_class AS (
    SELECT *,
        CASE
          WHEN percent_identity = 100 THEN '[100]'
          WHEN percent_identity >= 95 THEN '[95,100)'
          WHEN percent_identity >= 90 THEN '[90,95)'
          WHEN percent_identity >= 85 THEN '[85,90)'
          ELSE '[80,85)'
        END AS identity_interval
    FROM joined
),
right_bound_flags AS (
    SELECT *,
        (percent_identity >= 100)::int AS ge_100,
        (percent_identity >= 95)::int  AS ge_95,
        (percent_identity >= 90)::int  AS ge_90,
        (percent_identity >= 85)::int  AS ge_85,
        (percent_identity >= 80)::int  AS ge_80
    FROM interval_class
),
totals AS (
    SELECT
        tier, pdb_stid, pdb_sequence_type,
        COUNT(DISTINCT pdb_sid) AS total_sequences
    FROM right_bound_flags
    GROUP BY tier, pdb_stid, pdb_sequence_type
),
mapped_any AS (
    SELECT
        tier, chembl_type, pdb_stid, pdb_sequence_type,
        COUNT(DISTINCT pdb_sid) AS mapped_sequences_any
    FROM right_bound_flags
    GROUP BY tier, chembl_type, pdb_stid, pdb_sequence_type
),
interval_agg AS (
    SELECT
        tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval,
        COUNT(DISTINCT (pdb_sid, hit_sid)) AS raw_pairs_interval,
        COUNT(DISTINCT pdb_sid) AS distinct_sequences_interval
    FROM right_bound_flags
    GROUP BY tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval
),
right_bound_agg AS (
    SELECT
        tier, chembl_type, pdb_stid, pdb_sequence_type,
        COUNT(DISTINCT CASE WHEN ge_100=1 THEN (pdb_sid, hit_sid) END) AS raw_ge_100,
        COUNT(DISTINCT CASE WHEN ge_95 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_95,
        COUNT(DISTINCT CASE WHEN ge_90 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_90,
        COUNT(DISTINCT CASE WHEN ge_85 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_85,
        COUNT(DISTINCT CASE WHEN ge_80 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_80,
        COUNT(DISTINCT CASE WHEN ge_100=1 THEN pdb_sid END) AS distinct_ge_100,
        COUNT(DISTINCT CASE WHEN ge_95 =1 THEN pdb_sid END) AS distinct_ge_95,
        COUNT(DISTINCT CASE WHEN ge_90 =1 THEN pdb_sid END) AS distinct_ge_90,
        COUNT(DISTINCT CASE WHEN ge_85 =1 THEN pdb_sid END) AS distinct_ge_85,
        COUNT(DISTINCT CASE WHEN ge_80 =1 THEN pdb_sid END) AS distinct_ge_80
    FROM right_bound_flags
    GROUP BY tier, chembl_type, pdb_stid, pdb_sequence_type
)
SELECT
  ia.tier,
  ia.chembl_type,
  ia.pdb_stid,
  ia.pdb_sequence_type,
  ia.identity_interval,
  ia.raw_pairs_interval,
  ia.distinct_sequences_interval,
  t.total_sequences,
  ma.mapped_sequences_any,
  (ia.raw_pairs_interval::float / NULLIF(t.total_sequences,0)) AS interval_density,
  (ia.distinct_sequences_interval::float / NULLIF(t.total_sequences,0)) AS interval_occupancy,
  (ia.distinct_sequences_interval::float / NULLIF(ma.mapped_sequences_any,0)) AS interval_cond_occupancy,
  rba.raw_ge_100, rba.raw_ge_95, rba.raw_ge_90, rba.raw_ge_85, rba.raw_ge_80,
  rba.distinct_ge_100, rba.distinct_ge_95, rba.distinct_ge_90, rba.distinct_ge_85, rba.distinct_ge_80
FROM interval_agg ia
JOIN totals t ON t.tier = ia.tier AND t.pdb_stid = ia.pdb_stid AND t.pdb_sequence_type = ia.pdb_sequence_type
JOIN mapped_any ma ON ma.tier = ia.tier AND ma.pdb_stid = ia.pdb_stid AND ma.pdb_sequence_type = ia.pdb_sequence_type AND ma.chembl_type = ia.chembl_type
JOIN right_bound_agg rba ON rba.tier = ia.tier AND rba.chembl_type = ia.chembl_type
                        AND rba.pdb_stid = ia.pdb_stid AND rba.pdb_sequence_type = ia.pdb_sequence_type;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_alignment_sequence_metrics CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_alignment_sequence_metrics",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_alignment_sequence_metrics",
        dependencies=[
            "v_best_alignment_id",
            "v_best_alignment_id_query",
            "v_best_alignment_id_query_hit",
        ],
        description="Unified interval + right-bound sequence pair metrics per tier, chembl_type, and PDB sequence type.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # ------------------------------------------------------------------
    # NEW: Chain-level alignment sequence metrics MV (intersection of chains having all three PDB sequence types)
    # Mirrors mv_alignment_sequence_metrics but restricts to (pid,cid) that possess ATOMSEQ, SEQRES, SEQADV.
    # Produces identical column schema so reporting code can switch sources transparently.
    # ------------------------------------------------------------------
    "mv_alignment_sequence_metrics_chain": SQLArtifact(
        name="mv_alignment_sequence_metrics_chain",
        sql="""CREATE MATERIALIZED VIEW mv_alignment_sequence_metrics_chain AS
WITH intersection_chains AS (
    SELECT ps.pid, ps.cid
    FROM pdb_sequence ps
    WHERE ps.stid IN (1,4,7)  -- ATOMSEQ, SEQRES, SEQADV ids
    GROUP BY ps.pid, ps.cid
    HAVING COUNT(DISTINCT ps.stid) = 3
),
tier_union AS (
    SELECT 'identity-only'::text AS tier, * FROM v_best_alignment_id
    UNION ALL
    SELECT 'no-hit', * FROM v_best_alignment_id_query
    UNION ALL
    SELECT 'full', * FROM v_best_alignment_id_query_hit
),
pair_base AS (
    SELECT
        tu.tier,
        ps.sid AS pdb_sid,
        st.stid AS pdb_stid,
        st.type AS pdb_sequence_type,
        tu.hit_sid,
        tu.percent_identity
    FROM tier_union tu
    JOIN pdb_sequence ps ON ps.sid = tu.query_sid
    JOIN sequence_type st ON st.stid = ps.stid
    JOIN intersection_chains ic ON ic.pid = ps.pid AND ic.cid = ps.cid
),
chembl_type_map AS (
    SELECT ccs.sid AS hit_sid, 'COMPONENT'::text AS chembl_type FROM chembl_component_sequence ccs
    UNION ALL
    SELECT cvs.sid AS hit_sid, 'VARIANT'::text  AS chembl_type FROM chembl_variant_sequence cvs
),
joined AS (
    SELECT pb.*, ctm.chembl_type
    FROM pair_base pb
    JOIN chembl_type_map ctm ON ctm.hit_sid = pb.hit_sid
),
interval_class AS (
    SELECT *,
        CASE
          WHEN percent_identity = 100 THEN '[100]'
          WHEN percent_identity >= 95 THEN '[95,100)'
          WHEN percent_identity >= 90 THEN '[90,95)'
          WHEN percent_identity >= 85 THEN '[85,90)'
          ELSE '[80,85)'
        END AS identity_interval
    FROM joined
),
right_bound_flags AS (
    SELECT *,
        (percent_identity >= 100)::int AS ge_100,
        (percent_identity >= 95)::int  AS ge_95,
        (percent_identity >= 90)::int  AS ge_90,
        (percent_identity >= 85)::int  AS ge_85,
        (percent_identity >= 80)::int  AS ge_80
    FROM interval_class
),
totals AS (
    SELECT
        tier, pdb_stid, pdb_sequence_type,
        COUNT(DISTINCT pdb_sid) AS total_sequences
    FROM right_bound_flags
    GROUP BY tier, pdb_stid, pdb_sequence_type
),
mapped_any AS (
    SELECT
        tier, chembl_type, pdb_stid, pdb_sequence_type,
        COUNT(DISTINCT pdb_sid) AS mapped_sequences_any
    FROM right_bound_flags
    GROUP BY tier, chembl_type, pdb_stid, pdb_sequence_type
),
interval_agg AS (
    SELECT
        tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval,
        COUNT(DISTINCT (pdb_sid, hit_sid)) AS raw_pairs_interval,
        COUNT(DISTINCT pdb_sid) AS distinct_sequences_interval
    FROM right_bound_flags
    GROUP BY tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval
),
right_bound_agg AS (
    SELECT
        tier, chembl_type, pdb_stid, pdb_sequence_type,
        COUNT(DISTINCT CASE WHEN ge_100=1 THEN (pdb_sid, hit_sid) END) AS raw_ge_100,
        COUNT(DISTINCT CASE WHEN ge_95 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_95,
        COUNT(DISTINCT CASE WHEN ge_90 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_90,
        COUNT(DISTINCT CASE WHEN ge_85 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_85,
        COUNT(DISTINCT CASE WHEN ge_80 =1 THEN (pdb_sid, hit_sid) END) AS raw_ge_80,
        COUNT(DISTINCT CASE WHEN ge_100=1 THEN pdb_sid END) AS distinct_ge_100,
        COUNT(DISTINCT CASE WHEN ge_95 =1 THEN pdb_sid END) AS distinct_ge_95,
        COUNT(DISTINCT CASE WHEN ge_90 =1 THEN pdb_sid END) AS distinct_ge_90,
        COUNT(DISTINCT CASE WHEN ge_85 =1 THEN pdb_sid END) AS distinct_ge_85,
        COUNT(DISTINCT CASE WHEN ge_80 =1 THEN pdb_sid END) AS distinct_ge_80
    FROM right_bound_flags
    GROUP BY tier, chembl_type, pdb_stid, pdb_sequence_type
)
SELECT
  ia.tier,
  ia.chembl_type,
  ia.pdb_stid,
  ia.pdb_sequence_type,
  ia.identity_interval,
  ia.raw_pairs_interval,
  ia.distinct_sequences_interval,
  t.total_sequences,
  ma.mapped_sequences_any,
  (ia.raw_pairs_interval::float / NULLIF(t.total_sequences,0)) AS interval_density,
  (ia.distinct_sequences_interval::float / NULLIF(t.total_sequences,0)) AS interval_occupancy,
  (ia.distinct_sequences_interval::float / NULLIF(ma.mapped_sequences_any,0)) AS interval_cond_occupancy,
  rba.raw_ge_100, rba.raw_ge_95, rba.raw_ge_90, rba.raw_ge_85, rba.raw_ge_80,
  rba.distinct_ge_100, rba.distinct_ge_95, rba.distinct_ge_90, rba.distinct_ge_85, rba.distinct_ge_80
FROM interval_agg ia
JOIN totals t ON t.tier = ia.tier AND t.pdb_stid = ia.pdb_stid AND t.pdb_sequence_type = ia.pdb_sequence_type
JOIN mapped_any ma ON ma.tier = ia.tier AND ma.pdb_stid = ia.pdb_stid AND ma.pdb_sequence_type = ia.pdb_sequence_type AND ma.chembl_type = ia.chembl_type
JOIN right_bound_agg rba ON rba.tier = ia.tier AND rba.chembl_type = ia.chembl_type
                        AND rba.pdb_stid = ia.pdb_stid AND rba.pdb_sequence_type = ia.pdb_sequence_type;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_alignment_sequence_metrics_chain CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_alignment_sequence_metrics_chain",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_alignment_sequence_metrics_chain",
        dependencies=[
            "v_best_alignment_id",
            "v_best_alignment_id_query",
            "v_best_alignment_id_query_hit",
        ],
        description="Chain-level (intersection) interval + right-bound metrics per tier, chembl_type, and PDB sequence type.",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # Added indexes for chain-level alignment metrics MV
    "ix_mvasmc_tier_chembl_type": SQLArtifact(
        name="ix_mvasmc_tier_chembl_type",
        sql="CREATE INDEX IF NOT EXISTS ix_mvasmc_tier_chembl_type ON mv_alignment_sequence_metrics_chain (tier, chembl_type)",
        drop_sql="DROP INDEX IF EXISTS ix_mvasmc_tier_chembl_type",
        dependencies=["mv_alignment_sequence_metrics_chain"],
        description="Index for chain-level MV queries filtering by tier and chembl_type",
        artifact_type=ArtifactType.INDEX,
    ),
    "ix_mvasmc_dim_interval": SQLArtifact(
        name="ix_mvasmc_dim_interval",
        sql="CREATE INDEX IF NOT EXISTS ix_mvasmc_dim_interval ON mv_alignment_sequence_metrics_chain (tier, chembl_type, pdb_stid, pdb_sequence_type, identity_interval)",
        drop_sql="DROP INDEX IF EXISTS ix_mvasmc_dim_interval",
        dependencies=["mv_alignment_sequence_metrics_chain"],
        description="Composite dimensional + interval lookup for chain-level MV",
        artifact_type=ArtifactType.INDEX,
    ),
    "ix_mvasmc_dim": SQLArtifact(
        name="ix_mvasmc_dim",
        sql="CREATE INDEX IF NOT EXISTS ix_mvasmc_dim ON mv_alignment_sequence_metrics_chain (tier, pdb_stid, pdb_sequence_type)",
        drop_sql="DROP INDEX IF EXISTS ix_mvasmc_dim",
        dependencies=["mv_alignment_sequence_metrics_chain"],
        description="Support scans/joins on tier + pdb sequence type (chain-level MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_target_sequences": SQLArtifact(
        name="mv_ligand_target_sequences",
        sql="""CREATE MATERIALIZED VIEW mv_ligand_target_sequences AS
SELECT DISTINCT
    lat.clid,
    lat.chembl_tid,
    st.sid AS hit_sid
FROM mv_ligand_activity_target lat
JOIN sequence_target st ON st.chembl_tid = lat.chembl_tid;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_target_sequences CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_target_sequences",
        dependencies=["mv_ligand_activity_target"],
        description="Target sequences for ligand activities",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_best_alignment_identity80": SQLArtifact(
        name="mv_best_alignment_identity80",
        sql="""CREATE MATERIALIZED VIEW mv_best_alignment_identity80 AS
WITH ranked AS (
    SELECT
        bm.query_sid,
        bm.hit_sid,
        bm.percent_identity,
        bm.query_coverage,
        bm.hit_coverage,
        ROW_NUMBER() OVER (
            PARTITION BY bm.query_sid, bm.hit_sid
            ORDER BY bm.percent_identity DESC,
                     bm.query_coverage DESC NULLS LAST,
                     bm.hit_coverage DESC NULLS LAST
        ) AS rn
    FROM blast_match bm
    WHERE bm.percent_identity >= 80
)
SELECT
    query_sid,
    hit_sid,
    percent_identity,
    query_coverage,
    hit_coverage
FROM ranked
WHERE rn = 1;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_best_alignment_identity80 CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_best_alignment_identity80",
        dependencies=[],
        description="Best alignments with identity >= 80",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_ligand_pair_alignment_raw": SQLArtifact(
        name="mv_ligand_pair_alignment_raw",
        sql="""CREATE MATERIALIZED VIEW mv_ligand_pair_alignment_raw AS
WITH chain_map AS (
    SELECT ps.sid AS query_sid, ps.pid, ps.cid
    FROM pdb_sequence ps
),
ligand_bridge AS (
    SELECT plpi.pid, plpi.cid, plpi.pdblid
    FROM pdb_ligand_pdb_id plpi
),
structural_pairs AS (
    SELECT DISTINCT pdblid, clid, smmclid
    FROM pdb_ligand_chembl_ligand
)
SELECT
    sp.pdblid,
    sp.clid,
    sp.smmclid,
    ba.percent_identity
FROM mv_ligand_target_sequences lts
JOIN mv_best_alignment_identity80 ba ON ba.hit_sid = lts.hit_sid
JOIN chain_map cm ON cm.query_sid = ba.query_sid
JOIN ligand_bridge lb ON lb.pid = cm.pid AND lb.cid = cm.cid
JOIN structural_pairs sp ON sp.pdblid = lb.pdblid AND sp.clid = lts.clid;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_pair_alignment_raw CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_pair_alignment_raw",
        dependencies=["mv_ligand_target_sequences", "mv_best_alignment_identity80"],
        description="Raw pair alignment rows for ligand activity qualification",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_ligand_supported_pair_best_identity": SQLArtifact(
        name="mv_ligand_supported_pair_best_identity",
        sql="""CREATE MATERIALIZED VIEW mv_ligand_supported_pair_best_identity AS
WITH agg AS (
    SELECT
        pdblid,
        clid,
        MAX(percent_identity) AS max_percent_identity,
        MAX(smmclid)          AS best_smmclid,
        COUNT(*)              AS contributing_alignment_rows
    FROM mv_ligand_pair_alignment_raw
    GROUP BY pdblid, clid
)
SELECT
    pdblid,
    clid,
    best_smmclid,
    max_percent_identity,
    contributing_alignment_rows,
    CASE
        WHEN max_percent_identity = 100 THEN '[100]'
        WHEN max_percent_identity >= 95 THEN '[95,100)'
        WHEN max_percent_identity >= 90 THEN '[90,95)'
        WHEN max_percent_identity >= 85 THEN '[85,90)'
        ELSE '[80,85)'
    END AS identity_interval
FROM agg;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_supported_pair_best_identity CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_supported_pair_best_identity",
        dependencies=["mv_ligand_pair_alignment_raw"],
        description="Aggregated best identity per pair",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    "mv_ligand_supported_pair_identity_multibucket": SQLArtifact(
        name="mv_ligand_supported_pair_identity_multibucket",
        sql="""CREATE MATERIALIZED VIEW mv_ligand_supported_pair_identity_multibucket AS
WITH pair_intervals AS (
    SELECT DISTINCT
        r.pdblid,
        r.clid,
        r.smmclid,
        CASE
            WHEN r.percent_identity = 100 THEN '[100]'
            WHEN r.percent_identity >= 95 THEN '[95,100)'
            WHEN r.percent_identity >= 90 THEN '[90,95)'
            WHEN r.percent_identity >= 85 THEN '[85,90)'
            ELSE '[80,85)'
        END AS identity_interval
    FROM mv_ligand_pair_alignment_raw r
)
SELECT
    pdblid,
    clid,
    smmclid,
    identity_interval
FROM pair_intervals;""",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_supported_pair_identity_multibucket CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_supported_pair_identity_multibucket",
        dependencies=["mv_ligand_pair_alignment_raw"],
        description="Multi-bucket pair identity presence expansion",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # Views
    "v_ligand_sequence_target_fact_buckets": SQLArtifact(
        name="v_ligand_sequence_target_fact_buckets",
        sql="""CREATE OR REPLACE VIEW v_ligand_sequence_target_fact_buckets AS
SELECT
    mv.*,
    -- Disjoint identity buckets (non-overlapping intervals)
    CASE
        WHEN mv.percent_identity = 100::numeric THEN '100'
        WHEN mv.percent_identity >= 95::numeric THEN '95-<100'
        WHEN mv.percent_identity >= 90::numeric THEN '90-<95'
        WHEN mv.percent_identity >= 85::numeric THEN '85-<90'
        ELSE '80-<85'
    END AS identity_bucket_disjoint,
    -- Cumulative identity buckets (threshold-based)
    CASE
        WHEN mv.percent_identity = 100::numeric THEN '>=100'
        WHEN mv.percent_identity >= 95::numeric THEN '>=95'
        WHEN mv.percent_identity >= 90::numeric THEN '>=90'
        WHEN mv.percent_identity >= 85::numeric THEN '>=85'
        ELSE '>=80'
    END AS identity_bucket_cumulative,
    -- Additional useful calculations
    CASE
        WHEN mv.percent_identity = 100::numeric THEN 1
        WHEN mv.percent_identity >= 95::numeric THEN 2
        WHEN mv.percent_identity >= 90::numeric THEN 3
        WHEN mv.percent_identity >= 85::numeric THEN 4
        ELSE 5
    END AS identity_bucket_order
FROM mv_ligand_sequence_target_fact mv;""",
        drop_sql="DROP VIEW IF EXISTS v_ligand_sequence_target_fact_buckets CASCADE",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Identity bucket classification view",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_ligand_sequence_target_fact_activity": SQLArtifact(
        name="v_ligand_sequence_target_fact_activity",
        sql="""CREATE OR REPLACE VIEW v_ligand_sequence_target_fact_activity AS
SELECT
    f.query_sid,
    f.pdb_stid,
    f.hit_sid,
    f.chembl_type,
    f.chembl_tid,
    f.pdblid,
    f.clid,
    f.smmclid,
    f.percent_identity,
    f.query_coverage,
    f.hit_coverage,
    a.pchembl_value,
    a.aid,
    ay.asid
FROM mv_ligand_sequence_target_fact f
INNER JOIN chembl_ligand cl ON cl.clid = f.clid
INNER JOIN activity a ON a.molregno = cl.molregno
INNER JOIN assay ay ON ay.asid = a.asid
INNER JOIN target t ON t.tid = ay.tid AND t.chembl_tid = f.chembl_tid
WHERE a.pchembl_value IS NOT NULL;""",
        drop_sql="DROP VIEW IF EXISTS v_ligand_sequence_target_fact_activity CASCADE",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Activity-supported view (strict)",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_ligand_sequence_target_fact_activity_nohitcov": SQLArtifact(
        name="v_ligand_sequence_target_fact_activity_nohitcov",
        sql="""CREATE OR REPLACE VIEW v_ligand_sequence_target_fact_activity_nohitcov AS
SELECT
    f.query_sid,
    f.pdb_stid,
    f.hit_sid,
    f.chembl_type,
    f.chembl_tid,
    f.pdblid,
    f.clid,
    f.smmclid,
    f.percent_identity,
    f.query_coverage,
    f.hit_coverage,
    a.pchembl_value,
    a.aid,
    ay.asid
FROM mv_ligand_sequence_target_fact_nohitcov f
INNER JOIN chembl_ligand cl ON cl.clid = f.clid
INNER JOIN activity a ON a.molregno = cl.molregno
INNER JOIN assay ay ON ay.asid = a.asid
INNER JOIN target t ON t.tid = ay.tid AND t.chembl_tid = f.chembl_tid
WHERE a.pchembl_value IS NOT NULL;""",
        drop_sql="DROP VIEW IF EXISTS v_ligand_sequence_target_fact_activity_nohitcov CASCADE",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Activity-supported view (relaxed)",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_ligand_mapping_retention_summary": SQLArtifact(
        name="v_ligand_mapping_retention_summary",
        sql="""CREATE OR REPLACE VIEW v_ligand_mapping_retention_summary AS
WITH baseline AS (
    -- Level 1 baseline counts from structure-only mapping
    SELECT
        smmclid,
        COUNT(DISTINCT pdblid) AS baseline_pdblid_count,
        COUNT(DISTINCT (pdblid, clid)) AS baseline_pair_count
    FROM pdb_ligand_chembl_ligand
    GROUP BY smmclid
),
combined_data AS (
    -- Combine strict and relaxed data in one pass
    SELECT
        smmclid,
        CASE
            WHEN percent_identity = 100 THEN '>=100'
            WHEN percent_identity >= 95 THEN '>=95'
            WHEN percent_identity >= 90 THEN '>=90'
            WHEN percent_identity >= 85 THEN '>=85'
            ELSE '>=80'
        END AS identity_bucket_cumulative,
        pdblid,
        clid,
        'STRICT' AS variant
    FROM mv_ligand_sequence_target_fact
    UNION ALL
    SELECT
        smmclid,
        CASE
            WHEN percent_identity = 100 THEN '>=100'
            WHEN percent_identity >= 95 THEN '>=95'
            WHEN percent_identity >= 90 THEN '>=90'
            WHEN percent_identity >= 85 THEN '>=85'
            ELSE '>=80'
        END AS identity_bucket_cumulative,
        pdblid,
        clid,
        'RELAXED' AS variant
    FROM mv_ligand_sequence_target_fact_nohitcov
),
aggregated AS (
    SELECT
        smmclid,
        identity_bucket_cumulative,
        COUNT(DISTINCT pdblid) FILTER (WHERE variant = 'STRICT') AS surviving_pdblid_count_strict,
        COUNT(DISTINCT pdblid) FILTER (WHERE variant = 'RELAXED') AS surviving_pdblid_count_relaxed,
        COUNT(DISTINCT (pdblid, clid)) FILTER (WHERE variant = 'STRICT') AS surviving_pair_count_strict,
        COUNT(DISTINCT (pdblid, clid)) FILTER (WHERE variant = 'RELAXED') AS surviving_pair_count_relaxed
    FROM combined_data
    GROUP BY smmclid, identity_bucket_cumulative
)
SELECT
    a.smmclid,
    a.identity_bucket_cumulative,
    b.baseline_pdblid_count,
    b.baseline_pair_count,
    a.surviving_pdblid_count_strict,
    a.surviving_pdblid_count_relaxed,
    a.surviving_pair_count_strict,
    a.surviving_pair_count_relaxed,
    ROUND(a.surviving_pdblid_count_strict * 100.0 / NULLIF(b.baseline_pdblid_count, 0), 2) AS strict_retention_pdblid_pct,
    ROUND(a.surviving_pdblid_count_relaxed * 100.0 / NULLIF(b.baseline_pdblid_count, 0), 2) AS relaxed_retention_pdblid_pct,
    ROUND(a.surviving_pair_count_strict * 100.0 / NULLIF(b.baseline_pair_count, 0), 2) AS strict_retention_pair_pct,
    ROUND(a.surviving_pair_count_relaxed * 100.0 / NULLIF(b.baseline_pair_count, 0), 2) AS relaxed_retention_pair_pct,
    ROUND((a.surviving_pdblid_count_relaxed - a.surviving_pdblid_count_strict) * 100.0 / NULLIF(b.baseline_pdblid_count, 0), 2) AS retention_gain_pdblid_pp,
    ROUND((a.surviving_pair_count_relaxed - a.surviving_pair_count_strict) * 100.0 / NULLIF(b.baseline_pair_count, 0), 2) AS retention_gain_pair_pp
FROM aggregated a
JOIN baseline b ON a.smmclid = b.smmclid
ORDER BY
    a.smmclid,
    CASE a.identity_bucket_cumulative
        WHEN '>=100' THEN 1
        WHEN '>=95' THEN 2
        WHEN '>=90' THEN 3
        WHEN '>=85' THEN 4
        WHEN '>=80' THEN 5
    END;""",
        drop_sql="DROP VIEW IF EXISTS v_ligand_mapping_retention_summary CASCADE",
        dependencies=[
            "mv_ligand_sequence_target_fact",
            "mv_ligand_sequence_target_fact_nohitcov",
        ],
        description="Complete retention summary view",
        artifact_type=ArtifactType.VIEW,
    ),
    "v_ligand_mapping_partition_summary": SQLArtifact(
        name="v_ligand_mapping_partition_summary",
        sql="""CREATE OR REPLACE VIEW v_ligand_mapping_partition_summary AS
WITH combined_partitions AS (
    SELECT
        'STRICT' AS variant,
        smmclid,
        pdb_stid,
        chembl_type,
        COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid)) AS quadruplets,
        COUNT(DISTINCT (pdblid, clid, chembl_tid)) AS triplets,
        COUNT(DISTINCT (pdblid, clid)) AS pairs,
        COUNT(DISTINCT pdblid) AS pdblid_count
    FROM mv_ligand_sequence_target_fact
    GROUP BY smmclid, pdb_stid, chembl_type
    UNION ALL
    SELECT
        'RELAXED' AS variant,
        smmclid,
        pdb_stid,
        chembl_type,
        COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid)) AS quadruplets,
        COUNT(DISTINCT (pdblid, clid, chembl_tid)) AS triplets,
        COUNT(DISTINCT (pdblid, clid)) AS pairs,
        COUNT(DISTINCT pdblid) AS pdblid_count
    FROM mv_ligand_sequence_target_fact_nohitcov
    GROUP BY smmclid, pdb_stid, chembl_type
)
SELECT
    smmclid,
    pdb_stid,
    chembl_type,
    MAX(quadruplets) FILTER (WHERE variant = 'STRICT') AS strict_quadruplets,
    MAX(triplets) FILTER (WHERE variant = 'STRICT') AS strict_triplets,
    MAX(pairs) FILTER (WHERE variant = 'STRICT') AS strict_pairs,
    MAX(pdblid_count) FILTER (WHERE variant = 'STRICT') AS strict_pdblid,
    MAX(quadruplets) FILTER (WHERE variant = 'RELAXED') AS relaxed_quadruplets,
    MAX(triplets) FILTER (WHERE variant = 'RELAXED') AS relaxed_triplets,
    MAX(pairs) FILTER (WHERE variant = 'RELAXED') AS relaxed_pairs,
    MAX(pdblid_count) FILTER (WHERE variant = 'RELAXED') AS relaxed_pdblid,
    ROUND(MAX(quadruplets) FILTER (WHERE variant = 'RELAXED')::numeric /
          NULLIF(MAX(quadruplets) FILTER (WHERE variant = 'STRICT'), 0), 2) AS quadruplet_expansion_factor,
    ROUND(MAX(pdblid_count) FILTER (WHERE variant = 'RELAXED')::numeric /
          NULLIF(MAX(pdblid_count) FILTER (WHERE variant = 'STRICT'), 0), 2) AS pdblid_expansion_factor
FROM combined_partitions
GROUP BY smmclid, pdb_stid, chembl_type
ORDER BY smmclid, pdb_stid, chembl_type;""",
        drop_sql="DROP VIEW IF EXISTS v_ligand_mapping_partition_summary CASCADE",
        dependencies=[
            "mv_ligand_sequence_target_fact",
            "mv_ligand_sequence_target_fact_nohitcov",
        ],
        description="Partitioned quadruplet summary view",
        artifact_type=ArtifactType.VIEW,
    ),
    # Selected Indexes (most critical ones)
    "mv_ligand_sequence_target_fact_query_sid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_query_sid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_query_sid_idx ON mv_ligand_sequence_target_fact (query_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_query_sid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Primary index on query_sid",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_chembl_tid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_chembl_tid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_chembl_tid_idx ON mv_ligand_sequence_target_fact (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_chembl_tid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Primary index on chembl_tid",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_identity_only": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only",
        sql="CREATE MATERIALIZED VIEW mv_ligand_sequence_target_fact_identity_only AS SELECT * FROM f_ligand_sequence_target_base(80, FALSE, NULL, FALSE, NULL);",
        drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_ligand_sequence_target_fact_identity_only CASCADE",
        refresh_sql="REFRESH MATERIALIZED VIEW mv_ligand_sequence_target_fact_identity_only",
        refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_ligand_sequence_target_fact_identity_only",
        dependencies=["f_ligand_sequence_target_base"],
        description="Identity-only variant (no query or hit coverage filters) requiring non-null metrics",
        artifact_type=ArtifactType.MATERIALIZED_VIEW,
    ),
    # Indexes for identity-only MV
    "mv_ligand_sequence_target_fact_identity_only_query_sid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only_query_sid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_identity_only_query_sid_idx ON mv_ligand_sequence_target_fact_identity_only (query_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_identity_only_query_sid_idx",
        dependencies=["mv_ligand_sequence_target_fact_identity_only"],
        description="Primary index on query_sid (identity-only MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_identity_only_chembl_tid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only_chembl_tid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_identity_only_chembl_tid_idx ON mv_ligand_sequence_target_fact_identity_only (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_identity_only_chembl_tid_idx",
        dependencies=["mv_ligand_sequence_target_fact_identity_only"],
        description="Primary index on chembl_tid (identity-only MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_identity_only_quadruplet_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only_quadruplet_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_identity_only_quadruplet_idx ON mv_ligand_sequence_target_fact_identity_only (pdblid, query_sid, clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_identity_only_quadruplet_idx",
        dependencies=["mv_ligand_sequence_target_fact_identity_only"],
        description="Composite quadruplet index (identity-only MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_quadruplet_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_quadruplet_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_quadruplet_idx ON mv_ligand_sequence_target_fact (pdblid, query_sid, clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_quadruplet_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Composite index for quadruplet queries",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_query_sid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_query_sid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_query_sid_idx ON mv_ligand_sequence_target_fact_nohitcov (query_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_query_sid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Primary index on query_sid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_chembl_tid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_chembl_tid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_chembl_tid_idx ON mv_ligand_sequence_target_fact_nohitcov (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_chembl_tid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Primary index on chembl_tid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_active_pair_ligands_clid_uq": SQLArtifact(
        name="mv_active_pair_ligands_clid_uq",
        sql="CREATE UNIQUE INDEX mv_active_pair_ligands_clid_uq ON mv_active_pair_ligands (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_active_pair_ligands_clid_uq",
        dependencies=["mv_active_pair_ligands"],
        description="Unique index on clid",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_supported_pair_best_identity_pair_idx": SQLArtifact(
        name="mv_ligand_supported_pair_best_identity_pair_idx",
        sql="CREATE INDEX mv_ligand_supported_pair_best_identity_pair_idx ON mv_ligand_supported_pair_best_identity (pdblid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_supported_pair_best_identity_pair_idx",
        dependencies=["mv_ligand_supported_pair_best_identity"],
        description="Composite index on ligand pair",
        artifact_type=ArtifactType.INDEX,
    ),
    # --------------------------------------------------------------------------
    # Added: Full strict MV index coverage (ported from earlier DDL)
    # --------------------------------------------------------------------------
    "mv_ligand_sequence_target_fact_pdb_stid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_pdb_stid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_pdb_stid_idx ON mv_ligand_sequence_target_fact (pdb_stid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_pdb_stid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on pdb_stid (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_hit_sid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_hit_sid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_hit_sid_idx ON mv_ligand_sequence_target_fact (hit_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_hit_sid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on hit_sid (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_chembl_type_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_chembl_type_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_chembl_type_idx ON mv_ligand_sequence_target_fact (chembl_type);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_chembl_type_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on chembl_type (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_pdblid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_pdblid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_pdblid_idx ON mv_ligand_sequence_target_fact (pdblid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_pdblid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on pdblid (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_clid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_clid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_clid_idx ON mv_ligand_sequence_target_fact (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_clid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on clid (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_smmclid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_smmclid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_smmclid_idx ON mv_ligand_sequence_target_fact (smmclid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_smmclid_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on smmclid (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_percent_identity_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_percent_identity_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_percent_identity_idx ON mv_ligand_sequence_target_fact (percent_identity);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_percent_identity_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Index on percent_identity (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_triplet_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_triplet_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_triplet_idx ON mv_ligand_sequence_target_fact (pdblid, clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_triplet_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Composite triplet index (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_pair_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_pair_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_pair_idx ON mv_ligand_sequence_target_fact (pdblid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_pair_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Composite pair index (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_partition_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_partition_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_partition_idx ON mv_ligand_sequence_target_fact (smmclid, pdb_stid, chembl_type);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_partition_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Partition-style composite index (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_retention_covering_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_retention_covering_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_retention_covering_idx ON mv_ligand_sequence_target_fact (smmclid, percent_identity, pdblid, clid) INCLUDE (pdb_stid, chembl_type, query_sid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_retention_covering_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Covering index for retention / distribution queries (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_high_identity_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_high_identity_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_high_identity_idx ON mv_ligand_sequence_target_fact (smmclid, pdblid, clid) WHERE percent_identity >= 95;",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_high_identity_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Partial index for high-identity (>=95) filtering (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_component_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_component_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_component_idx ON mv_ligand_sequence_target_fact (smmclid, pdblid, clid) WHERE chembl_type = 'COMPONENT';",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_component_idx",
        dependencies=["mv_ligand_sequence_target_fact"],
        description="Partial index for component-only filtering (strict MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    # --------------------------------------------------------------------------
    # Added: Relaxed MV full index coverage
    # --------------------------------------------------------------------------
    "mv_ligand_sequence_target_fact_nohitcov_pdb_stid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_pdb_stid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_pdb_stid_idx ON mv_ligand_sequence_target_fact_nohitcov (pdb_stid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_pdb_stid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on pdb_stid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_hit_sid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_hit_sid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_hit_sid_idx ON mv_ligand_sequence_target_fact_nohitcov (hit_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_hit_sid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on hit_sid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_chembl_type_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_chembl_type_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_chembl_type_idx ON mv_ligand_sequence_target_fact_nohitcov (chembl_type);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_chembl_type_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on chembl_type (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_pdblid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_pdblid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_pdblid_idx ON mv_ligand_sequence_target_fact_nohitcov (pdblid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_pdblid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on pdblid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_clid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_clid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_clid_idx ON mv_ligand_sequence_target_fact_nohitcov (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_clid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on clid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_smmclid_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_smmclid_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_smmclid_idx ON mv_ligand_sequence_target_fact_nohitcov (smmclid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_smmclid_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on smmclid (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_percent_identity_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_percent_identity_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_percent_identity_idx ON mv_ligand_sequence_target_fact_nohitcov (percent_identity);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_percent_identity_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Index on percent_identity (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_quadruplet_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_quadruplet_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_quadruplet_idx ON mv_ligand_sequence_target_fact_nohitcov (pdblid, query_sid, clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_quadruplet_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Quadruplet composite index (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_triplet_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_triplet_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_triplet_idx ON mv_ligand_sequence_target_fact_nohitcov (pdblid, clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_triplet_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Triplet composite index (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_pair_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_pair_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_pair_idx ON mv_ligand_sequence_target_fact_nohitcov (pdblid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_pair_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Pair composite index (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_partition_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_partition_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_partition_idx ON mv_ligand_sequence_target_fact_nohitcov (smmclid, pdb_stid, chembl_type);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_partition_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Partition-style composite index (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_retention_covering_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_retention_covering_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_retention_covering_idx ON mv_ligand_sequence_target_fact_nohitcov (smmclid, percent_identity, pdblid, clid) INCLUDE (pdb_stid, chembl_type, query_sid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_retention_covering_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Covering index for retention queries (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_low_coverage_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_low_coverage_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_low_coverage_idx ON mv_ligand_sequence_target_fact_nohitcov (smmclid, pdblid, clid) WHERE hit_coverage < 80;",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_low_coverage_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Partial index for low hit_coverage (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_nohitcov_null_coverage_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_nohitcov_null_coverage_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_nohitcov_null_coverage_idx ON mv_ligand_sequence_target_fact_nohitcov (smmclid, pdblid, clid) WHERE hit_coverage IS NULL;",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_nohitcov_null_coverage_idx",
        dependencies=["mv_ligand_sequence_target_fact_nohitcov"],
        description="Partial index for NULL hit_coverage (relaxed MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    # --------------------------------------------------------------------------
    # Identity-only MV supplemental indexes
    # --------------------------------------------------------------------------
    "mv_ligand_sequence_target_fact_identity_only_percent_identity_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only_percent_identity_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_identity_only_percent_identity_idx ON mv_ligand_sequence_target_fact_identity_only (percent_identity);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_identity_only_percent_identity_idx",
        dependencies=["mv_ligand_sequence_target_fact_identity_only"],
        description="Index on percent_identity (identity-only MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_sequence_target_fact_identity_only_triplet_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only_triplet_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_identity_only_triplet_idx ON mv_ligand_sequence_target_fact_identity_only (pdblid, clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_identity_only_triplet_idx",
        dependencies=["mv_ligand_sequence_target_fact_identity_only"],
        description="Triplet composite index (identity-only MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    # ------------------------------------------------------------------
    # NEW: Indexes for ligand–assay convergence and presence MVs
    # ------------------------------------------------------------------
    "mv_pdb_ligand_assay_link_clid_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_clid_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_clid_idx ON mv_pdb_ligand_assay_link (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_clid_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Index on clid for ligand-assay link MV",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_pdb_ligand_assay_link_chembl_tid_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_chembl_tid_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_chembl_tid_idx ON mv_pdb_ligand_assay_link (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_chembl_tid_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Index on chembl_tid for ligand-assay link MV",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_pdb_ligand_assay_link_variant_id_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_variant_id_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_variant_id_idx ON mv_pdb_ligand_assay_link (variant_id);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_variant_id_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Index on variant_id for variant gating",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_pdb_ligand_assay_link_asid_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_asid_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_asid_idx ON mv_pdb_ligand_assay_link (asid);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_asid_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Index on assay id",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_pdb_ligand_assay_link_pid_cid_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_pid_cid_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_pid_cid_idx ON mv_pdb_ligand_assay_link (pid, cid);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_pid_cid_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Composite index for chain-specific ligand joins",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_pdb_ligand_assay_link_pdblid_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_pdblid_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_pdblid_idx ON mv_pdb_ligand_assay_link (pdblid);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_pdblid_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Index on structural ligand id",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_pdb_ligand_assay_link_route_composite_idx": SQLArtifact(
        name="mv_pdb_ligand_assay_link_route_composite_idx",
        sql="CREATE INDEX mv_pdb_ligand_assay_link_route_composite_idx ON mv_pdb_ligand_assay_link (chembl_tid, variant_id, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_pdb_ligand_assay_link_route_composite_idx",
        dependencies=["mv_pdb_ligand_assay_link"],
        description="Composite for (target, variant, ligand) lookups",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_component_presence_pdb_code_idx": SQLArtifact(
        name="mv_component_presence_pdb_code_idx",
        sql="CREATE INDEX mv_component_presence_pdb_code_idx ON mv_component_sequence_presence_identity_only (pdb_code);",
        drop_sql="DROP INDEX IF EXISTS mv_component_presence_pdb_code_idx",
        dependencies=["mv_component_sequence_presence_identity_only"],
        description="Index on pdb_code (component presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_component_presence_chembl_tid_idx": SQLArtifact(
        name="mv_component_presence_chembl_tid_idx",
        sql="CREATE INDEX mv_component_presence_chembl_tid_idx ON mv_component_sequence_presence_identity_only (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_component_presence_chembl_tid_idx",
        dependencies=["mv_component_sequence_presence_identity_only"],
        description="Index on chembl_tid (component presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_component_presence_clid_idx": SQLArtifact(
        name="mv_component_presence_clid_idx",
        sql="CREATE INDEX mv_component_presence_clid_idx ON mv_component_sequence_presence_identity_only (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_component_presence_clid_idx",
        dependencies=["mv_component_sequence_presence_identity_only"],
        description="Index on clid (component presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_component_presence_pdblid_idx": SQLArtifact(
        name="mv_component_presence_pdblid_idx",
        sql="CREATE INDEX mv_component_presence_pdblid_idx ON mv_component_sequence_presence_identity_only (pdblid);",
        drop_sql="DROP INDEX IF EXISTS mv_component_presence_pdblid_idx",
        dependencies=["mv_component_sequence_presence_identity_only"],
        description="Index on pdblid (component presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_component_presence_percent_identity_idx": SQLArtifact(
        name="mv_component_presence_percent_identity_idx",
        sql="CREATE INDEX mv_component_presence_percent_identity_idx ON mv_component_sequence_presence_identity_only (percent_identity DESC);",
        drop_sql="DROP INDEX IF EXISTS mv_component_presence_percent_identity_idx",
        dependencies=["mv_component_sequence_presence_identity_only"],
        description="Descending identity index (component presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_variant_presence_pdb_code_idx": SQLArtifact(
        name="mv_variant_presence_pdb_code_idx",
        sql="CREATE INDEX mv_variant_presence_pdb_code_idx ON mv_variant_sequence_presence_identity_only (pdb_code);",
        drop_sql="DROP INDEX IF EXISTS mv_variant_presence_pdb_code_idx",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Index on pdb_code (variant presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_variant_presence_chembl_tid_idx": SQLArtifact(
        name="mv_variant_presence_chembl_tid_idx",
        sql="CREATE INDEX mv_variant_presence_chembl_tid_idx ON mv_variant_sequence_presence_identity_only (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_variant_presence_chembl_tid_idx",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Index on chembl_tid (variant presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_variant_presence_clid_idx": SQLArtifact(
        name="mv_variant_presence_clid_idx",
        sql="CREATE INDEX mv_variant_presence_clid_idx ON mv_variant_sequence_presence_identity_only (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_variant_presence_clid_idx",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Index on clid (variant presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_variant_presence_pdblid_idx": SQLArtifact(
        name="mv_variant_presence_pdblid_idx",
        sql="CREATE INDEX mv_variant_presence_pdblid_idx ON mv_variant_sequence_presence_identity_only (pdblid);",
        drop_sql="DROP INDEX IF EXISTS mv_variant_presence_pdblid_idx",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Index on pdblid (variant presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_variant_presence_percent_identity_idx": SQLArtifact(
        name="mv_variant_presence_percent_identity_idx",
        sql="CREATE INDEX mv_variant_presence_percent_identity_idx ON mv_variant_sequence_presence_identity_only (percent_identity DESC);",
        drop_sql="DROP INDEX IF EXISTS mv_variant_presence_percent_identity_idx",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Descending identity index (variant presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_variant_presence_variant_id_idx": SQLArtifact(
        name="mv_variant_presence_variant_id_idx",
        sql="CREATE INDEX mv_variant_presence_variant_id_idx ON mv_variant_sequence_presence_identity_only (variant_id);",
        drop_sql="DROP INDEX IF EXISTS mv_variant_presence_variant_id_idx",
        dependencies=["mv_variant_sequence_presence_identity_only"],
        description="Index on variant_id (variant presence)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_component_variant_overlap_pdb_code_idx": SQLArtifact(
        name="mv_component_variant_overlap_pdb_code_idx",
        sql="CREATE INDEX mv_component_variant_overlap_pdb_code_idx ON mv_component_variant_presence_overlap (pdb_code);",
        drop_sql="DROP INDEX IF EXISTS mv_component_variant_overlap_pdb_code_idx",
        dependencies=["mv_component_variant_presence_overlap"],
        description="Index on pdb_code for overlap MV",
        artifact_type=ArtifactType.INDEX,
    ),
    # ------------------------------------------------------------------
    # Orchestration Function: Refresh new presence pipeline
    # ------------------------------------------------------------------
    "refresh_sequence_presence_mvs": SQLArtifact(
        name="refresh_sequence_presence_mvs",
        sql="""CREATE OR REPLACE FUNCTION refresh_sequence_presence_mvs(
    refresh_concurrently BOOLEAN DEFAULT FALSE,
    include_coverage BOOLEAN DEFAULT TRUE
) RETURNS TABLE(
    artifact TEXT,
    refreshed_at TIMESTAMP,
    rows BIGINT
) AS $$
DECLARE
    exec_stmt TEXT;
BEGIN
    -- 1. Foundational ligand/assay MVs
    IF refresh_concurrently THEN
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_pdb_ligand_chembl_activity_assay_base';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_pdb_ligand_assay_link';
    ELSE
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_pdb_ligand_chembl_activity_assay_base';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_pdb_ligand_assay_link';
    END IF;
    RETURN QUERY SELECT 'mv_pdb_ligand_chembl_activity_assay_base', clock_timestamp(), COUNT(*) FROM mv_pdb_ligand_chembl_activity_assay_base;
    RETURN QUERY SELECT 'mv_pdb_ligand_assay_link', clock_timestamp(), COUNT(*) FROM mv_pdb_ligand_assay_link;
    -- 2. Presence (identity-only)
    IF refresh_concurrently THEN
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_sequence_presence_identity_only';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_variant_sequence_presence_identity_only';
    ELSE
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_component_sequence_presence_identity_only';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_variant_sequence_presence_identity_only';
    END IF;
    RETURN QUERY SELECT 'mv_component_sequence_presence_identity_only', clock_timestamp(), COUNT(*) FROM mv_component_sequence_presence_identity_only;
    RETURN QUERY SELECT 'mv_variant_sequence_presence_identity_only', clock_timestamp(), COUNT(*) FROM mv_variant_sequence_presence_identity_only;
    -- 3. Optional coverage
    IF include_coverage THEN
        IF refresh_concurrently THEN
            EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_sequence_presence_cov';
            EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_variant_sequence_presence_cov';
        ELSE
            EXECUTE 'REFRESH MATERIALIZED VIEW mv_component_sequence_presence_cov';
            EXECUTE 'REFRESH MATERIALIZED VIEW mv_variant_sequence_presence_cov';
        END IF;
        RETURN QUERY SELECT 'mv_component_sequence_presence_cov', clock_timestamp(), COUNT(*) FROM mv_component_sequence_presence_cov;
        RETURN QUERY SELECT 'mv_variant_sequence_presence_cov', clock_timestamp(), COUNT(*) FROM mv_variant_sequence_presence_cov;
    END IF;
    -- 4. Overlap & counts
    IF refresh_concurrently THEN
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_variant_presence_overlap';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_component_presence_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_variant_presence_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_union_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_overlap_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_quadruplet_matrix';
        EXECUTE 'REFRESH MATERIALIZED VIEW CONCURRENTLY mv_presence_quadruplet_counts';
    ELSE
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_component_variant_presence_overlap';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_component_presence_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_variant_presence_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_presence_union_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_presence_overlap_counts';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_presence_quadruplet_matrix';
        EXECUTE 'REFRESH MATERIALIZED VIEW mv_presence_quadruplet_counts';
    END IF;
    RETURN QUERY SELECT 'mv_component_variant_presence_overlap', clock_timestamp(), COUNT(*) FROM mv_component_variant_presence_overlap;
    RETURN QUERY SELECT 'mv_component_presence_counts', clock_timestamp(), presence_rows FROM mv_component_presence_counts;
    RETURN QUERY SELECT 'mv_variant_presence_counts', clock_timestamp(), presence_rows FROM mv_variant_presence_counts;
    RETURN QUERY SELECT 'mv_presence_union_counts', clock_timestamp(), presence_rows FROM mv_presence_union_counts;
    RETURN QUERY SELECT 'mv_presence_overlap_counts', clock_timestamp(), overlap_rows FROM mv_presence_overlap_counts;
    RETURN QUERY SELECT 'mv_presence_quadruplet_matrix', clock_timestamp(), COUNT(*) FROM mv_presence_quadruplet_matrix;
    RETURN QUERY SELECT 'mv_presence_quadruplet_counts', clock_timestamp(), quadruplets FROM mv_presence_quadruplet_counts WHERE route = 'union';
END;
$$ LANGUAGE plpgsql;""",
        drop_sql="DROP FUNCTION IF EXISTS refresh_sequence_presence_mvs(BOOLEAN,BOOLEAN) CASCADE",
        dependencies=[
            "mv_pdb_ligand_chembl_activity_assay_base",
            "mv_pdb_ligand_assay_link",
            "mv_component_sequence_presence_identity_only",
            "mv_variant_sequence_presence_identity_only",
            "mv_component_sequence_presence_cov",
            "mv_variant_sequence_presence_cov",
            "mv_component_variant_presence_overlap",
            "mv_component_presence_counts",
            "mv_variant_presence_counts",
            "mv_presence_union_counts",
            "mv_presence_overlap_counts",
            "mv_presence_quadruplet_matrix",
            "mv_presence_quadruplet_counts",
        ],
        description="Refresh orchestrator for new ligand→assay→(component|variant) presence pipeline.",
        artifact_type=ArtifactType.FUNCTION,
    ),
    # ------------------------------------------------------------------
    # EXPLAIN TEMPLATE VIEW (reference queries for performance validation)
    # ------------------------------------------------------------------
    "v_sequence_presence_explain_templates": SQLArtifact(
        name="v_sequence_presence_explain_templates",
        sql="""CREATE OR REPLACE VIEW v_sequence_presence_explain_templates AS
SELECT * FROM (
    VALUES
      ('EXPLAIN (ANALYZE, BUFFERS) SELECT * FROM mv_component_sequence_presence_identity_only WHERE pdb_code = ''2W9S'';'),
      ('EXPLAIN (ANALYZE, BUFFERS) SELECT * FROM mv_variant_sequence_presence_identity_only WHERE chembl_tid = 12345 LIMIT 25;'),
      ('EXPLAIN (ANALYZE, BUFFERS) SELECT COUNT(*) FROM mv_component_variant_presence_overlap;'),
      ('EXPLAIN (ANALYZE, BUFFERS) SELECT * FROM mv_component_sequence_presence_cov WHERE percent_identity >= 90;'),
      ('EXPLAIN (ANALYZE, BUFFERS) SELECT pdb_code, COUNT(*) FROM mv_variant_sequence_presence_cov GROUP BY pdb_code ORDER BY COUNT(*) DESC LIMIT 10;')
) AS t(template_sql);""",
        drop_sql="DROP VIEW IF EXISTS v_sequence_presence_explain_templates CASCADE",
        dependencies=[
            "mv_component_sequence_presence_identity_only",
            "mv_variant_sequence_presence_identity_only",
            "mv_component_variant_presence_overlap",
            "mv_component_sequence_presence_cov",
            "mv_variant_sequence_presence_cov",
        ],
        description="Reference EXPLAIN templates for performance benchmarking of presence pipeline.",
        artifact_type=ArtifactType.VIEW,
    ),
    "mv_ligand_sequence_target_fact_identity_only_pair_idx": SQLArtifact(
        name="mv_ligand_sequence_target_fact_identity_only_pair_idx",
        sql="CREATE INDEX mv_ligand_sequence_target_fact_identity_only_pair_idx ON mv_ligand_sequence_target_fact_identity_only (pdblid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_sequence_target_fact_identity_only_pair_idx",
        dependencies=["mv_ligand_sequence_target_fact_identity_only"],
        description="Pair composite index (identity-only MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    # --------------------------------------------------------------------------
    # Pipeline MV index coverage (activity / alignment chain)
    # --------------------------------------------------------------------------
    "mv_active_pair_ligands_molregno_idx": SQLArtifact(
        name="mv_active_pair_ligands_molregno_idx",
        sql="CREATE INDEX mv_active_pair_ligands_molregno_idx ON mv_active_pair_ligands (molregno);",
        drop_sql="DROP INDEX IF EXISTS mv_active_pair_ligands_molregno_idx",
        dependencies=["mv_active_pair_ligands"],
        description="Index on molregno (active pair ligands)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_activity_target_clid_idx": SQLArtifact(
        name="mv_ligand_activity_target_clid_idx",
        sql="CREATE INDEX mv_ligand_activity_target_clid_idx ON mv_ligand_activity_target (clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_activity_target_clid_idx",
        dependencies=["mv_ligand_activity_target"],
        description="Index on clid (ligand activity target)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_activity_target_chembl_tid_idx": SQLArtifact(
        name="mv_ligand_activity_target_chembl_tid_idx",
        sql="CREATE INDEX mv_ligand_activity_target_chembl_tid_idx ON mv_ligand_activity_target (chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_activity_target_chembl_tid_idx",
        dependencies=["mv_ligand_activity_target"],
        description="Index on chembl_tid (ligand activity target)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_activity_target_aid_idx": SQLArtifact(
        name="mv_ligand_activity_target_aid_idx",
        sql="CREATE INDEX mv_ligand_activity_target_aid_idx ON mv_ligand_activity_target (aid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_activity_target_aid_idx",
        dependencies=["mv_ligand_activity_target"],
        description="Index on aid (ligand activity target)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_target_sequences_hit_sid_idx": SQLArtifact(
        name="mv_ligand_target_sequences_hit_sid_idx",
        sql="CREATE INDEX mv_ligand_target_sequences_hit_sid_idx ON mv_ligand_target_sequences (hit_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_target_sequences_hit_sid_idx",
        dependencies=["mv_ligand_target_sequences"],
        description="Index on hit_sid (target sequences)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_target_sequences_clid_tid_idx": SQLArtifact(
        name="mv_ligand_target_sequences_clid_tid_idx",
        sql="CREATE INDEX mv_ligand_target_sequences_clid_tid_idx ON mv_ligand_target_sequences (clid, chembl_tid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_target_sequences_clid_tid_idx",
        dependencies=["mv_ligand_target_sequences"],
        description="Composite (clid, chembl_tid) (target sequences)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_best_alignment_identity80_hit_sid_idx": SQLArtifact(
        name="mv_best_alignment_identity80_hit_sid_idx",
        sql="CREATE INDEX mv_best_alignment_identity80_hit_sid_idx ON mv_best_alignment_identity80 (hit_sid);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity80_hit_sid_idx",
        dependencies=["mv_best_alignment_identity80"],
        description="Index on hit_sid (best alignments)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_best_alignment_identity80_percent_idx": SQLArtifact(
        name="mv_best_alignment_identity80_percent_idx",
        sql="CREATE INDEX mv_best_alignment_identity80_percent_idx ON mv_best_alignment_identity80 (percent_identity DESC);",
        drop_sql="DROP INDEX IF EXISTS mv_best_alignment_identity80_percent_idx",
        dependencies=["mv_best_alignment_identity80"],
        description="Index on percent_identity DESC (best alignments)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_pair_alignment_raw_pair_idx": SQLArtifact(
        name="mv_ligand_pair_alignment_raw_pair_idx",
        sql="CREATE INDEX mv_ligand_pair_alignment_raw_pair_idx ON mv_ligand_pair_alignment_raw (pdblid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_pair_alignment_raw_pair_idx",
        dependencies=["mv_ligand_pair_alignment_raw"],
        description="Pair index (alignment raw)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_pair_alignment_raw_smmclid_idx": SQLArtifact(
        name="mv_ligand_pair_alignment_raw_smmclid_idx",
        sql="CREATE INDEX mv_ligand_pair_alignment_raw_smmclid_idx ON mv_ligand_pair_alignment_raw (smmclid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_pair_alignment_raw_smmclid_idx",
        dependencies=["mv_ligand_pair_alignment_raw"],
        description="Index on smmclid (alignment raw)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_pair_alignment_raw_identity_idx": SQLArtifact(
        name="mv_ligand_pair_alignment_raw_identity_idx",
        sql="CREATE INDEX mv_ligand_pair_alignment_raw_identity_idx ON mv_ligand_pair_alignment_raw (percent_identity DESC);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_pair_alignment_raw_identity_idx",
        dependencies=["mv_ligand_pair_alignment_raw"],
        description="Index on percent_identity DESC (alignment raw)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_supported_pair_best_identity_interval_idx": SQLArtifact(
        name="mv_ligand_supported_pair_best_identity_interval_idx",
        sql="CREATE INDEX mv_ligand_supported_pair_best_identity_interval_idx ON mv_ligand_supported_pair_best_identity (identity_interval);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_supported_pair_best_identity_interval_idx",
        dependencies=["mv_ligand_supported_pair_best_identity"],
        description="Index on identity_interval (best identity MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_supported_pair_best_identity_smmclid_idx": SQLArtifact(
        name="mv_ligand_supported_pair_best_identity_smmclid_idx",
        sql="CREATE INDEX mv_ligand_supported_pair_best_identity_smmclid_idx ON mv_ligand_supported_pair_best_identity (best_smmclid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_supported_pair_best_identity_smmclid_idx",
        dependencies=["mv_ligand_supported_pair_best_identity"],
        description="Index on best_smmclid (best identity MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_supported_pair_identity_multibucket_interval_idx": SQLArtifact(
        name="mv_ligand_supported_pair_identity_multibucket_interval_idx",
        sql="CREATE INDEX mv_ligand_supported_pair_identity_multibucket_interval_idx ON mv_ligand_supported_pair_identity_multibucket (identity_interval);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_supported_pair_identity_multibucket_interval_idx",
        dependencies=["mv_ligand_supported_pair_identity_multibucket"],
        description="Index on identity_interval (multibucket MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_supported_pair_identity_multibucket_pair_idx": SQLArtifact(
        name="mv_ligand_supported_pair_identity_multibucket_pair_idx",
        sql="CREATE INDEX mv_ligand_supported_pair_identity_multibucket_pair_idx ON mv_ligand_supported_pair_identity_multibucket (pdblid, clid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_supported_pair_identity_multibucket_pair_idx",
        dependencies=["mv_ligand_supported_pair_identity_multibucket"],
        description="Pair index (multibucket MV)",
        artifact_type=ArtifactType.INDEX,
    ),
    "mv_ligand_supported_pair_identity_multibucket_smmclid_idx": SQLArtifact(
        name="mv_ligand_supported_pair_identity_multibucket_smmclid_idx",
        sql="CREATE INDEX mv_ligand_supported_pair_identity_multibucket_smmclid_idx ON mv_ligand_supported_pair_identity_multibucket (smmclid);",
        drop_sql="DROP INDEX IF EXISTS mv_ligand_supported_pair_identity_multibucket_smmclid_idx",
        dependencies=["mv_ligand_supported_pair_identity_multibucket"],
        description="Index on smmclid (multibucket MV)",
        artifact_type=ArtifactType.INDEX,
    ),
}


# ======================================================================================
# Dependency Management
# ======================================================================================


class DependencyManager:
    def __init__(self, artifacts: Dict[str, SQLArtifact]):
        self.artifacts = artifacts
        self.graph = self._build_dependency_graph()

    def _build_dependency_graph(self) -> Dict[str, Set[str]]:
        """Build dependency graph: artifact -> set of artifacts that depend on it."""
        graph = {name: set() for name in self.artifacts.keys()}

        for name, artifact in self.artifacts.items():
            for dep in artifact.dependencies:
                if dep in graph:
                    graph[dep].add(name)

        return graph

    def get_creation_order(self, artifact_names: List[str]) -> List[str]:
        """Return artifacts in dependency order for creation (dependencies first)."""
        visited = set()
        result = []

        def visit(name: str):
            if name in visited or name not in self.artifacts:
                return
            visited.add(name)

            # Visit dependencies first
            for dep in self.artifacts[name].dependencies:
                visit(dep)

            result.append(name)

        for name in artifact_names:
            visit(name)

        return result

    def get_deletion_order(self, artifact_names: List[str]) -> List[str]:
        """Return artifacts in reverse dependency order for deletion (dependents first)."""
        visited = set()
        result = []

        def visit(name: str):
            if name in visited or name not in self.artifacts:
                return
            visited.add(name)

            # Visit dependents first
            for dependent in self.graph.get(name, set()):
                if dependent in artifact_names:  # Only include requested artifacts
                    visit(dependent)

            result.append(name)

        for name in artifact_names:
            visit(name)

        return result

    def get_dependents(self, artifact_name: str) -> Set[str]:
        """Get all artifacts that depend on this one."""
        return self.graph.get(artifact_name, set())

    def validate_dependencies(self, artifact_name: str) -> Tuple[bool, List[str]]:
        """Check if all dependencies exist. Returns (valid, missing_deps)."""
        if artifact_name not in self.artifacts:
            return False, [f"Artifact {artifact_name} not found"]

        missing = []
        for dep in self.artifacts[artifact_name].dependencies:
            if dep not in self.artifacts:
                missing.append(dep)

        return len(missing) == 0, missing


# ======================================================================================
# Database Operations
# ======================================================================================


class DatabaseExecutor:
    def __init__(self, db: UnifiedPostgresConnector):
        self.db = db

    def artifact_exists(self, artifact: SQLArtifact) -> bool:
        """Check if artifact exists in database."""
        try:
            if artifact.artifact_type == ArtifactType.FUNCTION:
                query = sql.SQL("SELECT 1 FROM pg_proc WHERE proname = %s")
            elif artifact.artifact_type == ArtifactType.MATERIALIZED_VIEW:
                query = sql.SQL("SELECT 1 FROM pg_matviews WHERE matviewname = %s")
            elif artifact.artifact_type == ArtifactType.VIEW:
                query = sql.SQL("SELECT 1 FROM pg_views WHERE viewname = %s")
            elif artifact.artifact_type == ArtifactType.INDEX:
                query = sql.SQL("SELECT 1 FROM pg_indexes WHERE indexname = %s")
            else:
                return False

            with self.db.raw_connection.cursor() as cur:
                cur.execute(query, (artifact.name,))
                return cur.fetchone() is not None
        except Exception as e:
            logger.warning(f"Error checking existence of {artifact.name}: {e}")
            return False

    def create_artifact(self, artifact: SQLArtifact, force: bool = False) -> bool:
        """Create a single artifact."""
        try:
            if force and self.artifact_exists(artifact):
                logger.info(f"Force dropping existing {artifact.name}")
                self.drop_artifact(artifact, cascade=True)

            logger.info(f"Creating {artifact.artifact_type.value} {artifact.name}")
            with self.db.raw_connection.cursor() as cur:
                cur.execute(artifact.sql)
                self.db.raw_connection.commit()

            logger.success(f"Successfully created {artifact.name}")
            return True

        except Exception as e:
            logger.error(f"Failed to create {artifact.name}: {e}")
            self.db.raw_connection.rollback()
            return False

    def drop_artifact(self, artifact: SQLArtifact, cascade: bool = False) -> bool:
        """Drop a single artifact."""
        try:
            drop_sql = artifact.drop_sql
            if cascade and "CASCADE" not in drop_sql.upper():
                drop_sql = drop_sql.replace(";", " CASCADE;")

            logger.info(f"Dropping {artifact.artifact_type.value} {artifact.name}")
            with self.db.raw_connection.cursor() as cur:
                cur.execute(drop_sql)
                self.db.raw_connection.commit()

            logger.success(f"Successfully dropped {artifact.name}")
            return True

        except Exception as e:
            logger.error(f"Failed to drop {artifact.name}: {e}")
            self.db.raw_connection.rollback()
            return False

    def refresh_materialized_view(
        self, artifact: SQLArtifact, concurrent: bool = False
    ) -> bool:
        """Refresh a materialized view."""
        if not artifact.can_refresh():
            logger.error(f"{artifact.name} is not a materialized view")
            return False

        try:
            # Determine refresh statement without shadowing psycopg2.sql module
            if concurrent and artifact.supports_concurrent_refresh():
                refresh_stmt = artifact.refresh_concurrent_sql
                logger.info(f"Refreshing {artifact.name} CONCURRENTLY")
            else:
                refresh_stmt = artifact.refresh_sql
                logger.info(f"Refreshing {artifact.name}")

            if not refresh_stmt:
                logger.error(f"No refresh SQL defined for {artifact.name}")
                return False

            with self.db.raw_connection.cursor() as cur:
                cur.execute(refresh_stmt)
                self.db.raw_connection.commit()

            logger.success(f"Successfully refreshed {artifact.name}")
            return True

        except Exception as e:
            logger.error(f"Failed to refresh {artifact.name}: {e}")
            self.db.raw_connection.rollback()
            return False

    def get_artifact_status(self, artifact: SQLArtifact) -> Dict:
        """Get detailed status information for an artifact."""
        status = {
            "name": artifact.name,
            "type": artifact.artifact_type.value,
            "exists": False,
            "size": None,
            "row_count": None,
        }

        try:
            status["exists"] = self.artifact_exists(artifact)

            if (
                status["exists"]
                and artifact.artifact_type == ArtifactType.MATERIALIZED_VIEW
            ):
                # Get size and row count for materialized views
                with self.db.raw_connection.cursor() as cur:
                    # Use sql.Identifier for safe table name interpolation
                    query = sql.SQL("SELECT COUNT(*) FROM {}").format(
                        sql.Identifier(artifact.name)
                    )
                    cur.execute(query)
                    fetch = cur.fetchone()
                    status["row_count"] = fetch[0] if fetch else None

                    # Use parameterized query for function argument
                    cur.execute(
                        "SELECT pg_size_pretty(pg_total_relation_size(%s))",
                        (artifact.name,),
                    )
                    fetch_size = cur.fetchone()
                    status["size"] = fetch_size[0] if fetch_size else None

        except Exception as e:
            logger.warning(f"Error getting status for {artifact.name}: {e}")

        return status


# ======================================================================================
# Connection Management (Following Repository Patterns)
# ======================================================================================


def _connect_activity_db(
    user: Optional[str] = None,
    password: Optional[str] = None,
    host: Optional[str] = None,
    port: Optional[int] = None,
    db_name: Optional[str] = None,
) -> UnifiedPostgresConnector:
    """Build ActivityDB connection using AppConfig.

    SSH tunneling is managed at the CLI layer (--manage-tunnel, --tunnel-mode).
    Configure host=127.0.0.1 and the tunnel port in appconfig.toml when using CLI-managed tunnel.
    """
    conn_params = AppConfig.set_activity_db_database_connection_parameters(
        cli_options={
            k: v
            for k, v in [
                ("name", db_name),
                ("user", user),
                ("password", password),
                ("host", host),
                ("port", port),
            ]
            if v is not None
        }
    )

    logger.info(
        f"Connecting to ActivityDB host={conn_params.host} db={conn_params.name}"
    )

    return build_unified_connector(
        db_name=conn_params.name,
        user=conn_params.user,
        password=conn_params.password,
        host=conn_params.host,
        port=conn_params.port,
        auto_connect=True,
    )


# ======================================================================================
# Typer CLI Application
# ======================================================================================

app = typer.Typer(
    name="sql-artifacts-management",
    help="Manage SQL artifacts for ligand mapping analysis",
    add_completion=False,
)


# Common database connection options
def _add_db_options(f):
    """Decorator to add common database connection options."""
    import functools

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    # Add type hints manually since decorators can't modify them easily
    wrapper.__annotations__ = f.__annotations__.copy()
    wrapper.__annotations__.update(
        {
            "user": Optional[str],
            "pwd": Optional[str],
            "host": Optional[str],
            "port": Optional[int],
            "database_name": Optional[str],
        }
    )

    return wrapper


@app.command()
def create(
    artifact_name: str = typer.Argument(..., help="Name of artifact to create"),
    force: bool = typer.Option(False, "--force", help="Drop existing artifact first"),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Create a single SQL artifact."""
    if artifact_name not in ARTIFACTS:
        logger.error(f"Artifact '{artifact_name}' not found")
        typer.echo(f"Available artifacts: {', '.join(sorted(ARTIFACTS.keys()))}")
        raise typer.Exit(1)

    artifact = ARTIFACTS[artifact_name]
    dep_mgr = DependencyManager(ARTIFACTS)

    # Check dependencies
    valid, missing = dep_mgr.validate_dependencies(artifact_name)
    if not valid:
        logger.error(f"Missing dependencies: {missing}")
        raise typer.Exit(1)

    # Get creation order (includes dependencies)
    creation_order = dep_mgr.get_creation_order([artifact_name])

    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        success_count = 0
        for name in creation_order:
            if (
                name == artifact_name or force
            ):  # Only create dependencies if they don't exist
                art = ARTIFACTS[name]
                if executor.create_artifact(
                    art, force=(name == artifact_name and force)
                ):
                    success_count += 1
                else:
                    logger.error(f"Failed to create {name}, stopping")
                    raise typer.Exit(1)
            elif not executor.artifact_exists(ARTIFACTS[name]):
                logger.error(
                    f"Dependency {name} does not exist. Use --force to create it."
                )
                raise typer.Exit(1)

        logger.success(f"Successfully created {success_count} artifact(s)")

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command()
def drop(
    artifact_name: str = typer.Argument(..., help="Name of artifact to drop"),
    cascade: bool = typer.Option(False, "--cascade", help="Use CASCADE option"),
    force: bool = typer.Option(False, "--force", help="Skip confirmation"),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Drop a single SQL artifact."""
    if artifact_name not in ARTIFACTS:
        logger.error(f"Artifact '{artifact_name}' not found")
        raise typer.Exit(1)

    artifact = ARTIFACTS[artifact_name]
    dep_mgr = DependencyManager(ARTIFACTS)

    # Check for dependents
    dependents = dep_mgr.get_dependents(artifact_name)
    if dependents and not cascade:
        logger.warning(f"Artifact {artifact_name} has dependents: {dependents}")
        logger.info("Use --cascade to drop dependents as well")
        if not force:
            if not typer.confirm("Continue anyway?"):
                raise typer.Exit(0)

    # Confirmation
    if not force:
        action = "CASCADE drop" if cascade else "drop"
        if not typer.confirm(f"Are you sure you want to {action} {artifact_name}?"):
            raise typer.Exit(0)

    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        if executor.drop_artifact(artifact, cascade):
            logger.success(f"Successfully dropped {artifact_name}")
        else:
            raise typer.Exit(1)

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command()
def refresh(
    artifact_name: str = typer.Argument(
        ..., help="Name of materialized view to refresh"
    ),
    concurrent: bool = typer.Option(
        False, "--concurrent", help="Use REFRESH CONCURRENTLY"
    ),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Refresh a materialized view."""
    if artifact_name not in ARTIFACTS:
        logger.error(f"Artifact '{artifact_name}' not found")
        raise typer.Exit(1)

    artifact = ARTIFACTS[artifact_name]
    if not artifact.can_refresh():
        logger.error(f"{artifact_name} is not a materialized view")
        raise typer.Exit(1)

    if concurrent and not artifact.supports_concurrent_refresh():
        logger.warning(
            f"{artifact_name} does not support concurrent refresh, using regular refresh"
        )
        concurrent = False

    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        if executor.refresh_materialized_view(artifact, concurrent):
            logger.success(f"Successfully refreshed {artifact_name}")
        else:
            raise typer.Exit(1)

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command()
def status(
    artifact_name: Optional[str] = typer.Argument(
        None, help="Name of specific artifact (optional)"
    ),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Show status of artifacts."""
    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        artifacts_to_check = (
            [artifact_name] if artifact_name else list(ARTIFACTS.keys())
        )

        if artifact_name and artifact_name not in ARTIFACTS:
            logger.error(f"Artifact '{artifact_name}' not found")
            raise typer.Exit(1)

        typer.echo(
            f"{'Name':<50} {'Type':<20} {'Exists':<10} {'Size':<15} {'Rows':<15}"
        )
        typer.echo("-" * 120)

        for name in sorted(artifacts_to_check):
            artifact = ARTIFACTS[name]
            status_info = executor.get_artifact_status(artifact)

            exists_str = "✓" if status_info["exists"] else "✗"
            size_str = status_info["size"] or "N/A"
            rows_str = (
                f"{status_info['row_count']:,}"
                if status_info["row_count"] is not None
                else "N/A"
            )

            typer.echo(
                f"{name:<50} {status_info['type']:<20} {exists_str:<10} {size_str:<15} {rows_str:<15}"
            )

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command(name="list")
def list_artifacts(
    artifact_type: Optional[ArtifactType] = typer.Option(
        None, "--type", help="Filter by artifact type"
    ),
):
    """List all available artifacts."""
    artifacts_to_show = ARTIFACTS

    if artifact_type:
        artifacts_to_show = {
            name: artifact
            for name, artifact in ARTIFACTS.items()
            if artifact.artifact_type == artifact_type
        }

    typer.echo(f"{'Name':<50} {'Type':<20} {'Dependencies':<30} {'Description'}")
    typer.echo("-" * 150)

    for name in sorted(artifacts_to_show.keys()):
        artifact = artifacts_to_show[name]
        deps_str = ", ".join(artifact.dependencies) if artifact.dependencies else "None"
        if len(deps_str) > 28:
            deps_str = deps_str[:25] + "..."

        desc_str = artifact.description
        if len(desc_str) > 48:
            desc_str = desc_str[:45] + "..."

        typer.echo(
            f"{name:<50} {artifact.artifact_type.value:<20} {deps_str:<30} {desc_str}"
        )

    typer.echo(f"\nTotal: {len(artifacts_to_show)} artifacts")


@app.command()
def create_all(
    artifact_type: Optional[ArtifactType] = typer.Option(
        None, "--type", help="Create only specific type"
    ),
    force: bool = typer.Option(False, "--force", help="Drop existing artifacts first"),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show what would be created without executing"
    ),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Create all artifacts in dependency order."""
    artifacts_to_create = list(ARTIFACTS.keys())

    if artifact_type:
        artifacts_to_create = [
            name
            for name, artifact in ARTIFACTS.items()
            if artifact.artifact_type == artifact_type
        ]

    dep_mgr = DependencyManager(ARTIFACTS)
    creation_order = dep_mgr.get_creation_order(artifacts_to_create)

    if dry_run:
        typer.echo("Would create artifacts in this order:")
        for i, name in enumerate(creation_order, 1):
            artifact = ARTIFACTS[name]
            typer.echo(f"{i:2d}. {name} ({artifact.artifact_type.value})")
        return

    if not force:
        typer.echo(f"Will create {len(creation_order)} artifacts in dependency order")
        if not typer.confirm("Continue?"):
            raise typer.Exit(0)

    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        success_count = 0
        for name in creation_order:
            artifact = ARTIFACTS[name]
            if executor.create_artifact(artifact, force):
                success_count += 1
            else:
                logger.error(f"Failed to create {name}, stopping")
                raise typer.Exit(1)

        logger.success(f"Successfully created {success_count} artifacts")

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command()
def drop_all(
    artifact_type: Optional[ArtifactType] = typer.Option(
        None, "--type", help="Drop only specific type"
    ),
    cascade: bool = typer.Option(False, "--cascade", help="Use CASCADE option"),
    force: bool = typer.Option(False, "--force", help="Skip confirmation"),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Show what would be dropped without executing"
    ),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Drop all artifacts in reverse dependency order."""
    artifacts_to_drop = list(ARTIFACTS.keys())

    if artifact_type:
        artifacts_to_drop = [
            name
            for name, artifact in ARTIFACTS.items()
            if artifact.artifact_type == artifact_type
        ]

    dep_mgr = DependencyManager(ARTIFACTS)
    deletion_order = dep_mgr.get_deletion_order(artifacts_to_drop)

    if dry_run:
        typer.echo("Would drop artifacts in this order:")
        for i, name in enumerate(deletion_order, 1):
            artifact = ARTIFACTS[name]
            typer.echo(f"{i:2d}. {name} ({artifact.artifact_type.value})")
        return

    # Extra confirmation for destructive operation
    if not force:
        action = "CASCADE drop" if cascade else "drop"
        typer.echo(
            f"Will {action} {len(deletion_order)} artifacts in reverse dependency order"
        )
        typer.echo("This is a destructive operation!")
        if not typer.confirm("Are you absolutely sure?"):
            raise typer.Exit(0)

    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        success_count = 0
        for name in deletion_order:
            artifact = ARTIFACTS[name]
            if executor.drop_artifact(artifact, cascade):
                success_count += 1
            else:
                logger.warning(f"Failed to drop {name}, continuing")

        logger.success(f"Successfully dropped {success_count} artifacts")

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command()
def refresh_all_mvs(
    concurrent: bool = typer.Option(
        False, "--concurrent", help="Use REFRESH CONCURRENTLY where supported"
    ),
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
):
    """Refresh all materialized views in dependency order."""
    mv_artifacts = [
        (name, artifact)
        for name, artifact in ARTIFACTS.items()
        if artifact.artifact_type == ArtifactType.MATERIALIZED_VIEW
    ]

    dep_mgr = DependencyManager(ARTIFACTS)
    mv_names = [name for name, _ in mv_artifacts]
    refresh_order = dep_mgr.get_creation_order(mv_names)  # Dependencies first

    typer.echo(f"Will refresh {len(refresh_order)} materialized views")

    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        success_count = 0
        for name in refresh_order:
            artifact = ARTIFACTS[name]
            use_concurrent = concurrent and artifact.supports_concurrent_refresh()

            if executor.refresh_materialized_view(artifact, use_concurrent):
                success_count += 1
            else:
                logger.error(f"Failed to refresh {name}")

        logger.success(f"Successfully refreshed {success_count} materialized views")

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


@app.command()
def dependencies(
    artifact_name: str = typer.Argument(..., help="Name of artifact to analyze"),
):
    """Show dependency information for an artifact."""
    if artifact_name not in ARTIFACTS:
        logger.error(f"Artifact '{artifact_name}' not found")
        raise typer.Exit(1)

    artifact = ARTIFACTS[artifact_name]
    dep_mgr = DependencyManager(ARTIFACTS)

    typer.echo(f"\nArtifact: {artifact_name}")
    typer.echo(f"Type: {artifact.artifact_type.value}")
    typer.echo(f"Description: {artifact.description}")

    # Direct dependencies
    if artifact.dependencies:
        typer.echo(f"\nDirect dependencies ({len(artifact.dependencies)}):")
        for dep in artifact.dependencies:
            if dep in ARTIFACTS:
                dep_type = ARTIFACTS[dep].artifact_type.value
                typer.echo(f"  - {dep} ({dep_type})")
            else:
                typer.echo(f"  - {dep} (NOT FOUND)")
    else:
        typer.echo("\nNo direct dependencies")

    # Things that depend on this artifact
    dependents = dep_mgr.get_dependents(artifact_name)
    if dependents:
        typer.echo(f"\nDependents ({len(dependents)}):")
        for dep in sorted(dependents):
            dep_type = ARTIFACTS[dep].artifact_type.value
            typer.echo(f"  - {dep} ({dep_type})")
    else:
        typer.echo("\nNo dependents")

    # Full creation order if this artifact were to be created
    creation_order = dep_mgr.get_creation_order([artifact_name])
    if len(creation_order) > 1:
        typer.echo(f"\nFull creation order ({len(creation_order)}):")
        for i, name in enumerate(creation_order, 1):
            marker = "👉" if name == artifact_name else "  "
            typer.echo(f"{marker} {i}. {name}")


@app.command()
def validate(
    user: Optional[str] = typer.Option(None, help="ActivityDB user override"),
    password: Optional[str] = typer.Option(None, "--password", help="ActivityDB password override"),
    host: Optional[str] = typer.Option(None, help="ActivityDB host override"),
    port: Optional[int] = typer.Option(None, help="ActivityDB port override"),
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="ActivityDB database name override"
    ),
) -> None:
    """Validate all artifacts and their dependencies."""
    typer.echo("Validating artifact definitions and database state...")

    dep_mgr = DependencyManager(ARTIFACTS)
    db = _connect_activity_db(user, password, host, port, db_name)
    executor = DatabaseExecutor(db)

    try:
        issues = []

        # Check for missing dependencies in definitions
        for name, artifact in ARTIFACTS.items():
            valid, missing = dep_mgr.validate_dependencies(name)
            if not valid:
                issues.append(
                    f"❌ {name}: Missing dependencies in definitions: {missing}"
                )

        # Check database state
        exists_count = 0
        missing_count = 0

        for name, artifact in ARTIFACTS.items():
            exists = executor.artifact_exists(artifact)
            if exists:
                exists_count += 1
                typer.echo(f"✓ {name}")
            else:
                missing_count += 1
                typer.echo(f"✗ {name}")

        typer.echo(f"\nSummary:")
        typer.echo(f"  Total artifacts: {len(ARTIFACTS)}")
        typer.echo(f"  Existing in DB: {exists_count}")
        typer.echo(f"  Missing from DB: {missing_count}")

        if issues:
            typer.echo(f"\nIssues found:")
            for issue in issues:
                typer.echo(f"  {issue}")
        else:
            typer.echo("\nNo definition issues found ✓")

    finally:
        if hasattr(db, "connection") and db.connection:
            db.connection.close()


if __name__ == "__main__":
    # Entry point
    app()
