"""
Ligand Mapping SQL Query Constants
===================================

Centralized SQL queries for ligand mapping analysis levels 1-4.
Source: reports/ligand_mapping_sql_artifacts_optimized.sql

This module contains all read-only analytical queries used by the ligand mapping
level apps. Queries are organized by level and analysis type.

Note: These queries assume the existence of materialized views and supporting
objects created by ligand_mapping_ddl_artifacts_optimized.sql.
"""

# =====================================================================================
# SECTION 2: LEVEL 1 — STRUCTURE-ONLY LIGAND MAPPING (BASELINE)
# =====================================================================================

LEVEL1_BASELINE_SQL = """
-- 2.1 Comprehensive Baseline Metrics (OPTIMIZED: Single pass)
WITH baseline_stats AS (
    SELECT
        COUNT(*) AS total_mapping_pairs,
        COUNT(DISTINCT pdblid) AS distinct_mapped_pdblid,
        COUNT(DISTINCT clid) AS distinct_mapped_clid
    FROM pdb_ligand_chembl_ligand
),
universe_stats AS (
    SELECT
        (SELECT COUNT(DISTINCT pdblid) FROM pdb_ligand) AS pdb_universe,
        (SELECT COUNT(DISTINCT clid) FROM chembl_ligand) AS chembl_universe
)
SELECT
    b.total_mapping_pairs,
    b.distinct_mapped_pdblid,
    b.distinct_mapped_clid,
    u.pdb_universe,
    u.chembl_universe,
    ROUND((b.distinct_mapped_pdblid * 100.0 / u.pdb_universe), 2) AS pdb_ligand_coverage_pct,
    ROUND((b.distinct_mapped_clid * 100.0 / u.chembl_universe), 2) AS chembl_ligand_coverage_pct,
    u.pdb_universe - b.distinct_mapped_pdblid AS unmapped_pdblid_count,
    u.chembl_universe - b.distinct_mapped_clid AS unmapped_clid_count
FROM baseline_stats b, universe_stats u
"""

LEVEL1_SMMCLID_SQL = """
-- 2.2 smmclid Distribution with Exclusivity Analysis (OPTIMIZED: Window functions)
WITH pdblid_smmclid_counts AS (
    SELECT
        pdblid,
        COUNT(DISTINCT smmclid) as pdblid_smmclid_count
    FROM pdb_ligand_chembl_ligand
    GROUP BY pdblid
),
clid_smmclid_counts AS (
    SELECT
        clid,
        COUNT(DISTINCT smmclid) as clid_smmclid_count
    FROM pdb_ligand_chembl_ligand
    GROUP BY clid
),
ligand_smmclid_analysis AS (
    SELECT
        plcl.smmclid,
        plcl.pdblid,
        plcl.clid,
        psc.pdblid_smmclid_count,
        csc.clid_smmclid_count
    FROM pdb_ligand_chembl_ligand plcl
    JOIN pdblid_smmclid_counts psc ON plcl.pdblid = psc.pdblid
    JOIN clid_smmclid_counts csc ON plcl.clid = csc.clid
)
SELECT
    smmclid,
    COUNT(*) AS pairs,
    COUNT(DISTINCT pdblid) AS distinct_pdblid,
    COUNT(DISTINCT clid) AS distinct_clid,
    COUNT(DISTINCT pdblid) FILTER (WHERE pdblid_smmclid_count = 1) AS pdblid_exclusives,
    COUNT(DISTINCT clid) FILTER (WHERE clid_smmclid_count = 1) AS clid_exclusives,
    ROUND(COUNT(DISTINCT pdblid) * 100.0 / SUM(COUNT(DISTINCT pdblid)) OVER(), 2) AS pdblid_pct,
    ROUND(COUNT(DISTINCT clid) * 100.0 / SUM(COUNT(DISTINCT clid)) OVER(), 2) AS clid_pct
FROM ligand_smmclid_analysis
GROUP BY smmclid
ORDER BY smmclid
"""

LEVEL1_FANOUT_SQL = """
-- 2.3 Multi-Level Occurrence and Fan-Out Analysis (OPTIMIZED: Combined)
WITH pdblid_analysis AS (
    SELECT
        pdblid,
        COUNT(DISTINCT smmclid) AS smmclid_count,
        COUNT(DISTINCT clid) AS clid_partners,
        STRING_AGG(smmclid::text, ', ') AS smmclid_list
    FROM pdb_ligand_chembl_ligand
    GROUP BY pdblid
),
fanout_by_smmclid AS (
    SELECT
        smmclid,
        pdblid,
        COUNT(DISTINCT clid) AS clid_partners
    FROM pdb_ligand_chembl_ligand
    GROUP BY smmclid, pdblid
)
SELECT
    'Multi-Level Occurrence' AS analysis_type,
    CASE
        WHEN smmclid_count = 1 THEN 'Single smmclid'
        WHEN smmclid_count = 2 THEN 'Two smmclid'
        ELSE 'Three+ smmclid'
    END AS category,
    COUNT(*) AS pdblid_count,
    ROUND(COUNT(*) * 100.0 / SUM(COUNT(*)) OVER(), 2) AS percentage,
    ROUND(AVG(clid_partners), 2) AS avg_clid_partners,
    MAX(clid_partners) AS max_clid_partners
FROM pdblid_analysis
GROUP BY CASE
    WHEN smmclid_count = 1 THEN 'Single smmclid'
    WHEN smmclid_count = 2 THEN 'Two smmclid'
    ELSE 'Three+ smmclid'
END
UNION ALL
SELECT
    'Fan-Out by smmclid' AS analysis_type,
    'smmclid=' || smmclid AS category,
    COUNT(DISTINCT pdblid) AS pdblid_count,
    NULL AS percentage,
    ROUND(AVG(clid_partners), 2) AS avg_clid_partners,
    MAX(clid_partners) AS max_clid_partners
FROM fanout_by_smmclid
GROUP BY smmclid
ORDER BY analysis_type, category
"""

LEVEL1_UNMAPPED_SQL = """
-- 2.4 Unmapped Ligand Characterization (OPTIMIZED)
WITH unmapped_analysis AS (
    SELECT
        pl.pdblid,
        pl.ligand_id,
        pl.standard_inchi_key,
        LENGTH(pl.ligand_id) AS ligand_id_length,
        CASE
            WHEN pl.standard_inchi_key IS NOT NULL THEN 'Has InChIKey'
            ELSE 'No InChIKey'
        END AS inchi_status
    FROM pdb_ligand pl
    WHERE NOT EXISTS (
        SELECT 1 FROM pdb_ligand_chembl_ligand plcl
        WHERE plcl.pdblid = pl.pdblid
    )
)
SELECT
    COUNT(*) AS unmapped_pdblid_count,
    COUNT(*) FILTER (WHERE inchi_status = 'Has InChIKey') AS unmapped_with_inchikey,
    ROUND(AVG(ligand_id_length), 1) AS avg_ligand_id_length,
    MODE() WITHIN GROUP (ORDER BY ligand_id_length) AS most_common_ligand_id_length,
    COUNT(*) FILTER (WHERE ligand_id_length = 3) AS three_char_ligands,
    COUNT(*) FILTER (WHERE ligand_id_length > 3) AS longer_ligands
FROM unmapped_analysis
"""

# =====================================================================================
# SECTION 3: LEVEL 2 — QUADRUPLET UNIVERSE (UNPARTITIONED BIOLOGICAL MAPPING)
# =====================================================================================

LEVEL2_QUADRUPLETS_SQL = """
-- 3.1 Core Quadruplet Counts with Expansion Factors (OPTIMIZED: Single query)
WITH combined_stats AS (
    SELECT
        'STRICT' AS variant,
        COUNT(*) AS total_quadruplets,
        COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid)) AS distinct_quadruplets,
        COUNT(DISTINCT (pdblid, clid, chembl_tid)) AS distinct_triplets,
        COUNT(DISTINCT (pdblid, clid)) AS distinct_pairs
    FROM mv_ligand_sequence_target_fact
    UNION ALL
    SELECT
        'RELAXED' AS variant,
        COUNT(*) AS total_quadruplets,
        COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid)) AS distinct_quadruplets,
        COUNT(DISTINCT (pdblid, clid, chembl_tid)) AS distinct_triplets,
        COUNT(DISTINCT (pdblid, clid)) AS distinct_pairs
    FROM mv_ligand_sequence_target_fact_nohitcov
),
with_ratios AS (
    SELECT
        variant,
        total_quadruplets,
        distinct_quadruplets,
        distinct_triplets,
        distinct_pairs,
        ROUND(distinct_quadruplets::numeric / NULLIF(distinct_pairs, 0), 2) AS quad_per_pair,
        ROUND(distinct_triplets::numeric / NULLIF(distinct_pairs, 0), 2) AS triplet_per_pair
    FROM combined_stats
)
SELECT
    w.*,
    CASE
        WHEN variant = 'RELAXED' THEN
            ROUND(distinct_quadruplets::numeric /
                  NULLIF(LAG(distinct_quadruplets) OVER (ORDER BY variant), 0), 2)
        ELSE NULL
    END AS expansion_factor_quad,
    CASE
        WHEN variant = 'RELAXED' THEN
            ROUND(distinct_triplets::numeric /
                  NULLIF(LAG(distinct_triplets) OVER (ORDER BY variant), 0), 2)
        ELSE NULL
    END AS expansion_factor_triplet,
    CASE
        WHEN variant = 'RELAXED' THEN
            ROUND(distinct_pairs::numeric /
                  NULLIF(LAG(distinct_pairs) OVER (ORDER BY variant), 0), 2)
        ELSE NULL
    END AS expansion_factor_pair
FROM with_ratios w
ORDER BY variant
"""

LEVEL2_EXPANSION_SQL = """
-- Level 2 expansion analysis is included in LEVEL2_QUADRUPLETS_SQL above
-- This placeholder remains for consistency with the strategy document
SELECT 'See LEVEL2_QUADRUPLETS_SQL for expansion factors' AS note
"""

# =====================================================================================
# SECTION 4: LEVEL 3 — QUADRUPLET PARTITION BY SEQUENCE TYPES & MAPPING LEVEL
# =====================================================================================

LEVEL3_PARTITIONS_SQL = """
-- 4.1 Complete Partitioned Breakdown (OPTIMIZED: Single pass with GROUPING SETS)
WITH combined_base AS (
    SELECT 'STRICT' AS variant, pdblid, query_sid, clid, chembl_tid, smmclid, pdb_stid, chembl_type
    FROM mv_ligand_sequence_target_fact
    UNION ALL
    SELECT 'RELAXED' AS variant, pdblid, query_sid, clid, chembl_tid, smmclid, pdb_stid, chembl_type
    FROM mv_ligand_sequence_target_fact_nohitcov
)
SELECT
    variant,
    smmclid,
    pdb_stid,
    chembl_type,
    COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid)) AS distinct_quadruplets,
    COUNT(DISTINCT (pdblid, clid, chembl_tid)) AS distinct_triplets,
    COUNT(DISTINCT (pdblid, clid)) AS distinct_pairs,
    COUNT(DISTINCT pdblid) AS distinct_pdblid,
    COUNT(DISTINCT clid) AS distinct_clid,
    -- Add ratios for analysis
    ROUND(COUNT(DISTINCT (pdblid, query_sid, clid, chembl_tid))::numeric /
          NULLIF(COUNT(DISTINCT (pdblid, clid)), 0), 2) AS quad_per_pair_ratio,
    ROUND(COUNT(DISTINCT pdblid)::numeric /
          NULLIF(COUNT(DISTINCT clid), 0), 2) AS pdblid_per_clid_ratio
FROM combined_base
GROUP BY variant, smmclid, pdb_stid, chembl_type
ORDER BY
    variant,
    COALESCE(smmclid, 999),
    COALESCE(pdb_stid, 999),
    COALESCE(chembl_type, 'ZZZZZ')
"""

LEVEL3_EXPANSION_DISTINCT_SQL = """
-- 4.2 Distinct Ligand/Pair Expansion by smmclid (OPTIMIZED: With percentage gains)
WITH expansion_stats AS (
    SELECT
        'STRICT' as variant,
        smmclid,
        COUNT(DISTINCT pdblid) as distinct_pdblid,
        COUNT(DISTINCT (pdblid, clid)) as distinct_pairs
    FROM mv_ligand_sequence_target_fact
    GROUP BY smmclid
    UNION ALL
    SELECT
        'RELAXED' as variant,
        smmclid,
        COUNT(DISTINCT pdblid) as distinct_pdblid,
        COUNT(DISTINCT (pdblid, clid)) as distinct_pairs
    FROM mv_ligand_sequence_target_fact_nohitcov
    GROUP BY smmclid
),
with_expansion AS (
    SELECT
        smmclid,
        MAX(CASE WHEN variant = 'STRICT' THEN distinct_pdblid END) AS strict_pdblid,
        MAX(CASE WHEN variant = 'RELAXED' THEN distinct_pdblid END) AS relaxed_pdblid,
        MAX(CASE WHEN variant = 'STRICT' THEN distinct_pairs END) AS strict_pairs,
        MAX(CASE WHEN variant = 'RELAXED' THEN distinct_pairs END) AS relaxed_pairs
    FROM expansion_stats
    GROUP BY smmclid
)
SELECT
    smmclid,
    strict_pdblid,
    relaxed_pdblid,
    ROUND(relaxed_pdblid::numeric / NULLIF(strict_pdblid, 0), 2) AS pdblid_expansion,
    strict_pairs,
    relaxed_pairs,
    ROUND(relaxed_pairs::numeric / NULLIF(strict_pairs, 0), 2) AS pair_expansion,
    relaxed_pdblid - strict_pdblid AS pdblid_gain,
    relaxed_pairs - strict_pairs AS pair_gain
FROM with_expansion
ORDER BY smmclid
"""

# =====================================================================================
# SECTION 5: LEVEL 4 — IDENTITY-QUALIFIED SURVIVAL (RETENTION)
# =====================================================================================

LEVEL4_DISJOINT_RETENTION_SQL = """
-- 5.2 Disjoint Bucket Retention (OPTIMIZED: Single query for all buckets)
WITH baseline AS (
    SELECT
        smmclid,
        COUNT(DISTINCT pdblid) AS baseline_pdblid_count,
        COUNT(DISTINCT (pdblid, clid)) AS baseline_pair_count
    FROM pdb_ligand_chembl_ligand
    GROUP BY smmclid
),
retention_data AS (
    SELECT
        smmclid,
        pdb_stid,
        chembl_type,
        CASE
            WHEN percent_identity = 100 THEN '100'
            WHEN percent_identity >= 95 THEN '95-<100'
            WHEN percent_identity >= 90 THEN '90-<95'
            WHEN percent_identity >= 85 THEN '85-<90'
            ELSE '80-<85'
        END AS identity_bucket_disjoint,
        pdblid,
        clid,
        'STRICT' AS variant
    FROM mv_ligand_sequence_target_fact
    UNION ALL
    SELECT
        smmclid,
        pdb_stid,
        chembl_type,
        CASE
            WHEN percent_identity = 100 THEN '100'
            WHEN percent_identity >= 95 THEN '95-<100'
            WHEN percent_identity >= 90 THEN '90-<95'
            WHEN percent_identity >= 85 THEN '85-<90'
            ELSE '80-<85'
        END AS identity_bucket_disjoint,
        pdblid,
        clid,
        'RELAXED' AS variant
    FROM mv_ligand_sequence_target_fact_nohitcov
),
aggregated AS (
    SELECT
        smmclid,
        pdb_stid,
        chembl_type,
        identity_bucket_disjoint,
        variant,
        COUNT(DISTINCT pdblid) AS surviving_pdblid_count,
        COUNT(DISTINCT (pdblid, clid)) AS surviving_pair_count
    FROM retention_data
    GROUP BY smmclid, pdb_stid, chembl_type, identity_bucket_disjoint, variant
)
SELECT
    a.smmclid,
    a.pdb_stid,
    a.chembl_type,
    a.identity_bucket_disjoint,
    b.baseline_pdblid_count,
    b.baseline_pair_count,
    MAX(CASE WHEN variant = 'STRICT' THEN surviving_pdblid_count END) AS surviving_pdblid_count_strict,
    MAX(CASE WHEN variant = 'RELAXED' THEN surviving_pdblid_count END) AS surviving_pdblid_count_relaxed,
    MAX(CASE WHEN variant = 'STRICT' THEN surviving_pair_count END) AS surviving_pair_count_strict,
    MAX(CASE WHEN variant = 'RELAXED' THEN surviving_pair_count END) AS surviving_pair_count_relaxed,
    ROUND(MAX(CASE WHEN variant = 'STRICT' THEN surviving_pdblid_count END) * 100.0 / b.baseline_pdblid_count, 2) AS retention_pdblid_pct_strict,
    ROUND(MAX(CASE WHEN variant = 'RELAXED' THEN surviving_pdblid_count END) * 100.0 / b.baseline_pdblid_count, 2) AS retention_pdblid_pct_relaxed,
    ROUND(MAX(CASE WHEN variant = 'STRICT' THEN surviving_pair_count END) * 100.0 / b.baseline_pair_count, 2) AS retention_pair_pct_strict,
    ROUND(MAX(CASE WHEN variant = 'RELAXED' THEN surviving_pair_count END) * 100.0 / b.baseline_pair_count, 2) AS retention_pair_pct_relaxed,
    ROUND((MAX(CASE WHEN variant = 'RELAXED' THEN surviving_pdblid_count END) -
           MAX(CASE WHEN variant = 'STRICT' THEN surviving_pdblid_count END)) * 100.0 / b.baseline_pdblid_count, 2) AS delta_pp_pdblid,
    ROUND((MAX(CASE WHEN variant = 'RELAXED' THEN surviving_pair_count END) -
           MAX(CASE WHEN variant = 'STRICT' THEN surviving_pair_count END)) * 100.0 / b.baseline_pair_count, 2) AS delta_pp_pair
FROM aggregated a
JOIN baseline b ON a.smmclid = b.smmclid
GROUP BY a.smmclid, a.pdb_stid, a.chembl_type, a.identity_bucket_disjoint, b.baseline_pdblid_count, b.baseline_pair_count
ORDER BY a.smmclid, a.pdb_stid, a.identity_bucket_disjoint
"""

LEVEL4_CUMULATIVE_RETENTION_SQL = """
WITH thresholds AS (
    SELECT 100 AS thr, '>=100' AS label, 1 AS ord UNION ALL
    SELECT 95, '>=95', 2 UNION ALL
    SELECT 90, '>=90', 3 UNION ALL
    SELECT 85, '>=85', 4 UNION ALL
    SELECT 80, '>=80', 5
),
baseline AS (
    SELECT
        smmclid,
        COUNT(DISTINCT pdblid) AS baseline_pdblid_count,
        COUNT(DISTINCT (pdblid, clid)) AS baseline_pair_count
    FROM pdb_ligand_chembl_ligand
    GROUP BY smmclid
),
combined AS (
    SELECT 'STRICT' AS variant, smmclid, pdblid, clid, percent_identity
    FROM mv_ligand_sequence_target_fact
    UNION ALL
    SELECT 'RELAXED' AS variant, smmclid, pdblid, clid, percent_identity
    FROM mv_ligand_sequence_target_fact_nohitcov
),
expanded AS (
    SELECT
        c.variant,
        c.smmclid,
        t.label AS identity_bucket_cumulative,
        t.ord,
        c.pdblid,
        c.clid
    FROM combined c
    JOIN thresholds t
      ON c.percent_identity >= t.thr
),
aggregated AS (
    SELECT
        variant,
        smmclid,
        identity_bucket_cumulative,
        ord,
        COUNT(DISTINCT pdblid) AS surviving_pdblid_count,
        COUNT(DISTINCT (pdblid, clid)) AS surviving_pair_count
    FROM expanded
    GROUP BY variant, smmclid, identity_bucket_cumulative, ord
),
pivoted AS (
    SELECT
        smmclid,
        identity_bucket_cumulative,
        ord,
        MAX(CASE WHEN variant='STRICT' THEN surviving_pdblid_count END) AS surviving_pdblid_count_strict,
        MAX(CASE WHEN variant='RELAXED' THEN surviving_pdblid_count END) AS surviving_pdblid_count_relaxed,
        MAX(CASE WHEN variant='STRICT' THEN surviving_pair_count END) AS surviving_pair_count_strict,
        MAX(CASE WHEN variant='RELAXED' THEN surviving_pair_count END) AS surviving_pair_count_relaxed
    FROM aggregated
    GROUP BY smmclid, identity_bucket_cumulative, ord
)
SELECT
    p.smmclid,
    p.identity_bucket_cumulative,
    b.baseline_pdblid_count,
    b.baseline_pair_count,
    p.surviving_pdblid_count_strict,
    p.surviving_pdblid_count_relaxed,
    p.surviving_pair_count_strict,
    p.surviving_pair_count_relaxed,
    ROUND(p.surviving_pdblid_count_strict * 100.0 / NULLIF(b.baseline_pdblid_count,0), 2)::float AS strict_retention_pdblid_pct,
    ROUND(p.surviving_pdblid_count_relaxed * 100.0 / NULLIF(b.baseline_pdblid_count,0), 2)::float AS relaxed_retention_pdblid_pct,
    ROUND(p.surviving_pair_count_strict * 100.0 / NULLIF(b.baseline_pair_count,0), 2)::float AS strict_retention_pair_pct,
    ROUND(p.surviving_pair_count_relaxed * 100.0 / NULLIF(b.baseline_pair_count,0), 2)::float AS relaxed_retention_pair_pct,
    ROUND((p.surviving_pdblid_count_relaxed - p.surviving_pdblid_count_strict) * 100.0 / NULLIF(b.baseline_pdblid_count,0), 2)::float AS delta_pp_pdblid,
    ROUND((p.surviving_pair_count_relaxed - p.surviving_pair_count_strict) * 100.0 / NULLIF(b.baseline_pair_count,0), 2)::float AS delta_pp_pair
FROM pivoted p
JOIN baseline b USING (smmclid)
ORDER BY p.smmclid, p.ord;
"""

LEVEL4_RELAXED_ONLY_SQL = """
-- 7.1 Relaxed-Only Composition (Rewritten for PostgreSQL 10 efficiency: LEFT JOIN anti-join + pre-bucketing)
WITH strict_keys AS (
    SELECT DISTINCT pdblid, query_sid, clid, chembl_tid
    FROM mv_ligand_sequence_target_fact
),
relaxed_only_quadruplets AS (
    SELECT r.*
    FROM mv_ligand_sequence_target_fact_nohitcov r
    LEFT JOIN strict_keys s
      ON r.pdblid = s.pdblid
     AND r.query_sid = s.query_sid
     AND r.clid = s.clid
     AND r.chembl_tid = s.chembl_tid
    WHERE s.pdblid IS NULL
),
relaxed_only_buckets AS (
    SELECT
        smmclid,
        CASE
            WHEN percent_identity >= 95 THEN '>=95'
            WHEN percent_identity >= 90 THEN '>=90'
            WHEN percent_identity >= 85 THEN '>=85'
            ELSE '>=80'
        END AS identity_bucket_cumulative,
        CASE
            WHEN hit_coverage IS NULL THEN 'NULL'
            WHEN hit_coverage >= 70 THEN '>=70'
            WHEN hit_coverage >= 40 THEN '40-<70'
            WHEN hit_coverage >= 10 THEN '10-<40'
            ELSE '<10'
        END AS hit_coverage_quality
    FROM relaxed_only_quadruplets
),
relaxed_only_agg AS (
    SELECT
        smmclid,
        identity_bucket_cumulative,
        hit_coverage_quality,
        COUNT(*) AS relaxed_only_quad_ct
    FROM relaxed_only_buckets
    GROUP BY smmclid, identity_bucket_cumulative, hit_coverage_quality
)
SELECT
    a.smmclid,
    a.identity_bucket_cumulative,
    a.hit_coverage_quality,
    a.relaxed_only_quad_ct,
    ROUND(a.relaxed_only_quad_ct * 100.0 /
          NULLIF(SUM(a.relaxed_only_quad_ct) OVER (PARTITION BY a.smmclid),0), 2)::float AS pct_of_relaxed_only_smmclid,
    ROUND(a.relaxed_only_quad_ct * 100.0 /
          NULLIF(SUM(a.relaxed_only_quad_ct) OVER (),0), 2)::float AS pct_of_total_relaxed_only
FROM relaxed_only_agg a
ORDER BY a.smmclid, a.identity_bucket_cumulative, a.hit_coverage_quality
"""

LEVEL4_HIT_COVERAGE_QUALITY_SQL = """
-- 6.1 Hit Coverage Quality Distribution (OPTIMIZED: With cumulative percentages)
WITH coverage_analysis AS (
    SELECT
        CASE
            WHEN hit_coverage IS NULL THEN 'NULL'
            WHEN hit_coverage >= 80 THEN '>=80 (Strict)'
            WHEN hit_coverage >= 70 THEN '70-<80'
            WHEN hit_coverage >= 40 THEN '40-<70'
            WHEN hit_coverage >= 10 THEN '10-<40'
            ELSE '<10'
        END AS hit_coverage_bucket,
        COUNT(*) as row_count
    FROM mv_ligand_sequence_target_fact_nohitcov
    GROUP BY CASE
        WHEN hit_coverage IS NULL THEN 'NULL'
        WHEN hit_coverage >= 80 THEN '>=80 (Strict)'
        WHEN hit_coverage >= 70 THEN '70-<80'
        WHEN hit_coverage >= 40 THEN '40-<70'
        WHEN hit_coverage >= 10 THEN '10-<40'
        ELSE '<10'
    END
)
SELECT
    hit_coverage_bucket,
    row_count,
    ROUND(row_count * 100.0 / SUM(row_count) OVER(), 2) as percentage,
    SUM(row_count) OVER (ORDER BY
        CASE hit_coverage_bucket
            WHEN 'NULL' THEN 0
            WHEN '>=80 (Strict)' THEN 1
            WHEN '70-<80' THEN 2
            WHEN '40-<70' THEN 3
            WHEN '10-<40' THEN 4
            ELSE 5
        END
    ) AS cumulative_count,
    ROUND(SUM(row_count) OVER (ORDER BY
        CASE hit_coverage_bucket
            WHEN 'NULL' THEN 0
            WHEN '>=80 (Strict)' THEN 1
            WHEN '70-<80' THEN 2
            WHEN '40-<70' THEN 3
            WHEN '10-<40' THEN 4
            ELSE 5
        END
    ) * 100.0 / SUM(row_count) OVER(), 2) AS cumulative_percentage
FROM coverage_analysis
ORDER BY
    CASE hit_coverage_bucket
        WHEN 'NULL' THEN 0
        WHEN '>=80 (Strict)' THEN 1
        WHEN '70-<80' THEN 2
        WHEN '40-<70' THEN 3
        WHEN '10-<40' THEN 4
        ELSE 5
    END
"""

# =====================================================================================
# SECTION 8: INTEGRITY CHECKS & VALIDATION QUERIES
# =====================================================================================

LEVEL4_INTEGRITY_STRICT_SUBSET_SQL = """
-- Check 1: Verify Strict ⊆ Relaxed
SELECT
    'Strict ⊆ Relaxed' AS check_name,
    COUNT(*) AS violation_count,
    CASE WHEN COUNT(*) = 0 THEN 'PASS' ELSE 'FAIL' END AS result
FROM (
    SELECT pdblid, query_sid, clid, chembl_tid
    FROM mv_ligand_sequence_target_fact
    EXCEPT
    SELECT pdblid, query_sid, clid, chembl_tid
    FROM mv_ligand_sequence_target_fact_nohitcov
) violations
"""

LEVEL4_INTEGRITY_CUM_MONO_SQL = """
WITH ordered AS (
    SELECT
        smmclid,
        identity_bucket_cumulative,
        strict_retention_pair_pct::float AS strict_retention_pair_pct,
        relaxed_retention_pair_pct::float AS relaxed_retention_pair_pct,
        ROW_NUMBER() OVER (PARTITION BY smmclid ORDER BY
            CASE identity_bucket_cumulative
                WHEN '>=100' THEN 1
                WHEN '>=95' THEN 2
                WHEN '>=90' THEN 3
                WHEN '>=85' THEN 4
                WHEN '>=80' THEN 5
            END
        ) AS rn
    FROM (
        SELECT
            p.smmclid,
            p.identity_bucket_cumulative,
            ROUND(p.surviving_pair_count_strict * 100.0 / NULLIF(b.baseline_pair_count,0), 2)::float AS strict_retention_pair_pct,
            ROUND(p.surviving_pair_count_relaxed * 100.0 / NULLIF(b.baseline_pair_count,0), 2)::float AS relaxed_retention_pair_pct
        FROM (
            WITH thresholds AS (
                SELECT 100 AS thr, '>=100' AS label, 1 AS ord UNION ALL
                SELECT 95, '>=95', 2 UNION ALL
                SELECT 90, '>=90', 3 UNION ALL
                SELECT 85, '>=85', 4 UNION ALL
                SELECT 80, '>=80', 5
            ),
            combined AS (
                SELECT 'STRICT' AS variant, smmclid, pdblid, clid, percent_identity
                FROM mv_ligand_sequence_target_fact
                UNION ALL
                SELECT 'RELAXED', smmclid, pdblid, clid, percent_identity
                FROM mv_ligand_sequence_target_fact_nohitcov
            ),
            expanded AS (
                SELECT
                    c.variant,
                    c.smmclid,
                    t.label AS identity_bucket_cumulative,
                    t.ord,
                    c.pdblid,
                    c.clid
                FROM combined c
                JOIN thresholds t ON c.percent_identity >= t.thr
            ),
            aggregated AS (
                SELECT
                    variant,
                    smmclid,
                    identity_bucket_cumulative,
                    ord,
                    COUNT(DISTINCT pdblid) AS surviving_pdblid_count,
                    COUNT(DISTINCT (pdblid, clid)) AS surviving_pair_count
                FROM expanded
                GROUP BY variant, smmclid, identity_bucket_cumulative, ord
            )
            SELECT
                a.smmclid,
                a.identity_bucket_cumulative,
                MAX(CASE WHEN a.variant='STRICT' THEN a.surviving_pair_count END) AS surviving_pair_count_strict,
                MAX(CASE WHEN a.variant='RELAXED' THEN a.surviving_pair_count END) AS surviving_pair_count_relaxed,
                a.ord
            FROM aggregated a
            GROUP BY a.smmclid, a.identity_bucket_cumulative, a.ord
        ) p
        JOIN (
            SELECT
                smmclid,
                COUNT(DISTINCT (pdblid, clid)) AS baseline_pair_count
            FROM pdb_ligand_chembl_ligand
            GROUP BY smmclid
        ) b ON p.smmclid = b.smmclid
    ) q
),
lagged AS (
    SELECT
        *,
        LAG(strict_retention_pair_pct) OVER (PARTITION BY smmclid ORDER BY rn) AS prev_strict,
        LAG(relaxed_retention_pair_pct) OVER (PARTITION BY smmclid ORDER BY rn) AS prev_relaxed
    FROM ordered
)
SELECT
    'Cumulative Monotonicity (non-decreasing)' AS check_name,
    COUNT(*) FILTER (
        WHERE (strict_retention_pair_pct < prev_strict)
           OR (relaxed_retention_pair_pct < prev_relaxed)
    ) AS violation_count,
    CASE
        WHEN COUNT(*) FILTER (
            WHERE (strict_retention_pair_pct < prev_strict)
               OR (relaxed_retention_pair_pct < prev_relaxed)
        ) = 0 THEN 'PASS' ELSE 'FAIL'
    END AS result
FROM lagged
WHERE prev_strict IS NOT NULL OR prev_relaxed IS NOT NULL;
"""

LEVEL4_INTEGRITY_LEVEL1_VS_LEVEL2_SQL = """
-- Check 3: Level 1 vs Level 2 Consistency
WITH level1_pairs AS (
    SELECT smmclid, COUNT(DISTINCT (pdblid, clid)) AS pair_count
    FROM pdb_ligand_chembl_ligand
    GROUP BY smmclid
),
level2_pairs AS (
    SELECT smmclid, COUNT(DISTINCT (pdblid, clid)) AS pair_count
    FROM mv_ligand_sequence_target_fact_nohitcov
    GROUP BY smmclid
)
SELECT
    'Level2 ≤ Level1 Pairs' AS check_name,
    COUNT(*) AS violation_count,
    CASE WHEN COUNT(*) = 0 THEN 'PASS' ELSE 'FAIL' END AS result
FROM (
    SELECT l1.smmclid
    FROM level1_pairs l1
    LEFT JOIN level2_pairs l2 ON l1.smmclid = l2.smmclid
    WHERE COALESCE(l2.pair_count, 0) > l1.pair_count
) violations
"""
# =====================================================================================
# SECTION 6: ACTIVITY-QUALIFIED PAIR IDENTITY DISTRIBUTION (COVERAGE-AGNOSTIC)
# =====================================================================================
#
# Counts distinct (pdblid, clid) ligand mapping pairs that are supported by at least one
# ChEMBL activity (ANY activity row, no pchembl_value restriction) for a target that is
# reachable via a BLAST best-hit sequence mapping (percent_identity >= 80) without
# applying any query_coverage or hit_coverage filters.
#
# Critical corrected join path (compared to earlier activity views):
#   sequence_target.sid = best.hit_sid
#   target.chembl_tid   = sequence_target.chembl_tid
#   assay.tid           = target.tid
#   activity.asid       = assay.asid
#   activity.molregno   = chembl_ligand.molregno
#
# Structural mapping connection:
#   pdb_sequence (query_sid) -> pdb_ligand_pdb_id (pid,cid) -> pdblid
#   pdblid + clid structural mapping table: pdb_ligand_chembl_ligand (for smmclid)
#
# Interval notation (mathematical):
#   [100], [95,100), [90,95), [85,90), [80,85)
#
# Output columns:
#   identity_interval  : interval bucket
#   smmclid            : best structural confidence (max smmclid per pair)
#   pair_count         : number of distinct pairs in that (interval, smmclid)
#
# NOTE: This query intentionally ignores coverage thresholds and includes ANY
#       activity row (no pchembl_value filter).
#
PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_SQL = """
-- Staged, ligand-first, activity-qualified best pair identity distribution.
-- Source MV: mv_ligand_supported_pair_best_identity (refresh the staged chain first).
-- Buckets are mutually exclusive: [100], [95,100), [90,95), [85,90), [80,85)
WITH bucketed AS (
    SELECT
        identity_interval,
        best_smmclid AS smmclid
    FROM mv_ligand_supported_pair_best_identity
)
SELECT
    identity_interval,
    smmclid,
    COUNT(*) AS pair_count
FROM bucketed
GROUP BY identity_interval, smmclid
ORDER BY
    CASE identity_interval
        WHEN '[100]' THEN 1
        WHEN '[95,100)' THEN 2
        WHEN '[90,95)' THEN 3
        WHEN '[85,90)' THEN 4
        ELSE 5
    END,
    smmclid;
"""

# --------------------------------------------------------------------------------------
# Cumulative (Multi-Bucket) Distribution
# --------------------------------------------------------------------------------------
# A pair can contribute to multiple identity intervals if it has at least one alignment
# (percent_identity >=80) in each interval. Uses the multi-bucket expansion MV
# mv_ligand_supported_pair_identity_multibucket. Counts are distinct (pdblid, clid)
# per (identity_interval, smmclid).
PAIR_IDENTITY_DISTRIBUTION_ACTIVITY_CUMULATIVE_SQL = """
-- Multi-bucket (cumulative) activity-qualified pair identity distribution.
-- Source MV: mv_ligand_supported_pair_identity_multibucket
WITH base AS (
    SELECT
        identity_interval,
        smmclid,
        pdblid,
        clid
    FROM mv_ligand_supported_pair_identity_multibucket
),
agg AS (
    SELECT
        identity_interval,
        smmclid,
        COUNT(DISTINCT (pdblid, clid)) AS pair_count
    FROM base
    GROUP BY identity_interval, smmclid
)
SELECT
    identity_interval,
    smmclid,
    pair_count
FROM agg
ORDER BY
    CASE identity_interval
        WHEN '[100]' THEN 1
        WHEN '[95,100)' THEN 2
        WHEN '[90,95)' THEN 3
        WHEN '[85,90)' THEN 4
        ELSE 5
    END,
    smmclid;
"""
