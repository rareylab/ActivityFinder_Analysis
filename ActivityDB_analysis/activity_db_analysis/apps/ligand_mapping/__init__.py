"""
Ligand Mapping Analysis Package
================================

This package contains modules for analyzing PDB-ChEMBL ligand mappings
at different analytical levels (1-4), producing standardized CSV and
LaTeX table outputs.

Modules:
--------
- ligand_mapping_sql: SQL query constants for all analysis levels
- latex_utils: LaTeX table generation utilities with auto-formatting
- LigandMappingLevel1: Structure-only mapping analysis
- LigandMappingLevel2: Quadruplet universe analysis
- LigandMappingLevel3: Partitioned quadruplets analysis
- LigandMappingLevel4: Identity-qualified retention analysis
"""

from .ligand_mapping_sql import (
    LEVEL1_BASELINE_SQL,
    LEVEL1_SMMCLID_SQL,
    LEVEL1_FANOUT_SQL,
    LEVEL1_UNMAPPED_SQL,
    LEVEL2_QUADRUPLETS_SQL,
    LEVEL2_EXPANSION_SQL,
    LEVEL3_PARTITIONS_SQL,
    LEVEL3_EXPANSION_DISTINCT_SQL,
    LEVEL4_DISJOINT_RETENTION_SQL,
    LEVEL4_CUMULATIVE_RETENTION_SQL,
    LEVEL4_RELAXED_ONLY_SQL,
    LEVEL4_HIT_COVERAGE_QUALITY_SQL,
    LEVEL4_INTEGRITY_STRICT_SUBSET_SQL,
    LEVEL4_INTEGRITY_CUM_MONO_SQL,
    LEVEL4_INTEGRITY_LEVEL1_VS_LEVEL2_SQL,
)

from .latex_utils import (
    render_latex_table,
    build_latex_table,
)

__all__ = [
    # SQL constants
    'LEVEL1_BASELINE_SQL',
    'LEVEL1_SMMCLID_SQL',
    'LEVEL1_FANOUT_SQL',
    'LEVEL1_UNMAPPED_SQL',
    'LEVEL2_QUADRUPLETS_SQL',
    'LEVEL2_EXPANSION_SQL',
    'LEVEL3_PARTITIONS_SQL',
    'LEVEL3_EXPANSION_DISTINCT_SQL',
    'LEVEL4_DISJOINT_RETENTION_SQL',
    'LEVEL4_CUMULATIVE_RETENTION_SQL',
    'LEVEL4_RELAXED_ONLY_SQL',
    'LEVEL4_HIT_COVERAGE_QUALITY_SQL',
    'LEVEL4_INTEGRITY_STRICT_SUBSET_SQL',
    'LEVEL4_INTEGRITY_CUM_MONO_SQL',
    'LEVEL4_INTEGRITY_LEVEL1_VS_LEVEL2_SQL',
    # LaTeX utilities
    'render_latex_table',
    'build_latex_table',
]
