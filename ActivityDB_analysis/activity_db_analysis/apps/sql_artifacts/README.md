# SQL Artifacts Management

Comprehensive management system for all SQL artifacts (functions, materialized views, views, and indexes) used by the ActivityDB analysis pipeline.

This app provides:
- Fine-grained control over creating, dropping, and refreshing artifacts
- Automatic dependency resolution and topological ordering
- Status monitoring and validation
- Batch operations for full pipeline management

## Overview

The SQL artifacts system manages over 100 database objects required for ligand mapping, sequence mapping, and activity qualification analyses. All artifact definitions are embedded in Python for full control, type safety, and version tracking.

## Why Managed SQL Artifacts?

This system provides significant advantages over managing raw SQL files:

| Benefit | Description |
|---------|-------------|
| **Dependency tracking** | Automatically creates dependencies in correct order—no manual ordering of scripts |
| **Version control** | All SQL embedded in Python, tracked with code changes in the same commit |
| **Type safety** | Python dataclasses (`SQLArtifact`) ensure artifact definitions are valid at import time |
| **Idempotency** | Safe to run create/refresh multiple times—existing objects are handled gracefully |
| **Status monitoring** | Easy to check what exists vs. what's missing with `sql-manager status` |
| **Reproducibility** | Anyone can recreate the full database schema from code with `create-all` |
| **Documentation** | Each artifact has description, dependencies, and type metadata inline |

Without this system, you would need to:
- Manually track which SQL scripts depend on others
- Run scripts in the correct order (or debug cryptic "relation does not exist" errors)
- Maintain separate documentation about what each script does
- Hope that SQL files stay in sync with the codebase

## Adding New Artifacts

To add a new materialized view, function, view, or index:

1. **Open** `sql_artifacts_management.py`
2. **Add** a new `SQLArtifact()` entry to the `ARTIFACTS` dict
3. **Specify**: `name`, `artifact_type`, `sql`, `dependencies`, `description`
4. The system will **automatically handle dependency ordering** when creating/dropping

### Example: Adding a New Materialized View

```python
"mv_my_new_analysis": SQLArtifact(
    name="mv_my_new_analysis",
    sql="""CREATE MATERIALIZED VIEW mv_my_new_analysis AS
SELECT
    t.chembl_tid,
    COUNT(DISTINCT l.clid) AS ligand_count
FROM mv_ligand_activity_target t
JOIN chembl_ligand l ON l.clid = t.clid
GROUP BY t.chembl_tid;""",
    drop_sql="DROP MATERIALIZED VIEW IF EXISTS mv_my_new_analysis CASCADE",
    refresh_sql="REFRESH MATERIALIZED VIEW mv_my_new_analysis",
    refresh_concurrent_sql="REFRESH MATERIALIZED VIEW CONCURRENTLY mv_my_new_analysis",
    dependencies=["mv_ligand_activity_target"],  # List MVs/functions this depends on
    description="Per-target ligand counts for my analysis",
    artifact_type=ArtifactType.MATERIALIZED_VIEW,
),
```

### Key Fields

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Artifact name (must match dict key) |
| `sql` | Yes | CREATE statement |
| `drop_sql` | Yes | DROP statement (use CASCADE if appropriate) |
| `dependencies` | Yes | List of artifact names this depends on (can be empty `[]`) |
| `description` | Yes | Human-readable description |
| `artifact_type` | Yes | One of `ArtifactType.FUNCTION`, `MATERIALIZED_VIEW`, `VIEW`, `INDEX` |
| `refresh_sql` | No | For MVs: standard refresh statement |
| `refresh_concurrent_sql` | No | For MVs: concurrent refresh (requires unique index) |

After adding, verify with:
```bash
uv run actfind sql-manager list | grep mv_my_new_analysis
uv run actfind sql-manager dependencies mv_my_new_analysis
```

### Artifact Types

| Type | Count | Description |
|------|-------|-------------|
| `function` | ~10 | SQL/PL-pgSQL functions for data transformation |
| `materialized_view` | ~30 | Pre-computed result sets for performance |
| `view` | ~10 | Logical views for convenience queries |
| `index` | ~60 | Performance indexes on materialized views |

## CLI Usage

### Single Artifact Operations

```bash
# Create a specific artifact (with dependencies)
uv run actfind sql-manager create mv_alignment_sequence_metrics

# Drop an artifact
uv run actfind sql-manager drop mv_alignment_sequence_metrics

# Drop with cascade (removes dependents)
uv run actfind sql-manager drop mv_alignment_sequence_metrics --cascade

# Refresh a materialized view
uv run actfind sql-manager refresh mv_alignment_sequence_metrics

# Concurrent refresh (requires unique index)
uv run actfind sql-manager refresh mv_alignment_sequence_metrics --concurrent
```

### Batch Operations

```bash
# Create all artifacts in dependency order
uv run actfind sql-manager create-all

# Create only materialized views
uv run actfind sql-manager create-all --type materialized_view

# Dry run (show what would be created)
uv run actfind sql-manager create-all --dry-run

# Refresh all materialized views
uv run actfind sql-manager refresh-all-mvs

# Refresh all concurrently where supported
uv run actfind sql-manager refresh-all-mvs --concurrent

# Drop all artifacts (use with caution!)
uv run actfind sql-manager drop-all --cascade
```

### Status and Validation

```bash
# Show status of all artifacts
uv run actfind sql-manager status

# Show status of specific artifact
uv run actfind sql-manager status mv_alignment_sequence_metrics

# List all available artifacts
uv run actfind sql-manager list

# List only functions
uv run actfind sql-manager list --type function

# Show dependencies for an artifact
uv run actfind sql-manager dependencies mv_alignment_sequence_metrics

# Validate all artifacts
uv run actfind sql-manager validate
```

## Key Materialized Views

### Sequence Mapping

| View | Purpose |
|------|---------|
| `mv_alignment_sequence_metrics` | Sequence-level alignment metrics for target mapping |
| `mv_alignment_sequence_metrics_chain` | Chain-level alignment metrics (intersection cohort) |
| `mv_sequence_totals_chain` | Chain counts for diagnostics |
| `mv_best_alignment_identity_floor` | Canonical best alignment per (query_sid, hit_sid) |

### Ligand Mapping

| View | Purpose |
|------|---------|
| `mv_ligand_sequence_target_fact` | Strict coverage tier (identity + query + hit ≥80%) |
| `mv_ligand_sequence_target_fact_nohitcov` | Relaxed tier (no hit coverage filter) |
| `mv_ligand_sequence_target_fact_identity_only` | Identity-only tier (no coverage filters) |
| `mv_ligand_supported_pair_best_identity` | Activity-qualified best pair identity |
| `mv_ligand_supported_pair_identity_multibucket` | Multi-bucket cumulative distribution |

### Activity Pipeline

| View | Purpose |
|------|---------|
| `mv_pdb_ligand_chembl_activity_assay_base` | Full ligand→activity→assay convergence |
| `mv_pdb_ligand_assay_link` | Collapsed ligand↔assay link |
| `mv_active_pair_ligands` | Active pair ligands for activity qualification |
| `mv_ligand_activity_target` | Ligand to activity to target mapping |

### Presence Analysis

| View | Purpose |
|------|---------|
| `mv_component_sequence_presence_identity_only` | Component route presence |
| `mv_variant_sequence_presence_identity_only` | Variant route presence |
| `mv_component_variant_presence_overlap` | Overlap between routes |
| `mv_presence_quadruplet_matrix` | Distinct quadruplets with membership flags |

## Key Functions

| Function | Purpose |
|----------|---------|
| `f_ligand_sequence_target_base` | Base function with parameterized coverage filters |
| `f_ligand_sequence_target_base_v2` | Unified-tier version using best alignment MV |
| `f_component_sequence_presence` | Component route presence function |
| `f_variant_sequence_presence` | Variant route presence function |
| `refresh_ligand_mapping_mvs` | Orchestrated refresh for mapping MVs |
| `refresh_sequence_presence_mvs` | Orchestrated refresh for presence pipeline |
| `get_ligand_mapping_mv_status` | Status reporting function |
| `analyze_ligand_mapping_performance` | Performance analysis function |

## Tier Views

The system provides three coverage tiers via views on `mv_best_alignment_identity_floor`:

| View | Filter |
|------|--------|
| `v_best_alignment_id` | Identity ≥80% only |
| `v_best_alignment_id_query` | + query coverage ≥80% |
| `v_best_alignment_id_query_hit` | + hit coverage ≥80% |

## Dependency Management

Artifacts are automatically created in dependency order. For example, creating `mv_alignment_sequence_metrics` will first create any missing dependencies.

### View Dependency Graph

```bash
# Show what would be created for a complex artifact
uv run actfind sql-manager dependencies mv_ligand_supported_pair_best_identity
```

Example output:
```
Artifact: mv_ligand_supported_pair_best_identity
Type: materialized_view
Description: Activity-qualified best pair identity distribution

Direct dependencies (3):
  - mv_active_pair_ligands (materialized_view)
  - mv_ligand_activity_target (materialized_view)
  - mv_best_alignment_identity80 (materialized_view)

Dependents (2):
  - mv_ligand_supported_pair_best_identity_interval_idx (index)
  - mv_ligand_supported_pair_best_identity_smmclid_idx (index)

Full creation order (5):
   1. mv_active_pair_ligands
   2. mv_ligand_activity_target
   3. mv_best_alignment_identity80
   4. mv_ligand_supported_pair_best_identity
👉 5. mv_ligand_supported_pair_best_identity
```

## Connection Options

All commands accept database connection parameters:

| Option | Description |
|--------|-------------|
| `--user` | ActivityDB user override |
| `--pwd` | ActivityDB password override |
| `--host` | ActivityDB host override |
| `--port` | ActivityDB port override |
| `--database-name` | Database name override |
SSH tunnels are managed automatically by the CLI. Use `--tunnel-mode` at the root level:
```bash
# Status with automatic tunnel management (default)
uv run actfind sql-manager status

# Explicit tunnel mode control  
uv run actfind --tunnel-mode required sql-manager status
```

## Programmatic Usage

```python
from activity_db_analysis.apps.sql_artifacts.sql_artifacts_management import (
    ARTIFACTS,
    DependencyManager,
    DatabaseExecutor,
    ArtifactType,
)
from activity_db_analysis.db.connector import build_unified_connector

# Get artifact info
artifact = ARTIFACTS["mv_alignment_sequence_metrics"]
print(f"Type: {artifact.artifact_type}")
print(f"Dependencies: {artifact.dependencies}")
print(f"Description: {artifact.description}")

# Check dependency order
dep_mgr = DependencyManager(ARTIFACTS)
order = dep_mgr.get_creation_order(["mv_alignment_sequence_metrics"])
print(f"Creation order: {order}")

# List all artifacts by type
mvs = [name for name, a in ARTIFACTS.items() 
       if a.artifact_type == ArtifactType.MATERIALIZED_VIEW]
```

## Common Workflows

### Fresh Database Setup

```bash
# Create all artifacts from scratch
uv run actfind sql-manager create-all

# Verify everything was created
uv run actfind sql-manager validate
```

### Refresh After Data Update

```bash
# Refresh all materialized views in dependency order
uv run actfind sql-manager refresh-all-mvs

# Or refresh specific views
uv run actfind sql-manager refresh mv_best_alignment_identity_floor
uv run actfind sql-manager refresh mv_alignment_sequence_metrics
```

### Sequence Mapping Setup

```bash
# Create only sequence mapping artifacts
uv run actfind sql-manager create v_best_alignment_id
uv run actfind sql-manager create v_best_alignment_id_query
uv run actfind sql-manager create v_best_alignment_id_query_hit
uv run actfind sql-manager create mv_alignment_sequence_metrics
uv run actfind sql-manager create mv_alignment_sequence_metrics_chain
```

### Ligand Mapping Pipeline

```bash
# Create the staged ligand pipeline
uv run actfind sql-manager create mv_pdb_ligand_chembl_activity_assay_base
uv run actfind sql-manager create mv_pdb_ligand_assay_link
uv run actfind sql-manager create mv_ligand_supported_pair_best_identity

# Refresh the pipeline
uv run actfind sql-manager refresh mv_pdb_ligand_chembl_activity_assay_base
uv run actfind sql-manager refresh mv_pdb_ligand_assay_link
uv run actfind sql-manager refresh mv_ligand_supported_pair_best_identity
```

## Justfile Recipes

Common workflows are available via the project justfile:

```bash
# Full pipeline for sequence mapping (creates MVs, extracts data, generates plots)
just sequence-mapping-single-mode identity-only ./artifacts

# Build all publication artifacts
just build-all-artifacts
```

## Troubleshooting

### Artifact Creation Fails

```bash
# Check if dependencies exist
uv run actfind sql-manager dependencies <artifact_name>

# Force recreate (drops existing first)
uv run actfind sql-manager create <artifact_name> --force
```

### Missing Dependencies

```bash
# Validate all definitions
uv run actfind sql-manager validate

# Create all dependencies first
uv run actfind sql-manager create-all --type function
uv run actfind sql-manager create-all --type materialized_view
```

### Concurrent Refresh Fails

Concurrent refresh requires a unique index. Check if the artifact supports it:

```bash
uv run actfind sql-manager list --type materialized_view
# Look for artifacts with concurrent refresh support
```

### Performance Issues

```bash
# Check MV sizes and row counts
uv run actfind sql-manager status

# Ensure indexes exist
uv run actfind sql-manager create-all --type index
```

## Status Output Format

```
Name                                               Type                 Exists     Size            Rows           
------------------------------------------------------------------------------------------------------------------------
f_ligand_sequence_target_base                      function             ✓          N/A             N/A            
mv_alignment_sequence_metrics                      materialized_view    ✓          1.2 GB          12,345,678     
mv_alignment_sequence_metrics_chain                materialized_view    ✗          N/A             N/A            
```

## License

This app is part of the ActivityDB Analysis repository and follows the repository's license.
