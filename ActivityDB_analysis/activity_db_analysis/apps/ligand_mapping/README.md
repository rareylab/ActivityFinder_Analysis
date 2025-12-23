# Ligand Mapping Analysis

Analysis and visualization of PDB-ChEMBL ligand pair identity distributions with activity qualification. This app generates two publication artifacts:
- **ligand_match_levels.tex** — PDB-ChEMBL ligand mapping confidence level table
- **pair_identity_distribution_activity_combined_panel** — 2×2 panel plot

## Quick Start

```bash
# Using justfile (recommended)
just ligand-match-level-counts                      # LaTeX table
just pair-identity-distribution-activity-data-pull  # Data extraction
just pair-identity-distribution-activity-plots      # Generate plots

# Build all artifacts at once
just build-all-artifacts
```

## Why Use This

This app answers: **How well can we link PDB ligands to ChEMBL compounds with known bioactivity data?**

- Quantifies ligand mapping success rates across sequence identity thresholds (80%–100%)
- Shows what proportion of structural ligands can be connected to activity databases
- Essential for understanding SAR coverage before running activity predictions
- Helps researchers assess data completeness in structure-activity analyses

## Features

- Pair identity distribution analysis across sequence identity intervals
- Activity-qualified ligand mapping (requires ChEMBL activity row support)
- Stacked bar visualizations with structural confidence level breakdown
- LaTeX table generation for publication-ready outputs
- Combined panel plots comparing disjoint vs cumulative distributions

## CLI Usage

### Generate Data (CSV + LaTeX)

```bash
# Generate disjoint pair identity distribution
uv run actfind pair-identity-distribution-activity --output-dir ./tables

# Generate cumulative (multi-bucket) distribution
uv run actfind pair-identity-distribution-activity --output-dir ./tables --cumulative
```

### Generate Ligand Match Level Counts

```bash
uv run actfind ligand-match-level-counts --output-dir ./tables
```

### Generate Plots

```bash
# Create counts and normalized stacked bar plots
uv run actfind plot-pair-identity-distribution-activity \
    --input ./tables/pair_identity_distribution_activity.csv \
    --output-dir ./figures

# Create combined 2×2 panel (disjoint + cumulative)
uv run actfind plot-pair-identity-distribution-activity \
    --input ./tables/pair_identity_distribution_activity.csv \
    --cumulative-input ./tables/pair_identity_distribution_activity_cumulative.csv \
    --output-dir ./figures \
    --combined --pdf
```

## Configuration

Database configuration via `appconfig.toml`:

```toml
[activity_db]
name = "activitydb_test"
user = "activityuser_test"
host = "127.0.0.1"
port = 5433
schema = "public"

[ssh]
host_alias = "your_ssh_alias"
remote_bind_host = "your.database.server.example.com"
remote_bind_port = 5432
local_bind_host = "127.0.0.1"
local_bind_port = 5433
```

Password via environment variable:
```bash
export ACTIVITY_DB_TEST_PWD="your_password"
```

## Prerequisites

The analysis requires materialized views managed by the `sql_artifacts` app:

```bash
uv run actfind sql-manager create mv_ligand_supported_pair_best_identity
uv run actfind sql-manager refresh mv_ligand_supported_pair_best_identity
```

## Output Files

| File | Description |
|------|-------------|
| `pair_identity_distribution_activity.csv` | Disjoint distribution data |
| `pair_identity_distribution_activity_cumulative.csv` | Cumulative mode data |
| `ligand_match_levels.tex` | Confidence level table |
| `*_counts.png/.pdf` | Absolute counts stacked bar |
| `*_percents.png/.pdf` | Normalized composition stacked bar |
| `*_combined_panel.png/.pdf` | 2×2 panel (disjoint + cumulative) |

## LaTeX Requirements

Generated `.tex` files require:
```latex
\usepackage{booktabs}    % Professional table formatting
\usepackage{siunitx}     % Number formatting with S columns
```

## Part of Publication Artifacts

This app generates 2 of the 6 publication artifacts. Run all with:
```bash
just build-all-artifacts
```
