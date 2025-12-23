# HET Code Variants Analysis

This module analyzes and visualizes the distribution of variants per HET code in the ActivityFinder database. It generates the **variants_per_het_code_histogram** publication artifact.

## Quick Start

```bash
# Using justfile (recommended)
just plot-het

# Using CLI directly
uv run actfind plot-het-code-variants --output-dir ./output
```

## Overview

The analysis examines how many different variants (unique SMILES representations) exist for each HET code in the PDB. This helps understand the chemical diversity represented by each HET code identifier.

## Why HET Codes Have Multiple Variants

A single HET code (e.g., ATP) can have multiple SMILES representations due to:

- **Stereoisomers**: R/S configurations, cis/trans isomerism
- **Tautomers**: Proton shifts between atoms (keto-enol, heterocycle shifts)
- **Protonation states**: Different ionization states at varying pH
- **Representational differences**: Explicit vs implicit hydrogens, aromatic notation

Understanding variant counts helps assess chemical complexity, matching challenges, and data quality.

## CLI Usage

```bash
uv run actfind plot-het-code-variants \
    --output-dir ./output \
    --output-prefix variants_analysis
```

### Command Options

| Option | Description |
|--------|-------------|
| `--output-dir` | Output directory (default: `./output`) |
| `--output-prefix` | Base name for output files (default: `variants_per_het_code`) |
| `--latex/--no-latex` | Generate LaTeX table (default: true) |
| `--create-plots/--no-create-plots` | Create histogram plots (default: true) |

## Configuration

Database connection via `appconfig.toml`:

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
export ACT_DB_TEST_PWD="your_password"
```

## Output Files

| File | Description |
|------|-------------|
| `{prefix}_data.csv` | Raw data (num_variants, num_ligands, percentage) |
| `{prefix}_histogram.png/.pdf` | Distribution plot with log-scale y-axis |
| `{prefix}_table.tex` | Publication-ready LaTeX table |

## LaTeX Requirements

Generated `.tex` files require:
```latex
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{siunitx}
```

## Part of Publication Artifacts

This app generates `variants_per_het_code_histogram`, one of the 6 publication artifacts. Run all with:

```bash
just build-all-artifacts
```
