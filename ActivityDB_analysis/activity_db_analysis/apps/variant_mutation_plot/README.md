# Variant Mutation Analysis

This module provides functionality to analyze and visualize the distribution of mutation counts in ChEMBL variant sequences. It generates the **mutation_distribution_plot** publication artifact.

## Quick Start

```bash
# Using justfile (recommended)
just plot-variant

# Using CLI directly
uv run actfind plot-variant-mutations --output-dir ./output
```

## Overview

The analysis examines how many mutations are present in each variant sequence in the ChEMBL database. This helps understand the distribution of genetic variants and their complexity in terms of mutation counts.

## ChEMBL Variant Sequences Explained

ChEMBL stores protein target sequences, including natural variants and mutations that differ from canonical reference sequences.

- **Variant sequences table**: The `variant_sequences` table contains protein sequences that differ from the canonical/reference sequence for a given target
- **Mutations recorded**: Amino acid substitutions are recorded in these variant sequences (e.g., "V600E" indicates valine→glutamic acid at position 600)
- **Purpose of this analysis**: Understanding how many mutations typically occur per variant helps characterize the complexity of genetic variation in drug targets
- **Pharmacogenomics relevance**: Drug response can vary significantly based on target mutations—a drug effective against wild-type may be ineffective against a mutant, or vice versa

## Usage

### Using Justfile (Recommended)

```bash
# Create variant mutation plots with default settings
just plot-variant
```

### Using the CLI Command

```bash
uv run actfind plot-variant-mutations \
    --output-dir ./output \
    --output-prefix mutation_analysis
```

By default, the command uses an embedded SQL query. The command uses the ChemblInteraction module to connect to the ChEMBL PostgreSQL database, which provides proper connection handling and supports SSH tunneling if configured.

**LaTeX Table Requirements**: The generated `.tex` files require these packages in your LaTeX document:
```latex
\usepackage{booktabs}   % For professional table formatting
\usepackage{siunitx}    % For enhanced number formatting
```

### Command Options

- `--output-dir`: Directory where output files will be saved (default: `./output`)
- `--output-prefix`: Base name for output files (default: `mutation_distribution`)
- `--latex/--no-latex`: Generate LaTeX table (default: true)
- `--create-plots/--no-create-plots`: Create bar plots (default: true)

## Configuration

Database connection parameters are specified in `appconfig.toml`:

```toml
[chembl_db]
name = "chembl"
user = "chembl_user"
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

Database password is set via environment variable:
```bash
export CHEMBL_PWD="your_password"
```

## Output Files

The command generates:

| File | Description |
|------|-------------|
| `{prefix}_data.csv` | Raw data (num_mutations, num_sequences, percentage) |
| `{prefix}_plot.png/.pdf` | Bar plot with value labels |
| `{prefix}_table.tex` | Publication-ready LaTeX table |

## SQL Query

The default SQL query is embedded in `variant_mutation_sql.py`. It extracts mutation counts from ChEMBL's `variant_sequences` table, handling NULL/empty fields and multiple mutations separated by semicolons or commas.

## Interpreting Results

The analysis typically shows:
- Most sequences have 0-1 mutations (wild-type or single variants)
- Decreasing frequency with increasing mutation counts
- A long tail of sequences with many mutations

## Part of Publication Artifacts

This app generates `mutation_distribution_plot`, one of the 6 publication artifacts for the ActivityFinder paper. Run all artifacts with:

```bash
just build-all-artifacts
```