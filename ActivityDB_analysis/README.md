# ActivityDB Analysis

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Companion Python package for the C++ tool **ActivityFinder** developed in the Rarey lab at the ZBH (Center for Bioinformatics) in Hamburg, Germany. This package generates the figures and tables for the ActivityFinder publication.

## Table of Contents

- [Publication Artifacts](#publication-artifacts)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
  - [Application Configuration](#application-configuration)
  - [SSH Tunnel Setup](#ssh-tunnel-setup)
  - [Environment Variables](#environment-variables)
- [Generating Artifacts](#generating-artifacts)
  - [Build All Artifacts](#build-all-artifacts)
  - [Individual Artifact Recipes](#individual-artifact-recipes)
- [Development](#development)
- [Project Structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Contact](#contact)

## Publication Artifacts

This repository generates 6 specific artifacts for the ActivityFinder publication:

| Artifact | Description | Recipe |
|----------|-------------|--------|
| `amino_acid_distribution_colored` | Amino acid identifier mapping bar chart (from PDB CCD data) | `just aa-plot` |
| `mutation_distribution_plot` | Variant mutation distribution bar plot | `just plot-variant` |
| `variants_per_het_code_histogram` | HET code variants histogram | `just plot-het` |
| `target_matching_interval_stacked_bars_identity_only` | Sequence mapping stacked bar plot | `just sequence-mapping-single-mode` |
| `ligand_match_levels.tex` | PDB-ChEMBL ligand mapping confidence level table (LaTeX) | `just ligand-match-level-counts` |
| `pair_identity_distribution_activity_combined_panel` | Pair identity distribution 2×2 panel plot | `just pair-identity-distribution-activity-plots` |

## Requirements

- Python 3.11 or higher
- PostgreSQL database access (ActivityDB and ChEMBL)
- SSH access to the database server (for remote connections)
- [uv](https://github.com/astral-sh/uv) (recommended) or pip

## Installation

### Using uv (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd ActivityDB_analysis

# Install dependencies
uv sync

# Verify installation
uv run actfind --help
```

### Using pip

```bash
# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in editable mode
pip install -e .

# With development dependencies
pip install -e ".[dev]"
```

## Configuration

### Application Configuration

The package uses `appconfig.toml` for database connection settings:

```toml
[activity_db]
name = "activitydb_test"
user = "activityuser_test"
host = "127.0.0.1"
port = 5433  # Local tunnel port
schema = "public"

[chembl_db]
name = "chembl"
user = "chembl_user"
host = "127.0.0.1"
port = 5433  # Local tunnel port
schema = "public"

[ssh]
host_alias = "your_ssh_alias"
remote_bind_host = "your.database.server.example.com"
remote_bind_port = 5432
local_bind_host = "127.0.0.1"
local_bind_port = 5433

[tool_version]
activity_finder = "0.31.1"
```

### SSH Tunnel Setup

The CLI automatically manages SSH tunnels to connect to the remote database. Configure your SSH client by adding an entry to `~/.ssh/config`:

```ssh
Host your_ssh_alias
    HostName your.database.server.example.com
    User your_username
    Port 22
    IdentityFile ~/.ssh/your_key
    # Optional: Add ProxyJump if needed
    # ProxyJump jump_host
    
    # Optional: Keep connection alive
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

**Tunnel Management Modes:**

| Flag | Mode | Description |
|------|------|-------------|
| `--manage-tunnel` | auto (default) | Start tunnel if needed, preserve external tunnels |
| `--no-tunnel` | skip | No tunnel management (for local development) |
| `--tunnel-mode required` | required | Fail if tunnel cannot be established |
| `--tunnel-mode optional` | optional | Continue without tunnel if unavailable |

Check tunnel status:

```bash
uv run actfind tunnel-status
```

### Environment Variables

Database passwords must be set via environment variables (never stored in config files):

```bash
# ActivityDB password (based on username suffix in config)
export ACTIVITY_DB_PRODUCTION_PWD="your_password"  # for *_production users
export ACTIVITY_DB_DEV_PWD="your_password"         # for *_dev users
export ACTIVITY_DB_TEST_PWD="your_password"        # for *_test users

# ChEMBL database password
export CHEMBL_PWD="your_password"
```

Add these to your shell profile (`~/.bashrc`, `~/.zshrc`, etc.) for persistence.

## Generating Artifacts

### Build All Artifacts

To generate all publication artifacts at once:

```bash
just build-all-artifacts
```

This creates all figures and tables in `./output/figures/` and `./output/tables/`.

### Individual Artifact Recipes

#### 1. Amino Acid Distribution Plot

Creates a colored bar chart showing amino acid identifier mappings from the PDB Chemical Component Dictionary.

```bash
just aa-plot <csv_file> [output_dir] [filename]
```

**Note:** Requires an external CSV file from the PDB CCD (not generated from the database). See `just help` for CSV format details.

#### 2. Mutation Distribution Plot

Creates a bar plot showing variant mutation distributions.

```bash
just plot-variant
```

#### 3. HET Code Variants Histogram

Creates a histogram of variants per HET code.

```bash
just plot-het
```

#### 4. Sequence Mapping Stacked Bars

Creates stacked bar plots showing target matching intervals. This recipe also creates/refreshes the required database materialized views.

```bash
just sequence-mapping-single-mode [coverage_mode] [output_dir] [config_file]
```

Default coverage mode is `identity-only`.

#### 5. Ligand Match Level Counts Table

Generates a LaTeX table of PDB-ChEMBL ligand mapping confidence levels.

```bash
just ligand-match-level-counts [output_dir]
```

#### 6. Pair Identity Distribution Panel

Creates pair identity distribution plots including a combined 2×2 panel.

```bash
# First, pull the data
just pair-identity-distribution-activity-data-pull [output_dir]

# Then create the plots
just pair-identity-distribution-activity-plots [figures_dir] [tables_dir]
```

### Utility Recipes

```bash
just check      # Check dependencies and configuration
just clean      # Clean generated files
just config     # Show current configuration
just help       # Show detailed usage information
just quickstart # Quick setup for new users
```

## Development

### Development Setup

```bash
# Install with dev dependencies
uv sync --extra dev

# Or with pip
pip install -e ".[dev]"
```

### Running Tests

```bash
uv run pytest
uv run pytest --cov=activity_db_analysis --cov-report=html
```

### Code Style

```bash
uv run black activity_db_analysis/
uv run ruff check activity_db_analysis/
```

## Project Structure

```
ActivityDB_analysis/
├── activity_db_analysis/           # Main package
│   ├── cli.py                      # Main CLI entry point
│   ├── apps/                       # Analysis applications
│   │   ├── amino_acid_mapping/     # Amino acid distribution plots
│   │   ├── het_code_plot/          # HET code variants analysis
│   │   ├── variant_mutation_plot/  # Variant mutation analysis
│   │   ├── ligand_mapping/         # PDB-ChEMBL ligand mapping
│   │   ├── sequence_mapping/       # Sequence mapping reports
│   │   └── sql_artifacts/          # SQL artifact management
│   ├── db/                         # Database layer
│   └── siteconfig/                 # Configuration
├── artifacts/                      # Generated publication artifacts
├── output/                         # Default output directory
├── appconfig.toml                  # Database configuration
├── desired_artifacts.txt           # List of target artifacts
├── justfile                        # Task runner recipes
├── pyproject.toml                  # Package metadata
└── uv.lock                         # Locked dependencies
```

## Troubleshooting

### SSH Tunnel Issues

**Port already in use:**
```bash
# Check what's using the tunnel port
lsof -i :5433

# Or use --no-tunnel for local testing
uv run actfind --no-tunnel <command>
```

**SSH connection fails:**
```bash
# Test SSH connection directly
ssh -v your_ssh_alias

# Verify SSH config
cat ~/.ssh/config | grep -A5 your_ssh_alias
```

### Database Connection Issues

**Authentication failed:**
```bash
# Verify environment variables are set
echo $ACTIVITY_DB_TEST_PWD
echo $CHEMBL_PWD

# Check with verbose logging
uv run actfind --logging-level DEBUG <command>
```

**Missing materialized views:**
```bash
# Check artifact status
uv run actfind sql-manager status

# Create missing artifacts
uv run actfind sql-manager create <artifact-name>
```

### General Debugging

```bash
# Enable trace-level logging
uv run actfind --logging-level TRACE <command>

# Verify Python environment
uv run python -c "import activity_db_analysis; print(activity_db_analysis.__file__)"
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please open a GitHub issue or contact the maintainers.

---

*This package is part of the ActivityFinder project for mapping PDB ligands to ChEMBL bioactivity data.*