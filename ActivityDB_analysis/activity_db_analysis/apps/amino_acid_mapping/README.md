# Amino Acid Mapping

Tools to analyze and visualize the mapping from three-letter amino acid identifiers to one-letter codes. This app generates the **amino_acid_distribution_colored** publication artifact.

No database connection is required—this app works with external CSV data from the PDB Chemical Component Dictionary.

## Quick Start

```bash
# Using justfile (recommended)
just aa-plot ./input/amino_acid_mapping.csv ./artifacts amino_acid_distribution_colored

# Using CLI directly
uv run actfind plot-amino-acid-mapping \
    --input /path/to/amino_acid_mapping.csv \
    --output-dir ./artifacts \
    --output-prefix amino_acid_distribution_colored \
    --pdf --no-show
```

## Features

- Colored bar chart (canonical vs non-canonical amino acids)
- High-resolution PNG and vector PDF outputs
- Optional CSV summary for the bar chart
- LaTeX longtable with chemistry-aware formatting

## Input Data Source

The input CSV is obtained from the **PDB Chemical Component Dictionary (CCD)**, not from the database.

**Where to obtain:**
- Download CCD from [wwPDB](https://www.wwpdb.org/data/ccd)
- Filter for amino acid-type components (`_chem_comp.type` contains "PEPTIDE")
- Extract columns: `mon_id`, `one_letter_code`, `name`

**Expected CSV format (semicolon-delimited):**
```csv
"Three Letter Code";"One Letter  Code";"Name"
"ALA";"A";"ALANINE"
"MSE";"M";"SELENOMETHIONINE"
```

## CLI Usage

```bash
uv run actfind plot-amino-acid-mapping \
    --input ./amino_acid_mapping.csv \
    --output-dir ./artifacts \
    --output-prefix amino_acid_distribution_colored \
    --delimiter ";" \
    --pdf \
    --no-show
```

### Options

| Option | Description |
|--------|-------------|
| `--input` | Path to input CSV file (required) |
| `--output-dir` | Output directory (default: `./output`) |
| `--output-prefix` | Base filename (default: `amino_acid_distribution`) |
| `--delimiter` | CSV delimiter (default: `;`) |
| `--pdf` | Generate PDF output |
| `--no-show` | Don't display plot interactively |

## Output Files

| File | Description |
|------|-------------|
| `{prefix}.png` | Bar chart (PNG) |
| `{prefix}.pdf` | Bar chart (PDF) |
| `{prefix}_data.csv` | Summary counts (optional) |

## Understanding the Output

- **X-axis**: One-letter amino acid codes (A, C, D, ... Y)
- **Y-axis**: Count of distinct three-letter codes mapping to each one-letter code
- **Blue bars**: Canonical amino acids (standard 20)
- **Orange bars**: Non-canonical amino acids (modified, rare, synthetic)

## LaTeX Table Generation

```python
from activity_db_analysis.apps.amino_acid_mapping import generate_latex_table
from pathlib import Path

tex_path = generate_latex_table(
    input_csv=Path("./amino_acid_mapping.csv"),
    output_path=Path("./amino_acid_mapping_table.tex"),
    delimiter=";",
)
```

LaTeX prerequisites:
```latex
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{siunitx}
```

## Part of Publication Artifacts

This app generates `amino_acid_distribution_colored`, one of the 6 publication artifacts. Run all with:

```bash
just build-all-artifacts
```

**Note:** Unlike other artifacts, this requires an external CSV file from the PDB CCD.
