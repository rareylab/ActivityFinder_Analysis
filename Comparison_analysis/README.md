# ActivityFinder and BioChemGraph comparison

This project analyzes ActivityFinder results and compares them with BioChemGraph,
including the effects of target and molecule matching. It also contains supporting
sequence, mutation, database, and visualization analyses. The scripts use
ActivityDB and ChEMBL data in addition to input exports from both systems.

## Project layout

```text
Comparison_analysis/
├── pyproject.toml
└── src/
    ├── definitions.py               # Input/output paths and local DB settings
    ├── extract_json_from_activityfinder.py
    ├── extract_json_from_biochemgraph.py
    ├── compare_jsons.py
    ├── create_sankey_diagram.py
    ├── analyse_notfound_targets.py
    ├── sequence_analysis/
    ├── mutation_analysis/
    ├── molecule_analysis/
    └── utils/                        # Database, data, and visualization helpers
```

`src/ActivityFinderRestApi.ipynb` demonstrates querying the ActivityFinder REST
API for individual PDB entries. Other analysis scripts are collected directly
under `src/` and in the topical subdirectories above.

## Setup and configuration

The project dependencies and optional development dependencies are listed in
[`pyproject.toml`](pyproject.toml). From this directory, create an environment and
install the project (including the optional development tools if needed):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

Before running scripts, review `src/definitions.py`. It defines the `input_files/`
and `state_presubmission/` locations, required input filenames, and local database
connection settings. The data exports and database services are not included in
this project. As checked in, `src/definitions.py` has blank environment-specific
assignments and cannot be imported as-is; complete the local configuration before
running scripts. Do not commit credentials. Some analyses require additional
input files beyond those checked by `definitions.py`.

## Typical comparison workflow

Run from this directory. The scripts import modules both from `src/` and from the
project root, so add both to `PYTHONPATH`:

```bash
export PYTHONPATH=".:src${PYTHONPATH:+:$PYTHONPATH}"
```

1. Place the required ActivityFinder exports and BioChemGraph report in the
   configured `input_files/` directory.
2. Extract the BioChemGraph and ActivityFinder records into JSON files:

   ```bash
   python src/extract_json_from_biochemgraph.py
   python src/extract_json_from_activityfinder.py
   ```

3. Run `src/compare_jsons.py` to compare the extracted records. This analysis also
   queries ActivityDB and ChEMBL, so working database access and local connection
   settings are required.
4. Run `src/create_sankey_diagram.py` for the more detailed breakdown of results
   and reasons for differences. It depends on generated data and database access
   as well.

The extraction scripts can take time on the full datasets. The other scripts are
individual analyses and may depend on outputs created by earlier steps; check each
script and its input paths before running it.

## Script reference

These descriptions summarize the current code. Many files perform work as soon as
they are run (or imported), and are analysis scripts rather than command-line tools
with configurable arguments. Check their paths and inputs before execution.

### Main analysis scripts (`src/`)

- **`definitions.py`** — Centralizes paths, matching-confidence levels, and database settings;
  checks required inputs and creates output directories.
- **`extract_json_from_activityfinder.py`** — Reads ActivityFinder activity and target
  exports and writes JSON summaries at several molecule-matching and sequence-identity
  thresholds. Also writes a quality summary and processes a historical activity export.
- **`extract_json_from_biochemgraph.py`** — Reads the BioChemGraph TSV report and writes
  JSON summaries of PDB–target–molecule records and target–UniProt mappings.
- **`compare_jsons.py`** — Compares extracted ActivityFinder and BioChemGraph records,
  including PDB, target, and molecule overlaps. Also uses ActivityDB/ChEMBL data and
  the ChEMBL API.
- **`create_sankey_diagram.py`** — Breaks down shared and unmatched records and their
  likely causes, then produces Sankey and Sunburst visualizations using extracted
  data and database-backed helpers.
- **`analyse_notfound_targets.py`** — Investigates targets not found by ActivityFinder,
  joining ligand, database, and sequence information to characterize missing cases.
- **`cheatsheet.py`** — Calculates overlap statistics and creates Venn diagrams for
  PDBs, ChEMBL ligands and targets, HET codes, and combined records. It contains
  machine-specific paths and data assumptions.
- **`check_inconsistencies.py`** — Compares PDB–ligand–target–HET-code combinations in
  an enriched presence file with the raw ActivityFinder export; writes CSVs for
  records unique to either source and their intersection.
- **`2d_gridplots.py`** — Creates faceted tile plots of molecule-matching categories
  against sequence-identity intervals for PDB, component, or variant data. Expects
  prepared CSV inputs.
- **`find_moldiff_examples_db.py`** — Explores molecule-matching differences between
  ActivityFinder and BioChemGraph, including overlap, HET-code, and example cases.
- **`paper_molmatching_plots.py`** — Generates selected molecule-comparison figures
  for paper examples. The script uses fixed structures and an environment-specific
  output path.
- **`get_assays_with_variant_sequence.py`** — Finds assays with variant sequences,
  associates them with PDB records, and compares mutation annotations between
  component and variant sequences. Produces intermediate and difference datasets.

> **Visualization tip:** The Sunburst view from `create_sankey_diagram.py` is especially helpful for exploring the breakdown.

### Sequence analysis (`src/sequence_analysis/`)

- **`sequence_matching.py`** — Aligns PDB SEQRES chains with UniProt sequences for
  BioChemGraph PDB–UniProt pairs, selects the best identity, and saves the results.
  Retrieves UniProt sequences online and uses a local SEQRES input.
- **`sequence_identity_plots.py`** — Reads the sequence-identity results and writes histogram plots, including a log-scale version.
- **`find_slightly_different_target_from_biochem.py`** — Finds ActivityFinder records
  with exact target-sequence identity and molecule matching where the PDB is absent
  from BioChemGraph's target set; saves all and Gold-level subsets.

### Molecule analysis (`src/molecule_analysis/`)

- **`generate_moldifference_example_plots.py`** — Renders hard-coded ChEMBL/PDB InChI
  examples as side-by-side molecule images, illustrating matching differences and
  structural discrepancies.

### Mutation analysis (`src/mutation_analysis/`)

- **`extract_mutation_data.py`** — Queries mutation/alignment details for a PDB entry,
  split into component- and variant-sequence records.
- **`analyse_mutation_data.py`** — Converts mutation CSVs into nested dictionaries;
  when run directly, combines files from the configured mutation-data directory into
  a JSON output.
- **`variant_data_linker.py`** — Provides lookups linking assays, targets, HET codes,
  ligands, and component/variant sequences; its analysis investigates variant-related
  matches.
- **`remove_false_positives_from_variant_sequence.py`** — Filters generated
  variant-mutation differences to remove gap/unknown-residue cases and writes a
  cleaned CSV.
- **`analyse_component_only_mutations.py`** — Summarizes cleaned mutation-difference
  data, writes mutation-count histograms, and builds an assay/PDB/ligand mapping from
  the raw activity export.

### Shared utilities (`src/utils/`)

These modules support the analysis scripts and are not intended as standalone
workflows:

- **`db.py`** — ActivityDB and ChEMBL query helpers, including cached query results
  and activity lookups.
- **`getter.py`** — Retrieves commonly used ActivityDB/ChEMBL tables and mappings,
  reusing saved output files where available.
- **`general.py`** — Shared helpers for splitting combinations, retrieving molecule
  structures and activities, and enriching molecule-match results.
- **`visualisation.py`** — Venn/SuperVenn charts, molecule depictions, and
  molecule-array visualization helpers.

## Inputs and outputs

- `input_files/` contains the source exports. `definitions.py` currently expects
  ActivityFinder activity, target, and enriched-activity CSVs, plus BioChemGraph
  full and simple TSV reports. Some scripts also read other files from this
  directory, such as ligand and sequence data. Published input data are available
  in the [StrAcTable Supporting Data](https://www.fdr.uni-hamburg.de/record/18244).
- `state_presubmission/` is the configured output root. The definitions module
  creates its output subdirectories, including `output_files/`, `plots/`,
  `automated_analysis_outputs/`, `molecule_difference_data/`, and
  `mutation_data/` when it is imported.
- Extracted JSON, intermediate tables, plots, and other analysis products are
  written under those output directories. They are generated artifacts, not
  source inputs.

The published Supporting Data covers the released inputs; some analyses may also
require additional local exports, database access, and environment-specific
configuration. Reproducibility depends on these inputs and services being available.
