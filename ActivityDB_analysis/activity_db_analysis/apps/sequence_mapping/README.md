# Sequence Mapping Reporting & Metrics

This app generates the **target_matching_interval_stacked_bars_identity_only** publication artifact—a key visualization showing PDB-ChEMBL sequence alignment coverage.

---

## 1. What This App Does

Generates publication-ready figures showing PDB-ChEMBL sequence alignment coverage. It quantifies how many unique PDB sequences can be mapped to ChEMBL target sequences at different identity thresholds (≥80%, ≥85%, ≥90%, ≥95%, ≥100%). Output includes stacked bar charts, CSV data files, and LaTeX tables.

---

## 2. Quick Start

```bash
# Using justfile (recommended) - full pipeline for identity-only mode
just sequence-mapping-single-mode identity-only ./artifacts

# Or run individual steps with CLI
uv run actfind sql-manager create mv_alignment_sequence_metrics
uv run actfind sequence-mapper generate-all --output-dir ./output --coverage-modes identity-only
uv run actfind sequence-mapper plot --input ./output/target_matches_*.csv --output ./figures/
```

Build all publication artifacts including sequence mapping:
```bash
just build-all-artifacts
```

---

## 3. Sample Output

| Output Type | Description | Example Filename |
|-------------|-------------|------------------|
| **Stacked Bar Chart** | Sequence pairs by identity interval | `target_matching_interval_stacked_bars_identity_only.png` |
| **CSV Data** | Detailed metrics per coverage mode | `target_matches_by_sequences_sequence_identity_ranges_identity_only.csv` |
| **LaTeX Tables** | Publication-ready tables | `target_matching_interval_counts_identity_only.tex` |

---

## 4. CLI Reference

### Generate Reports

```bash
# Sequence-level reports
uv run actfind extract-pdb-chembl-sequence-mappings ./output --coverage-mode identity-only

# Chain-level reports
uv run actfind extract-pdb-chembl-sequence-mappings ./output --coverage-mode identity-only --chain-level
```

### Create Plots

```bash
uv run actfind sequence-mapper plot \
  --input ./output/target_matches_by_sequences_sequence_identity_ranges_identity_only.csv \
  --output ./figures/target_matching_interval_stacked_bars_identity_only.png \
  --coverage-mode identity-only
```

---

## 5. Objectives

We quantify and visualize how many unique PDB sequence ↔ ChEMBL sequence pairings (“unique PDB–ChEMBL sequence pairs”) exist at different sequence identity similarity granularities.

Breakdowns:

- ChEMBL sequence classification: `COMPONENT` vs `VARIANT`
- PDB sequence type: `ATOMSEQ`, `SEQRES`, `SEQADV` (the three internal sequence representations)
- Identity thresholds (cumulative “right‑bound”): ≥100, ≥95, ≥90, ≥85, ≥80
- Disjoint identity intervals (mutually exclusive buckets): `[100]`, `[95,100)`, `[90,95)`, `[85,90)`, `[80,85)`
- Coverage mode tiers (filters on percent identity + coverage): `full`, `no-hit`, `identity-only`
- Optional: chain‑level intersection cohort vs global (sequence‑level) cohort

Outputs support:
1. Analytical tables (wide CSV + LaTeX).
2. “Simplified” publication tables (counts + column‑wise proportions).
3. Stacked bar plots (absolute + normalized composition).

---

## 7. Data Sources (Materialized Views)

| MV | Level | Purpose |
|----|-------|---------|
| `mv_alignment_sequence_metrics` | Sequence-level | Pre-aggregated right-bound + interval metrics across tiers |
| `mv_alignment_sequence_metrics_chain` | Chain-level | Same shape; restricted to chains having **all three** sequence types (intersection) |
| `mv_sequence_totals_chain` | Chain-level (aux) | Provides total & mapped chain counts (used for possible future diagnostics) |

### 2.1 Coverage Modes (Tier Mapping)

| Coverage Mode (CLI) | Tier Label (MV) | Semantics |
|---------------------|-----------------|-----------|
| `identity-only` | `identity-only` | Only identity ≥ 80; no coverage filters |
| `no-hit` | `no-hit` | Identity ≥ 80; applies query coverage ≥ 80; hit coverage not required |
| `full` | `full` | Identity ≥ 80; query coverage ≥ 80; hit coverage ≥ 80 |

The tier union is materialized; report generation filters by `tier = <mapped>`.

---

## 8. Chain-Level Semantics

Chain-level reporting restricts the universe to `(pid, cid)` (PDB entry + chain) pairs that possess **all three** sequence types (ATOMSEQ, SEQRES, SEQADV). This intersection ensures cross-type comparability within a consistent structural cohort.

Key implications:

- Counts are (necessarily) ≤ sequence-level counts for RAW and distinct sequences per threshold / interval.
- DENSITY / OCC / COND_OCC may increase or decrease relative to sequence-level depending on enrichment vs exclusion.
- Filenames are prefixed with `chain_level_` for chain-level outputs.
- The chain-level MV has its own indexes to keep retrieval fast (`ix_mvasmc_*`).

---

## 9. Metric Definitions

All metrics are computed after best-alignment selection per (query_sid, hit_sid) in the tier views feeding the MV.

Per identity **interval row** or **right-bound threshold row**:

| Column Suffix | Name | Definition | Range / Notes |
|---------------|------|------------|---------------|
| `_RAW` | Raw pair count | Number of qualifying alignment-derived unique PDB–ChEMBL sequence pair rows (multi-mapping retained) | Integer ≥ 0 |
| `_DENSITY` | Pair density | `_RAW / total_sequences` | Float ∈ [0, ∞) (often ≤ 1, but can exceed 1 because a single PDB sequence can map to >1 target sequence in different contexts) |
| `_OCC` | Distinct occupancy | `distinct_sequences_with_pair / total_sequences` | Float ∈ [0, 1] |
| `_COND_OCC` | Conditional occupancy | `distinct_sequences_with_pair / mapped_sequences_any` | Float ∈ [0, 1]; denominator is the number of sequences with at least one match at the **base floor (≥80)** across any interval/threshold for that (chembl_type, sequence type) |

Additional rows:

- `TOTAL_SEQUENCES` (placed in *_RAW columns for discoverability): number of distinct PDB sequences (or intersection chain sequences in chain-level MV) available for that PDB sequence type, regardless of identity interval.
- `MAPPED_SEQUENCES` (in *_RAW columns): number of sequences that have **at least one** match at identity ≥ 80 (foundation for `mapped_sequences_any` denominator in `_COND_OCC`).

### 9.1 Distinguishing Denominators

| Metric | Numerator | Denominator | Intent |
|--------|-----------|-------------|--------|
| DENSITY | Raw pair rows in bucket | Total sequences | Mapping multiplicity magnitude |
| OCC | Distinct sequences in bucket | Total sequences | Coverage breadth across base population |
| COND_OCC | Distinct sequences in bucket | Mapped sequences (≥80 at any bucket) | Normalizes out sequences that never participate, highlighting distribution within the mapped sub-cohort |

---

## 10. Right-Bound vs Interval Logic

- **Right-Bound Thresholds** (`>=100`, `>=95`, ... , `>=80`): cumulative perspective. A sequence pair counted at ≥95 also appears in ≥90, ≥85, ≥80.
- **Disjoint Intervals**: mutually exclusive partition of the ≥80 universe.

The MV stores both forms:
- Cumulative flag aggregations: `raw_ge_100`, `distinct_ge_100`, ... `raw_ge_80`, `distinct_ge_80`
- Interval rows: `identity_interval`, `raw_pairs_interval`, `distinct_sequences_interval`

Python transformation builds DataFrames for each side (COMPONENT / VARIANT) and PDB sequence type.

---

## 11. Coverage Mode Influence

Coverage modes alter the **input** set of best alignments passing coverage filters; all subsequent per-threshold / interval metrics rest atop that filtered universe. Metric definitions themselves do not change.

Typical analytic flows:
1. Run `full` (strict) for high-confidence narrative.
2. Compare `no-hit` vs `full` for incremental gain from including low hit coverage alignments.
3. Compare `identity-only` vs `no-hit` for breadth vs coverage quality trade-offs.

---

## 12. Output Artifact Naming

Pattern components:

| Component | Value |
|-----------|-------|
| Optional prefix | `chain_level_` (only for chain-level) |
| Base stem (right-bound) | `target_matches_by_sequences_right_bound` |
| Base stem (interval) | `target_matches_by_sequences_sequence_identity_ranges` |
| Coverage suffix | `_full` `_no_hit` `_identity_only` (underscored) |
| Simplified tables | `target_matching_interval_counts_*.tex`, `target_matching_interval_proportions_*.tex` |
| Extensions | `.csv` / `.tex` |

Examples:
```
target_matches_by_sequences_right_bound_full.csv
chain_level_target_matches_by_sequences_sequence_identity_ranges_no_hit.csv
target_matching_interval_counts_identity_only.tex
chain_level_target_matching_interval_proportions_full_chain_level.tex
```

---

## 13. Simplified Publication Tables

Two LaTeX tables per coverage mode (and per chain-level flag):
1. Counts (`*_interval_counts_...`): RAW counts only (interval rows + total).
2. Proportions (`*_interval_proportions_...`): Column-wise percentage composition within each (side, sequence type) summing to 100.00%.

They:
- Power the normalized stacked bar figure.
- Provide compact thesis/table-ready summaries.
- Use `siunitx` formatting macros (assumed in the broader document).

---

## 14. Plot Generation

Script: `plot_sequence_matching_interval_stacked_bars.py`

Input: Interval CSV (sequence-level or chain-level) with coverage suffix.

Figure layout (2 × 2):
- Rows: `Component`, `Variant`
- Columns: Absolute RAW stacked counts | Normalized composition
- Stacks: Identity intervals (color-coded)

If chain-level CSV is used, title/legend should be manually set to reflect chain-level context (script does not infer automatically).

---

## 15. Parity & Regression Safety

Script: `sequence_mapping_parity_check.py`

Purpose:
- Recompute legacy “best hit + coverage predicate” logic (for *sequence-level*) and compare to the tier view / MV-driven approach.
- Validates `pair_count_raw` and `distinct_sequences` at both right-bound and interval levels across coverage modes.

Chain-Level parity is not implemented (legacy did not have a dedicated chain-level MV). Future extension can add query filtering by intersection chains to parity logic if needed.

Recommended workflow:
1. Refresh MVs.
2. Run parity checker for `full`, `no-hit`, `identity-only`.
3. Investigate any diffs before proceeding to publication outputs.

---

## 16. Indexing (Chain-Level MV)

To safeguard performance parity for interactive reporting:
- `ix_mvasmc_tier_chembl_type`: Filters combining tier + side.
- `ix_mvasmc_dim_interval`: Full dimensional + interval predicate.
- `ix_mvasmc_dim`: Tier + sequence type splits.

These mirror common WHERE clauses used when isolating tiers and building Python DataFrames.

---

## 17. Adding / Modifying Metrics

To add a new metric (e.g., coverage span class):
1. Extend MV CTEs with the additional derived column(s).
2. Add columns to SELECT projection in both sequence-level and chain-level MVs (keep schemas parallel).
3. Update transformation helpers in `sequence_mapping.py` to include new metric suffix groups (if exposing in CSV/TeX).
4. Adjust simplified table logic only if the new metric should appear (currently only RAW is used).
5. Add tests to assert presence of new column tokens in MV SQL.

---

## 18. Known Limitations / Future Enhancements

| Area | Current State | Possible Enhancement |
|------|---------------|----------------------|
| Parity (chain-level) | Not implemented | Add chain-level parity variant |
| Multi-mode batch | One mode per invocation (unless CLI helper loops) | Provide `--all-coverage-modes` switch |
| Identity floor configurability | Fixed at ≥80 (MV design) | Parameterize MV (would increase complexity) |
| Plot title automation | Manual | Infer chain-level & coverage from filename |
| Historical trend | Not tracked | Append date-stamped snapshots or log table |
| Additional metrics | Identity-focused | Add coverage quantile distributions / per-target stratifications |

---

## 19. Parameter Caveat

`min_identity_floor` argument in the reporting function is **inert**: the MVs are already anchored at identity ≥ 80. Passing a different value will not change results; retained only for backward compatibility. A warning is recommended if non‑80 is supplied.

---

## 20. Glossary

| Term | Meaning |
|------|--------|
| PDB sequence type | Representation flavor (ATOMSEQ: atomic coordinates derived, SEQRES: author-declared residues, SEQADV: SEQRES plus SEQADV-corrected) |
| Component / Variant | ChEMBL target sequence classification (canonical component vs curated variant) |
| Pair (sequence-level) | PDB sequence (query_sid) ↔ ChEMBL sequence (hit_sid) best alignment |
| Chain-level intersection | (pid, cid) chain having all three sequence types present |
| RAW pair count | Count of per-pair alignment records satisfying tier + bucket constraints |
| Occupancy | Distinct PDB sequences with ≥1 qualifying pair divided by total sequences |
| Conditional occupancy | Distinct PDB sequences with ≥1 qualifying pair divided by those sequences that map at identity ≥80 at all |
| Density | Raw pair multiplicity per total sequence (exposes “fan-out”) |

---

## 21. Validation Checklist (Manual QA)

| Check | Expected |
|-------|----------|
| Distinct count monotonicity (right-bound) | `>=100 ≤ >=95 ≤ ... ≤ >=80` for RAW & distinct counts |
| Interval sum consistency | Sum of interval RAW counts == RAW at >=80 threshold |
| Sequence vs chain-level | Chain-level RAW counts ≤ sequence-level for same tier & bucket |
| Simplified proportions | Each column sums to 100.00 (allow minor rounding) |
| Coverage tiers | Identity-only ≥ no-hit ≥ full (RAW counts & distinct) |

---

## 22. Support & Extension Notes

When extending:
- Keep sequence-level and chain-level MV schemas synchronized.
- Update tests (add tokens for new columns).
- Maintain consistent naming patterns to avoid breaking plot scripts.
- Document any new denominators clearly here.

---

*This README is the authoritative description of the sequence mapping reporting layer and its metrics, complementing the SQL artifact registry documentation.* 

---

## 23. Justfile Recipe Reference

The primary recipe for sequence mapping is `sequence-mapping-single-mode`, which runs the full pipeline for a single coverage mode.

### Available Recipe

```bash
# Full pipeline: create MVs → extract data → generate plots
just sequence-mapping-single-mode [coverage_mode] [output_dir] [config_file]

# Examples:
just sequence-mapping-single-mode                           # identity-only, ./artifacts
just sequence-mapping-single-mode full ./output             # full coverage mode
just sequence-mapping-single-mode identity-only ./artifacts appconfig.toml
```

### What the Recipe Does

1. **Phase 1**: Drops and recreates MVs (`mv_alignment_sequence_metrics`, `mv_alignment_sequence_metrics_chain`)
2. **Phase 2**: Generates interval CSV data for the specified coverage mode
3. **Phase 3**: Creates PNG and PDF stacked bar plots

### Build All Artifacts

To build all publication artifacts including sequence mapping:

```bash
just build-all-artifacts
```

### Manual CLI Workflow

For finer control, use the CLI directly:

```bash
# 1. Create SQL artifacts
uv run actfind sql-manager create mv_alignment_sequence_metrics

# 2. Generate data for all coverage modes
uv run actfind sequence-mapper generate-all --output-dir ./output --coverage-modes identity-only

# 3. Create plots
uv run actfind sequence-mapper plot \
    --input ./output/target_matches_by_sequences_sequence_identity_ranges_identity_only.csv \
    --output ./figures/target_matching_interval_stacked_bars_identity_only.png
```

### File Naming

| Cohort | Example CSV | Example Plot |
|--------|-------------|--------------|
| Sequence-level | `target_matches_by_sequences_sequence_identity_ranges_identity_only.csv` | `target_matching_interval_stacked_bars_identity_only.png` |
| Chain-level | `chain_level_target_matches_*_identity_only.csv` | `target_matching_interval_stacked_bars_chain_level_identity_only.png` |

---

## 24. Creating / Refreshing SQL Artifacts

The sequence mapping reports depend on SQL artifacts managed by `sql-manager`.

### Quick Setup

```bash
# Create all artifacts (recommended for fresh setup)
uv run actfind sql-manager create-all

# Or create only sequence mapping artifacts
uv run actfind sql-manager create mv_alignment_sequence_metrics
uv run actfind sql-manager create mv_alignment_sequence_metrics_chain
```

### Refresh After Data Updates

```bash
uv run actfind sql-manager refresh mv_alignment_sequence_metrics
uv run actfind sql-manager refresh mv_alignment_sequence_metrics_chain
```

### Check Status

```

---
