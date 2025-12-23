"""
amino_acid_mapping app

Utilities to analyze and visualize the mapping from three-letter amino acid
identifiers to one-letter codes.

This package provides:
- plot_amino_acid_distribution: Create colored bar charts distinguishing canonical
  from non-canonical amino acids, and optionally export a CSV summary.
- generate_latex_table: Produce a publication-ready LaTeX longtable from the
  mapping CSV.

Typical usage is via the CLI commands that wrap these functions. They accept a
CSV input and write outputs (PNG/PDF/TEX/CSV) to a chosen output directory.
"""

from typing import Set

# Canonical one-letter amino acid codes (standard 20)
DEFAULT_CANONICAL_AMINO_ACIDS: Set[str] = {
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
}

try:
    # These will be provided by the sibling module implementing the logic.
    from .amino_acid_mapping import (
        plot_amino_acid_distribution,
        generate_latex_table,
    )
except Exception:
    # Keep the package importable even if the implementation module is missing
    # during incremental development. Provide stubs that inform the user.
    def plot_amino_acid_distribution(*args, **kwargs):  # type: ignore[empty-body]
        raise RuntimeError(
            "plot_amino_acid_distribution is not available. "
            "Ensure 'amino_acid_mapping.py' is present and importable."
        )

    def generate_latex_table(*args, **kwargs):  # type: ignore[empty-body]
        raise RuntimeError(
            "generate_latex_table is not available. "
            "Ensure 'amino_acid_mapping.py' is present and importable."
        )

__all__ = [
    "DEFAULT_CANONICAL_AMINO_ACIDS",
    "plot_amino_acid_distribution",
    "generate_latex_table",
]
