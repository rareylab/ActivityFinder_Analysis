"""
Amino acid mapping utilities:
- Plot distribution of three-letter identifiers mapped to single-letter amino acid codes.
- Generate a LaTeX longtable from the mapping CSV.

The plotting function colors canonical (standard 20) vs non-canonical amino acids.
The LaTeX table escapes special characters and applies simple chemistry-aware formatting.

Expected CSV columns (semicolon separated by default):
- "Three Letter Code"
- "One Letter  Code" (note: source has two spaces, we normalize this)
- "Name"

Both functions are file-system side-effecting by design: they read an input CSV and write
outputs (PNG/PDF/TEX/CSV) to the given output directory/path.
"""

from __future__ import annotations

import re
import csv
import math
import pathlib
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from loguru import logger
from activity_db_analysis.plotting.palettes import CANONICAL_AMINO_COLORS


# Default canonical one-letter amino acid codes (standard 20)
DEFAULT_CANONICAL_AMINO_ACIDS: Set[str] = {
    "A", "R", "N", "D", "C",
    "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V",
}


# --- CSV parsing helpers ---------------------------------------------------------------------------

def _normalize_header(header: str) -> str:
    """
    Normalize CSV header names by:
    - Stripping surrounding whitespace
    - Collapsing multiple spaces to single spaces
    - Lowercasing
    """
    if header is None:
        return ""
    # Remove enclosing quotes and trim
    h = str(header).strip().strip('"').strip("'")
    # Collapse runs of whitespace
    h = re.sub(r"\s+", " ", h)
    return h.lower()


def _read_mapping_csv(
    csv_path: pathlib.Path,
    delimiter: str = ";",
) -> pd.DataFrame:
    """
    Read the amino acid mapping CSV and normalize to a canonical schema with columns:
    - three_letter (str)
    - one_letter (str)
    - name (str)
    """
    csv_path = pathlib.Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    logger.info(f"Reading amino acid mapping CSV: {csv_path}")
    df = pd.read_csv(csv_path, sep=delimiter, dtype=str, quoting=csv.QUOTE_MINIMAL, keep_default_na=False)

    # Normalize headers
    original_cols = list(df.columns)
    norm_cols = [_normalize_header(c) for c in df.columns]

    # Build a mapping from normalized to original
    col_map = dict(zip(norm_cols, original_cols))

    # Identify relevant columns (tolerate spacing differences)
    # Common variants observed:
    #   "Three Letter Code"
    #   "One Letter  Code"  <-- mind the double space
    #   "One Letter Code"
    expected_three = ("three letter code", "three-letter code", "3-letter code")
    expected_one = ("one letter code", "one letter  code", "1-letter code")
    expected_name = ("name", "description", "chemical name")

    def _find_col(candidates: Sequence[str]) -> Optional[str]:
        for cand in candidates:
            if cand in col_map:
                return col_map[cand]
        # Fallback: fuzzy contains match
        for nc, oc in zip(norm_cols, original_cols):
            for cand in candidates:
                if cand in nc:
                    return oc
        return None

    three_col = _find_col(expected_three)
    one_col = _find_col(expected_one)
    name_col = _find_col(expected_name)

    if not three_col or not one_col or not name_col:
        raise ValueError(
            f"Unable to detect expected columns in CSV.\n"
            f"Detected columns: {original_cols}\n"
            f"Expected something like: {expected_three}, {expected_one}, {expected_name}"
        )

    # Standardize to canonical schema and strip quotes/whitespace from values
    out = pd.DataFrame({
        "three_letter": df[three_col].astype(str).str.strip().str.strip('"').str.strip("'"),
        "one_letter": df[one_col].astype(str).str.strip().str.strip('"').str.strip("'"),
        "name": df[name_col].astype(str).str.strip().str.strip('"').str.strip("'"),
    })

    logger.debug(f"Loaded {len(out)} rows from mapping CSV.")
    return out


# --- Plotting --------------------------------------------------------------------------------------

def plot_amino_acid_distribution(
    input_csv: pathlib.Path,
    output_dir: pathlib.Path,
    output_filename_base: str = "amino_acid_distribution",
    delimiter: str = ";",
    write_csv_summary: bool = True,
    canonical_amino_acids: Optional[Set[str]] = None,
) -> List[pathlib.Path]:
    """
    Parse the amino acid mapping CSV and generate a histogram (bar chart) showing how
    many distinct three-letter identifiers map to each single-letter code.

    Args:
        input_csv: Path to the input CSV file.
        output_dir: Directory to write the output files (created if not existing).
        output_filename_base: Basename for the output files (no extension).
        delimiter: CSV delimiter (default ';').

        write_csv_summary: If True, write a summary CSV with counts and type.
        canonical_amino_acids: Optional override set of canonical one-letter codes.

    Returns:
        List of paths to the generated files (PNG, PDF, and optionally CSV summary).
    """
    if canonical_amino_acids is None:
        canonical_amino_acids = DEFAULT_CANONICAL_AMINO_ACIDS

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = _read_mapping_csv(pathlib.Path(input_csv), delimiter=delimiter)

    # Count how many 3-letter codes map to each 1-letter code
    # We count rows grouped by `one_letter` (not unique 3-letter solely), because
    # the CSV often already represents each 3-letter as one row.
    counts = df.groupby("one_letter", as_index=False)["three_letter"].count()
    counts.rename(columns={"three_letter": "count"}, inplace=True)

    # Sort descending by count
    counts.sort_values(by="count", ascending=False, inplace=True, kind="stable")
    codes = counts["one_letter"].tolist()
    values = counts["count"].tolist()

    # Color canonical vs non-canonical
    # Viridis-based colors (bright = canonical, dark = non-canonical)
    canonical_color = CANONICAL_AMINO_COLORS["Canonical"]
    non_canonical_color = CANONICAL_AMINO_COLORS["Non-canonical"]
    colors = [canonical_color if c in canonical_amino_acids else non_canonical_color for c in codes]

    # Plot
    plt.figure(figsize=(14, 8))
    bars = plt.bar(codes, values, color=colors, edgecolor="black", alpha=0.85, linewidth=0.6)

    # Label values atop bars
    for bar, count in zip(bars, values):
        # Offset text a bit above the bar height
        y = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2.0,
            y + max(0.02 * max(values), 0.5),
            f"{int(count)}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )

    # Titles and labels
    plt.xlabel("Single-letter amino acid code (sorted by frequency)", fontsize=12, fontweight="bold")
    plt.ylabel("Count of three-letter identifiers", fontsize=12, fontweight="bold")
    plt.title(
        "Distribution of Three-Letter Amino Acid Codes\nMapped to Single-Letter Codes",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    # Legend
    legend_elements = [
        Patch(facecolor=canonical_color, edgecolor="black", label="Canonical (20 standard AA)"),
        Patch(facecolor=non_canonical_color, edgecolor="black", label="Non-canonical"),
    ]
    plt.legend(handles=legend_elements, loc="upper right", fontsize=10)

    # Grid, ticks, layout
    plt.grid(axis="y", alpha=0.3, linestyle="--")
    plt.xticks(rotation=0, ha="center")
    plt.tight_layout()

    # Save outputs
    figures_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    png_path = figures_dir / f"{output_filename_base}.png"
    pdf_path = figures_dir / f"{output_filename_base}.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    logger.info(f"Wrote plot: {png_path}")
    logger.info(f"Wrote plot: {pdf_path}")

    written: List[pathlib.Path] = [png_path, pdf_path]

    # Optional CSV summary
    if write_csv_summary:
        # Add type column (canonical / non-canonical)
        summary = counts.copy()
        summary.insert(
            loc=2,
            column="type",
            value=["Canonical" if c in canonical_amino_acids else "Non-canonical" for c in summary["one_letter"]],
        )
        tables_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        csv_summary_path = tables_dir / f"{output_filename_base}_data.csv"
        summary.to_csv(csv_summary_path, index=False)
        logger.info(f"Wrote summary CSV: {csv_summary_path}")
        written.append(csv_summary_path)

    plt.close()

    # Simple console stats
    canonical_codes = [c for c in codes if c in canonical_amino_acids]
    non_canonical_codes = [c for c in codes if c not in canonical_amino_acids]
    canonical_total = sum(counts[counts["one_letter"].isin(canonical_codes)]["count"])
    non_canonical_total = sum(counts[counts["one_letter"].isin(non_canonical_codes)]["count"])
    if values:
        logger.info(
            "Summary: total three-letter entries=%d, unique single-letter codes=%d, "
            "canonical total=%d, non-canonical total=%d, most frequent=%s (%d), least frequent=%s (%d)",
            sum(values),
            len(values),
            int(canonical_total),
            int(non_canonical_total),
            codes[0],
            int(values[0]),
            codes[-1],
            int(values[-1]),
        )

    return written


# --- LaTeX table generation ------------------------------------------------------------------------

def _escape_latex(text: str) -> str:
    """
    Escape special LaTeX characters in text.
    Order matters: escape backslash first to avoid double-escaping.
    """
    if text is None:
        return ""
    out = str(text)
    out = out.replace("\\", r"\textbackslash{}")
    out = out.replace("&", r"\&")
    out = out.replace("%", r"\%")
    out = out.replace("$", r"\$")
    out = out.replace("#", r"\#")
    out = out.replace("^", r"\textasciicircum{}")
    out = out.replace("_", r"\_")
    out = out.replace("{", r"\{")
    out = out.replace("}", r"\}")
    out = out.replace("~", r"\textasciitilde{}")
    return out


def _format_chemical_name_for_latex(name: str) -> str:
    """
    Apply simple chemistry-aware formatting after escaping:
    - Subscripts for digits following letters: e.g., C6H6 -> C$_{6}$H$_{6}$
    - Superscripts for charges: e.g., 2+ -> $^{2+}$
    """
    s = _escape_latex(name)

    # Subscripts for digits immediately following letters (simple heuristic)
    s = re.sub(r"([A-Za-z])(\d+)", r"\1$_{\2}$", s)

    # Superscripts for charges like 2+ or 3-
    s = re.sub(r"(\d+)([+-])", r"$^{\1\2}$", s)

    return s


def generate_latex_table(
    input_csv: pathlib.Path,
    output_path: pathlib.Path,
    delimiter: str = ";",
    table_caption: str = "Mapping from three-letter amino acid identifiers to one-letter codes",
    table_label: str = "tab:amino_acid_mapping",
) -> pathlib.Path:
    """
    Generate a publication-ready LaTeX longtable from the amino acid mapping CSV.

    The output requires the following LaTeX packages in your document preamble:
      - \\usepackage{booktabs}
      - \\usepackage{longtable}
      - \\usepackage{siunitx}

    Args:
        input_csv: Path to input CSV.
        output_path: Path to write the LaTeX file (.tex).
        delimiter: CSV delimiter (default ';').
        table_caption: Caption for the table.
        table_label: LaTeX label for referencing the table.

    Returns:
        The path to the written LaTeX file.
    """
    df = _read_mapping_csv(pathlib.Path(input_csv), delimiter=delimiter)

    # Sort by one_letter then three_letter for readability
    df_sorted = df.sort_values(by=["one_letter", "three_letter"], kind="stable")

    # Build LaTeX content
    lines: List[str] = []
    lines.append("% Requires: \\usepackage{booktabs}, \\usepackage{longtable}, \\usepackage{siunitx}")
    lines.append("\\begin{longtable}{lll}")
    lines.append(f"  \\caption{{{_escape_latex(table_caption)}}} \\label{{{_escape_latex(table_label)}}} \\\\")
    lines.append("  \\toprule")
    lines.append("  \\textbf{Three-letter} & \\textbf{One-letter} & \\textbf{Name} \\\\")
    lines.append("  \\midrule")
    lines.append("  \\endfirsthead")
    lines.append("  \\multicolumn{3}{c}{\\textit{Continued from previous page}} \\\\")
    lines.append("  \\toprule")
    lines.append("  \\textbf{Three-letter} & \\textbf{One-letter} & \\textbf{Name} \\\\")
    lines.append("  \\midrule")
    lines.append("  \\endhead")
    lines.append("  \\midrule")
    lines.append("  \\multicolumn{3}{c}{\\textit{Continued on next page}} \\\\")
    lines.append("  \\endfoot")
    lines.append("  \\bottomrule")
    lines.append("  \\endlastfoot")

    for _, row in df_sorted.iterrows():
        three = _escape_latex(str(row.get("three_letter", "")).strip())
        one = _escape_latex(str(row.get("one_letter", "")).strip())
        name = _format_chemical_name_for_latex(str(row.get("name", "")).strip())
        lines.append(f"  {three} & {one} & {name} \\\\")

    lines.append("\\end{longtable}")
    latex_content = "\n".join(lines) + "\n"

    output_path = pathlib.Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(latex_content, encoding="utf-8")
    logger.info(f"Wrote LaTeX table: {output_path}")

    return output_path


__all__ = [
    "DEFAULT_CANONICAL_AMINO_ACIDS",
    "plot_amino_acid_distribution",
    "generate_latex_table",
]
