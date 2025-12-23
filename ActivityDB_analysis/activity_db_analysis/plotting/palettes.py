"""
Centralized Viridis-based palette utilities.

This module standardizes all plotting color usage across the project on a
6‑color discrete sampling of the perceptually uniform *viridis* colormap,
plus helper functions for deriving smaller (or larger) discrete sets,
continuous colormaps, and canonical semantic mappings used by multiple
plotting modules.

Color Philosophy
----------------
The user requirement: "brighter means warmer or more contrast."
Accordingly we treat the *brightest / yellow* end of viridis as the
high / hot / maximum side and the *dark purple* end as the low / cold /
minimum side.

Provided / Locked 6-Color Set (bright → dark):
    0: #fde725  (bright yellow)
    1: #7ad151  (yellow‑green)
    2: #22a884  (teal‑green)
    3: #2a788e  (blue‑teal)  ← use for single‑color (monochrome) plots
    4: #414487  (indigo)
    5: #440154  (dark purple)

For ordered / quantitative categories, map increasing magnitude or
quality to earlier indices (lower index = brighter). For reversed
interpretations you can simply reverse the list.

Utilities exported here avoid ad‑hoc color declarations scattered
throughout the codebase and make future palette swaps trivial.

Usage Examples
--------------
    from activity_db_analysis.plotting.palettes import (
        VIRIDIS_6_BRIGHT_TO_DARK, mono_color, get_viridis_discrete
    )

    # 3 ordered categories (high, mid, low):
    colors = get_viridis_discrete(3)              # ['#fde725', '#22a884', '#440154']

    # 3 ordered categories (low, mid, high):
    colors_rev = list(reversed(get_viridis_discrete(3)))

    # Single-color plot (bars / outline emphasis):
    base_color = mono_color()

    # Semantic mappings (e.g. hit coverage buckets):
    from activity_db_analysis.plotting import palettes
    bucket_color = palettes.HIT_COVERAGE_BUCKET_COLORS["40-<70"]

"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence

import matplotlib.cm as cm
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Fixed Discrete Anchor Colors
# ---------------------------------------------------------------------------

# Bright (high) → dark (low) ordering (locked to user-specified hex values)
VIRIDIS_6_BRIGHT_TO_DARK: List[str] = [
    "#fde725",
    "#7ad151",
    "#22a884",
    "#2a788e",
    "#414487",
    "#440154",
]

# Convenience reversed ordering (dark → bright)
VIRIDIS_6_DARK_TO_BRIGHT: List[str] = list(reversed(VIRIDIS_6_BRIGHT_TO_DARK))

# Monochrome default (specified by user: use #2a788e)
MONOCHROME_COLOR: str = "#2a788e"


# ---------------------------------------------------------------------------
# Core Functions
# ---------------------------------------------------------------------------

def mono_color() -> str:
    """
    Return the standard single (monochrome) plot color.

    Chosen mid-tone (#2a788e) provides good contrast against both light
    and dark backgrounds and sits centrally in perceptual luminance.
    """
    return MONOCHROME_COLOR


def get_viridis_discrete(n: int, bright_high: bool = True) -> List[str]:
    """
    Return a list of `n` discrete viridis-based colors.

    Behavior:
      - If n == 1: returns [MONOCHROME_COLOR]
      - If 2 ≤ n ≤ 6: returns an optimally spaced subset of the locked 6-color list.
      - If n  > 6: samples evenly from the continuous matplotlib viridis colormap.

    Args
    ----
    n : int
        Number of discrete colors requested (>=1).
    bright_high : bool
        If True (default) preserve bright → dark ordering with the first
        element being the brightest (interpreted as "highest" category).
        If False, the list is returned dark → bright.

    Returns
    -------
    list[str]
        List of hex color strings.

    Notes
    -----
    For n < 6 we pick indices that maximize perceptual separation while
    respecting ordering. Example (bright_high=True):
        n=2 => [0, 5]
        n=3 => [0, 2, 5]
        n=4 => [0, 2, 3, 5]  (keeps both mid tones separated)
        n=5 => [0, 1, 2, 4, 5]
        n=6 => full set
    """
    if n < 1:
        raise ValueError("n must be >= 1")

    base = VIRIDIS_6_BRIGHT_TO_DARK if bright_high else VIRIDIS_6_DARK_TO_BRIGHT

    if n == 1:
        return [MONOCHROME_COLOR]

    if n <= 6:
        index_map = {
            2: [0, 5],
            3: [0, 2, 5],
            4: [0, 2, 3, 5],
            5: [0, 1, 2, 4, 5],
            6: list(range(6)),
        }
        return [base[i] for i in index_map[n]]

    # n > 6 → sample the continuous viridis colormap evenly
    continuous = cm.get_cmap("viridis", n)
    return [mcolors.to_hex(continuous(i)) for i in range(n)]


def get_gradient_colormap(name: str = "viridis_custom") -> mcolors.LinearSegmentedColormap:
    """
    Build a LinearSegmentedColormap spanning exactly the 6 locked anchor colors
    (bright → dark). Useful when you want to force interpolation only through
    the provided hex values rather than matplotlib's default continuous viridis.

    Parameters
    ----------
    name : str
        Name for the resulting colormap (namespace only).

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    return mcolors.LinearSegmentedColormap.from_list(name, VIRIDIS_6_BRIGHT_TO_DARK)


# ---------------------------------------------------------------------------
# Semantic / Domain-Specific Color Mappings
# ---------------------------------------------------------------------------

# Amino acid canonical vs non-canonical
CANONICAL_AMINO_COLORS = {
    "Canonical": VIRIDIS_6_BRIGHT_TO_DARK[0],       # bright / emphasis
    "Non-canonical": VIRIDIS_6_BRIGHT_TO_DARK[-1],  # dark
}

# Sequence type (3 categories) – spaced bright, mid, dark
SEQUENCE_TYPE_COLORS = {
    "ATOMSEQ": VIRIDIS_6_BRIGHT_TO_DARK[0],
    "SEQRES":  VIRIDIS_6_BRIGHT_TO_DARK[2],
    "SEQADV":  VIRIDIS_6_BRIGHT_TO_DARK[-1],
}

# Expansion metric ordering (quadruplet highest salience → pdblid)
EXPANSION_METRIC_COLORS = {
    "quadruplet_expansion_factor": VIRIDIS_6_BRIGHT_TO_DARK[0],
    "triplet_expansion_factor":    VIRIDIS_6_BRIGHT_TO_DARK[1],
    "pair_expansion_factor":       VIRIDIS_6_BRIGHT_TO_DARK[2],
    "pdblid_expansion_factor":     VIRIDIS_6_BRIGHT_TO_DARK[4],
}

# Hit coverage / quality buckets (6 including NULL)
# Brightest = highest quality; darkest = NULL/unknown.
HIT_COVERAGE_BUCKET_COLORS = {
    ">=80 (Strict)": VIRIDIS_6_BRIGHT_TO_DARK[0],
    "70-<80":         VIRIDIS_6_BRIGHT_TO_DARK[1],
    "40-<70":         VIRIDIS_6_BRIGHT_TO_DARK[2],
    "10-<40":         VIRIDIS_6_BRIGHT_TO_DARK[3],
    "<10":            VIRIDIS_6_BRIGHT_TO_DARK[4],
    "NULL":           VIRIDIS_6_BRIGHT_TO_DARK[5],
}

# Relaxed-only composition quality (no NULL distinguished separately here)
RELAXED_ONLY_QUALITY_COLORS = {
    ">=70":   VIRIDIS_6_BRIGHT_TO_DARK[0],
    "40-<70": VIRIDIS_6_BRIGHT_TO_DARK[2],
    "10-<40": VIRIDIS_6_BRIGHT_TO_DARK[3],
    "<10":    VIRIDIS_6_BRIGHT_TO_DARK[4],
    "NULL":   VIRIDIS_6_BRIGHT_TO_DARK[5],
}

# Representation (smmclid 1..6) – assuming 1 is "better" (more specific) -> brightest.
REPRESENTATION_COLORS_BY_SMMCLID = {
    1: VIRIDIS_6_BRIGHT_TO_DARK[0],
    2: VIRIDIS_6_BRIGHT_TO_DARK[1],
    3: VIRIDIS_6_BRIGHT_TO_DARK[2],
    4: VIRIDIS_6_BRIGHT_TO_DARK[3],
    5: VIRIDIS_6_BRIGHT_TO_DARK[4],
    6: VIRIDIS_6_BRIGHT_TO_DARK[5],
}

# Same mapping but keyed by the verbose labels (after translation)
def build_representation_label_colors(label_sequence: Sequence[str]) -> dict[str, str]:
    """
    Construct a mapping from full representation labels (already ordered
    by smmclid ascending) to viridis colors. This helper lets callers pass
    in whatever label list they actually use (after translation).

    Example
    -------
        labels = [
            "InChI-Str Connection Layer",
            "InChI-Str Connection & Hydrogen Layer",
            "SMI w/o Stereo",
            "SMI w/ Stereo",
            "InChI-Str",
            "InChI-Key",
        ]
        rep_colors = build_representation_label_colors(labels)
    """
    colors = get_viridis_discrete(len(label_sequence))
    return {label: colors[i] for i, label in enumerate(label_sequence)}


# ---------------------------------------------------------------------------
# Convenience Data Structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PaletteSet:
    """
    Bundle of semantic palette dictionaries for convenient star import usage
    or IDE exploration.
    """
    canonical_amino: dict
    sequence_type: dict
    expansion_metrics: dict
    hit_coverage_buckets: dict
    relaxed_only_quality: dict
    representation_by_smmclid: dict


PALETTES = PaletteSet(
    canonical_amino=CANONICAL_AMINO_COLORS,
    sequence_type=SEQUENCE_TYPE_COLORS,
    expansion_metrics=EXPANSION_METRIC_COLORS,
    hit_coverage_buckets=HIT_COVERAGE_BUCKET_COLORS,
    relaxed_only_quality=RELAXED_ONLY_QUALITY_COLORS,
    representation_by_smmclid=REPRESENTATION_COLORS_BY_SMMCLID,
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    # Fixed lists
    "VIRIDIS_6_BRIGHT_TO_DARK",
    "VIRIDIS_6_DARK_TO_BRIGHT",
    "MONOCHROME_COLOR",
    # Functions
    "mono_color",
    "get_viridis_discrete",
    "get_gradient_colormap",
    "build_representation_label_colors",
    # Semantic mappings
    "CANONICAL_AMINO_COLORS",
    "SEQUENCE_TYPE_COLORS",
    "EXPANSION_METRIC_COLORS",
    "HIT_COVERAGE_BUCKET_COLORS",
    "RELAXED_ONLY_QUALITY_COLORS",
    "REPRESENTATION_COLORS_BY_SMMCLID",
    # Aggregate container
    "PALETTES",
]
