#!/usr/bin/env python3
"""
Create a 2x2 stacked bar figure for target mappings (interval-binned CSV).

Layout:
- Rows: Component, Variant
- Columns: left = absolute counts, right = normalized (0..1)
- X-axis categories: PDB sequence types (ATOM, SEQRES, SEQADV)
- Stacks: percent identity intervals (five bins)

Input CSV
- Must be the consolidated interval-binned report written by:
  activity_db_analysis.apps.sequence_mapping.sequence_mapping.generate_target_matching_reports
  -> coverage-mode–specific filenames such as:
     target_matches_by_sequences_sequence_identity_ranges_full.csv
     target_matches_by_sequences_sequence_identity_ranges_no_hit.csv
     target_matches_by_sequences_sequence_identity_ranges_identity_only.csv
     (Underscore naming convention enforced: full | no_hit | identity_only)
- Expected columns (examples):
  COMPONENT_ATOMSEQ_RAW, COMPONENT_SEQRES_RAW, COMPONENT_SEQADV_RAW,
  COMPONENT_ATOMSEQ_DENSITY, COMPONENT_SEQRES_DENSITY, COMPONENT_SEQADV_DENSITY,
  COMPONENT_ATOMSEQ_OCC, COMPONENT_SEQRES_OCC, COMPONENT_SEQADV_OCC,
  COMPONENT_ATOMSEQ_COND_OCC, COMPONENT_SEQRES_COND_OCC, COMPONENT_SEQADV_COND_OCC,
  VARIANT_ATOMSEQ_RAW,   VARIANT_SEQRES_RAW,   VARIANT_SEQADV_RAW,
  VARIANT_ATOMSEQ_DENSITY, VARIANT_SEQRES_DENSITY, VARIANT_SEQADV_DENSITY,
  VARIANT_ATOMSEQ_OCC, VARIANT_SEQRES_OCC, VARIANT_SEQADV_OCC,
  VARIANT_ATOMSEQ_COND_OCC, VARIANT_SEQRES_COND_OCC, VARIANT_SEQADV_COND_OCC
- Row index contains the interval labels:
  [100], [95,100), [90,95), [85,90), [80,85)
  plus optional aggregate rows:
    TOTAL_SEQUENCES, MAPPED_SEQUENCES (these are removed before plotting;
    TOTAL_SEQUENCES values are only used to annotate subtitles under x‑ticks).

Primary usage now occurs via the Typer CLI sub-app:
    uv run actfind sequence-mapper plot --input <interval_csv> --output figure.png
or batch:
    uv run actfind sequence-mapper plot-batch --input-dir reports --output-dir figures

Notes / Important clarification:
- This module ONLY CONSUMES *_RAW columns. “Normalized” panels are derived by row-wise normalization
  within each PDB sequence type.
- *_DENSITY, *_OCC, *_COND_OCC columns are ignored for these plots.
- File coverage mode fragments use underscores; human-readable titles may show hyphenated forms.
"""

from __future__ import annotations

import pathlib
from typing import Dict, List, Sequence, Tuple, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.axes import Axes

# Global font scaling (130%)
FONT_SCALE = 1.30


def _fs(size: float) -> float:
    """Scale a base font size by the global FONT_SCALE."""
    return size * FONT_SCALE


plt.rcParams.update({"font.size": _fs(10)})


def _format_number_latex_style(value: int | float) -> str:
    """
    Format a number with thin space separators (LaTeX style) instead of commas.
    Uses Unicode thin space (U+2009) to separate thousands.

    Examples:
        1234 -> "1 234"
        1234567 -> "1 234 567"
    """
    if np.isnan(value):
        return ""
    int_val = int(value)
    # Format with comma separator first, then replace with thin space
    formatted = f"{int_val:,}"
    # Replace comma with Unicode thin space (U+2009)
    return formatted.replace(",", "\u2009")


# Identity interval order and colors (viridis scheme - dark to bright; inverted from original)
INTERVAL_ORDER: List[str] = ["[100]", "[95,100)", "[90,95)", "[85,90)", "[80,85)"]
INTERVAL_COLORS: Dict[str, str] = {
    "[100]": "#440154",
    "[95,100)": "#2a788e",
    "[90,95)": "#22a884",
    "[85,90)": "#7ad151",
    "[80,85)": "#fde725",
}

# PDB sequence types (column suffixes in CSV) and display labels
SEQ_TYPES: List[str] = ["ATOMSEQ", "SEQRES", "SEQADV"]
SEQ_DISPLAY: Dict[str, str] = {
    "ATOMSEQ": "ATOM",
    "SEQRES": "SEQRES",
    "SEQADV": "SEQADV",
}

# Sides (ChEMBL types)
SIDES: List[str] = ["COMPONENT", "VARIANT"]


def _read_csv(csv_path: pathlib.Path) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Read interval-binned report CSV and return (interval_df, totals_by_seqtype).
    """
    df_full = pd.read_csv(csv_path, index_col=0)
    df_full.index = df_full.index.map(lambda s: str(s).strip())

    totals: Dict[str, int] = {}
    if "TOTAL_SEQUENCES" in df_full.index:
        totals_row = df_full.loc["TOTAL_SEQUENCES"]
        for st in SEQ_TYPES:
            for side in ("COMPONENT", "VARIANT"):
                col = f"{side}_{st}_RAW"
                if col in totals_row and pd.notna(totals_row[col]):
                    totals[st] = int(totals_row[col])
                    break
            totals.setdefault(st, 0)

    for _agg in ("TOTAL_SEQUENCES", "MAPPED_SEQUENCES"):
        if _agg in df_full.index:
            df_full = df_full.drop(index=_agg)
    df = df_full
    present = [i for i in INTERVAL_ORDER if i in df.index]
    if not present:
        raise ValueError("Input CSV does not contain any expected interval labels.")
    df = df.reindex(present).fillna(0.0)
    return df, totals


def _extract_side_matrices(
    df: pd.DataFrame, side: str
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    intervals_present = [i for i in INTERVAL_ORDER if i in df.index]
    if not intervals_present:
        raise ValueError(f"No expected intervals present for side: {side}")

    raw_rows: List[List[float]] = []
    for st in SEQ_TYPES:
        col = f"{side}_{st}_RAW"
        raw_vals = [
            float(df.loc[i, col]) if col in df.columns else 0.0
            for i in intervals_present
        ]
        raw_rows.append(raw_vals)

    raw_matrix = np.asarray(raw_rows, dtype=float)
    row_sums = raw_matrix.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        norm_matrix = np.divide(raw_matrix, row_sums, where=(row_sums != 0.0))
    norm_matrix = np.nan_to_num(norm_matrix, nan=0.0, posinf=0.0, neginf=0.0)
    return raw_matrix, norm_matrix, intervals_present


def _stacked_bars_by_intervals(
    ax: Axes,
    matrix: np.ndarray,
    intervals: Sequence[str],
    title: str,
    ylabel: str,
    normalized: bool,
    ylim: Tuple[float, float] | None = None,
    xtick_subtitles: Sequence[str] | None = None,
) -> None:
    n_seq = matrix.shape[0]
    x = np.arange(n_seq)
    bottoms = np.zeros(n_seq, dtype=float)

    for j in range(matrix.shape[1]):
        heights = matrix[:, j]
        label = intervals[j]
        color = INTERVAL_COLORS.get(label, "#414487")
        ax.bar(
            x,
            heights,
            bottom=bottoms,
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.6,
            alpha=0.95,
        )
        bottoms = bottoms + heights

    if not normalized:
        totals = matrix.sum(axis=1)
        for xi, total in zip(x, totals):
            ax.text(
                xi,
                float(total) * 1.005
                if float(total) > 0
                else 0.02 * (ylim[1] if ylim else 1.0),
                _format_number_latex_style(int(round(total))),
                ha="center",
                va="bottom",
                fontsize=_fs(9),
                fontweight="bold",
            )

    ax.set_title(title, fontsize=_fs(13), pad=8)
    ax.set_xticks(x)
    base_labels = [SEQ_DISPLAY.get(st, st) for st in SEQ_TYPES]
    if xtick_subtitles:
        labels = [f"{lbl}\n{sub}" for lbl, sub in zip(base_labels, xtick_subtitles)]
    else:
        labels = base_labels
    ax.set_xticklabels(labels)
    ax.tick_params(axis="x", which="both", pad=6, labelsize=_fs(10))
    ax.tick_params(axis="y", which="both", labelsize=_fs(10))
    ax.set_ylabel(ylabel, fontsize=_fs(11))
    ax.grid(axis="y", linestyle="--", alpha=0.35)

    if normalized:
        ax.set_ylim(0.0, 1.0)
    elif ylim is not None:
        ax.set_ylim(*ylim)


def plot_target_matching_interval_stacked_bars(
    csv_path: pathlib.Path,
    out_path: pathlib.Path,
    figure_title: str | None = None,
    additional_output_paths: Iterable[pathlib.Path] | None = None,
) -> None:
    """
    Build a 2x2 figure:
      Row 0: Component — [absolute | normalized]
      Row 1: Variant   — [absolute | normalized]
    If additional_output_paths provided, the same in-memory figure is saved
    to each of those paths (avoids full re-rendering).
    """
    df, totals = _read_csv(csv_path)
    comp_raw, comp_norm, comp_intervals = _extract_side_matrices(df, "COMPONENT")
    var_raw, var_norm, var_intervals = _extract_side_matrices(df, "VARIANT")

    comp_totals = comp_raw.sum(axis=1)
    var_totals = var_raw.sum(axis=1)
    global_max = float(max(comp_totals.max(initial=0.0), var_totals.max(initial=0.0)))
    abs_ylim = (0.0, global_max * 1.05 if global_max > 0 else 1.0)
    var_max = float(var_totals.max(initial=0.0))
    var_abs_ylim = (0.0, var_max * 1.12 if var_max > 0 else 1.0)

    ROW_SPACING = 0.3
    LEGEND_Y = 0.915
    SUPTITLE_Y = 0.95

    fig, axes = plt.subplots(
        nrows=2, ncols=2, figsize=(14, 10), constrained_layout=False
    )
    fig.subplots_adjust(
        top=0.80, hspace=ROW_SPACING, wspace=0.20, left=0.08, right=0.98, bottom=0.08
    )

    _stacked_bars_by_intervals(
        ax=axes[0, 0],
        matrix=comp_raw,
        intervals=comp_intervals,
        title="Component — absolute counts",
        ylabel="Absolute count",
        normalized=False,
        ylim=abs_ylim,
        xtick_subtitles=[
            _format_number_latex_style(int(totals.get(st, 0))) for st in SEQ_TYPES
        ],
    )
    _stacked_bars_by_intervals(
        ax=axes[0, 1],
        matrix=comp_norm,
        intervals=comp_intervals,
        title="Component — normalized",
        ylabel="Proportion within sequence type",
        normalized=True,
        xtick_subtitles=[
            _format_number_latex_style(int(totals.get(st, 0))) for st in SEQ_TYPES
        ],
    )
    _stacked_bars_by_intervals(
        ax=axes[1, 0],
        matrix=var_raw,
        intervals=var_intervals,
        title="Variant — absolute counts",
        ylabel="Absolute count",
        normalized=False,
        ylim=var_abs_ylim,
        xtick_subtitles=[
            _format_number_latex_style(int(totals.get(st, 0))) for st in SEQ_TYPES
        ],
    )
    _stacked_bars_by_intervals(
        ax=axes[1, 1],
        matrix=var_norm,
        intervals=var_intervals,
        title="Variant — normalized",
        ylabel="Proportion within sequence type",
        normalized=True,
        xtick_subtitles=[
            _format_number_latex_style(int(totals.get(st, 0))) for st in SEQ_TYPES
        ],
    )

    legend_labels = [
        i for i in INTERVAL_ORDER if (i in comp_intervals) or (i in var_intervals)
    ]
    legend_handles = [
        Patch(facecolor=INTERVAL_COLORS.get(i, "#414487"), edgecolor="black", label=i)
        for i in legend_labels
    ]
    if legend_handles:
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, LEGEND_Y),
            ncol=min(5, len(legend_handles)),
            frameon=True,
            fontsize=_fs(11),
            title="Percent identity interval",
            title_fontsize=_fs(12),
        )

    if figure_title:
        fig.suptitle(figure_title, fontsize=_fs(16), y=SUPTITLE_Y)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    if additional_output_paths:
        for extra in additional_output_paths:
            extra.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(extra, dpi=300, bbox_inches="tight")
    plt.close(fig)
