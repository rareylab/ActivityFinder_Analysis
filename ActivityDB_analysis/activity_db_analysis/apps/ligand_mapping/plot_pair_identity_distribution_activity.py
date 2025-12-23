#!/usr/bin/env python3
"""
plot_pair_identity_distribution_activity
=======================================

Activity-qualified (ANY ChEMBL activity row), coverage-agnostic ligand pair
identity distribution plotting utilities.

Input CSV
---------
Expected columns produced by PairIdentityDistributionActivity app:
    identity_interval     (categorical, one of: [100], [95,100), [90,95), [85,90), [80,85))
    smmclid               (int; raw structural matching id: 1,2,3,4,6 with 5 retired)
    confidence_level      (int; compressed contiguous scale 1..5 where smmclid 6 → confidence_level 5)
    pair_count            (int)
    interval_total        (int)  -- total pairs within that identity interval (all confidence levels)
    pct_within_interval   (float) -- pair_count / interval_total * 100
    pct_global            (float) -- pair_count / global_total * 100

Plots
-----
1. Stacked bar (absolute counts) per identity interval, stacked by smmclid.
2. Stacked bar (percent of interval) per identity interval, stacked by smmclid.
   (Optional annotations: interval_total on top of each bar.)

Color Scheme
------------
Uses 5-color viridis subset (light→dark) mapped by confidence_level (1..5):
    confidence_level 1: "#fde725"
    confidence_level 2: "#5ec962"
    confidence_level 3: "#21918c"
    confidence_level 4: "#3b528b"
    confidence_level 5: "#440154"

CLI
---
Examples:
    python -m activity_db_analysis.apps.ligand_mapping.plot_pair_identity_distribution_activity \
        --input pair_identity_distribution_activity.csv \
        --output-dir ./figures \
        --counts --percents --annotate

Programmatic:
    from activity_db_analysis.apps.ligand_mapping.plot_pair_identity_distribution_activity import (
        create_pair_identity_distribution_activity_plots
    )
    create_pair_identity_distribution_activity_plots(Path("pair_identity_distribution_activity.csv"), Path("./figs"))

"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger

# Optional seaborn styling if available (graceful fallback)
try:
    import seaborn as sns  # type: ignore

    _HAS_SEABORN = True
except ImportError:  # pragma: no cover
    _HAS_SEABORN = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INTERVAL_ORDER = ["[100]", "[95,100)", "[90,95)", "[85,90)", "[80,85)"]

SMMCLID_COLORS: Dict[int, str] = {
    # Keys now represent confidence_level (compressed 1..5)
    1: "#fde725",
    2: "#5ec962",
    3: "#21918c",
    4: "#3b528b",
    5: "#440154",
}

DEFAULT_DPI = 250

# Central text scaling factor (130% of original sizes)
TEXT_SCALE = 1.3


# ---------------------------------------------------------------------------
# Formatting Helpers
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Data Loading & Preparation
# ---------------------------------------------------------------------------


def _read_input(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    if df.empty:
        return df
    expected_cols = {
        "identity_interval",
        "smmclid",
        "confidence_level",
        "pair_count",
        "interval_total",
        "pct_within_interval",
        "pct_global",
    }
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected columns in input CSV: {sorted(missing)}")
    df["pair_count"] = (
        pd.to_numeric(df["pair_count"], errors="coerce").fillna(0).astype(int)
    )
    df["interval_total"] = (
        pd.to_numeric(df["interval_total"], errors="coerce").fillna(0).astype(int)
    )
    df["pct_within_interval"] = pd.to_numeric(
        df["pct_within_interval"], errors="coerce"
    )
    df["pct_global"] = pd.to_numeric(df["pct_global"], errors="coerce")
    df["smmclid"] = pd.to_numeric(df["smmclid"], errors="coerce").astype(int)
    # confidence_level should already be present (written by data extraction step);
    # if absent (older CSV) synthesize it by compressing smmclid (6→5).
    if "confidence_level" not in df.columns:
        df["confidence_level"] = (
            df["smmclid"].map({1: 1, 2: 2, 3: 3, 4: 4, 6: 5}).astype(int)
        )
    df["identity_interval"] = pd.Categorical(
        df["identity_interval"], categories=INTERVAL_ORDER, ordered=True
    )
    df = df.sort_values(["identity_interval", "smmclid"])
    return df


def _pivot_counts(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    pivot = (
        df.pivot_table(
            index="identity_interval",
            columns="confidence_level",
            values="pair_count",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(INTERVAL_ORDER)
        .fillna(0)
    )
    desired_cols = [1, 2, 3, 4, 5, 6]
    for col in desired_cols:
        if col not in pivot.columns:
            pivot[col] = 0
    pivot = pivot[desired_cols]

    recomputed_totals = pivot.sum(axis=1)
    provided_totals = (
        df.groupby("identity_interval", observed=True)["interval_total"]
        .max()
        .reindex(INTERVAL_ORDER)
    )
    if not provided_totals.equals(recomputed_totals):
        diffs = (recomputed_totals - provided_totals).fillna(0)
        diff_str = ", ".join(
            f"{idx}: provided={int(provided_totals.get(idx, 0))} recomputed={int(recomputed_totals.get(idx, 0))}"
            for idx in INTERVAL_ORDER
            if not pd.isna(diffs.get(idx, 0)) and diffs.get(idx, 0) != 0
        )
        if diff_str:
            logger.warning(
                "Interval total mismatch detected; using recomputed totals. Differences: {}",
                diff_str,
            )
    interval_totals = recomputed_totals
    return pivot, interval_totals


def _derive_percents_from_counts(pivot_counts: pd.DataFrame) -> pd.DataFrame:
    """
    Derive normalized (fractional) composition within each identity interval.

    Returns values in range 0..1 (NOT percentages). These feed 'normalized'
    panels which display full-height stacked bars analogous to target mapping plots.
    """
    totals = pivot_counts.sum(axis=1).replace(0, np.nan)
    pivot_norm = pivot_counts.div(totals, axis=0).fillna(0)  # fractions 0..1
    return pivot_norm


# ---------------------------------------------------------------------------
# Plotting Helpers
# ---------------------------------------------------------------------------


def _apply_theme():
    if _HAS_SEABORN:
        sns.set_theme(style="whitegrid", context="notebook")
    else:
        plt.style.use("ggplot")


def _get_plot_order(cols) -> List[int]:
    # Return levels in descending order so level 6 is stacked at the bottom and level 1 at the top
    return sorted([int(c) for c in cols], reverse=True)


def _annotate_totals(ax, x_positions, totals: pd.Series):
    for x, (idx, val) in zip(x_positions, totals.items()):
        if np.isnan(val):
            continue
        ax.text(
            x,
            val * 1.01 if val != 0 else 0.05,
            _format_number_latex_style(val),
            ha="center",
            va="bottom",
            fontsize=9 * TEXT_SCALE,
            fontweight="bold",
            color="#3b528b",
        )


def _annotate_pct_totals(ax, x_positions):
    for x in x_positions:
        ax.text(
            x,
            101.0,
            "100%",
            ha="center",
            va="bottom",
            fontsize=9 * TEXT_SCALE,
            fontweight="bold",
            color="#3b528b",
        )


def _stacked_bar(
    pivot: pd.DataFrame,
    output_path: Path,
    interval_totals: pd.Series,
    title: str,
    ylabel: str,
    annotate: bool,
    pct_mode: bool = False,
    dpi: int = DEFAULT_DPI,
    save_pdf: bool = False,
):
    _apply_theme()
    fig, ax = plt.subplots(figsize=(11, 5))  # Keep figure size; only scale text
    plt.rcParams.update(
        {
            "axes.titlesize": 18 * TEXT_SCALE,
            "axes.labelsize": 16 * TEXT_SCALE,
            "xtick.labelsize": 14 * TEXT_SCALE,
            "ytick.labelsize": 14 * TEXT_SCALE,
            "legend.fontsize": 14 * TEXT_SCALE,
        }
    )
    x = np.arange(len(pivot.index))
    bottoms = np.zeros(len(pivot.index), dtype=float)

    # Filter to active structural confidence levels (level 5 retired)
    allowed_levels = {1, 2, 3, 4, 5}
    smmclid_levels = [
        lvl for lvl in _get_plot_order(pivot.columns) if lvl in allowed_levels
    ]
    for smmclid in smmclid_levels:
        values = pivot[smmclid].to_numpy()
        ax.bar(
            x,
            values,
            bottom=bottoms,
            color=SMMCLID_COLORS.get(smmclid, "#3b528b"),
            edgecolor="black",
            linewidth=0.5,
            label=f"Level {smmclid}",
        )
        bottoms += values

    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, fontsize=14 * TEXT_SCALE)
    ax.set_xlabel("Maximum Sequence Identity Interval", fontsize=16 * TEXT_SCALE)
    ax.set_ylabel(ylabel, fontsize=16 * TEXT_SCALE)
    ax.set_title(title, fontsize=18 * TEXT_SCALE)
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        frameon=False,
        fontsize=9 * TEXT_SCALE,
        title="smmclid",
        ncol=1,
    )
    if pct_mode:
        # Normalized (fraction) panel
        ax.set_ylim(0, 1.0)
        # No top annotations for normalized panels (bars fill to 1)
    else:
        y_max = bottoms.max() * 1.12 if bottoms.size else 1
        ax.set_ylim(0, y_max if y_max > 0 else 1)
        if annotate:
            _annotate_totals(ax, x, interval_totals)
    fig.subplots_adjust(top=0.85)
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _combined_panel(
    disjoint_counts: pd.DataFrame,
    disjoint_norm: pd.DataFrame,
    disjoint_interval_totals: pd.Series,
    cumulative_counts: pd.DataFrame,
    cumulative_norm: pd.DataFrame,
    cumulative_interval_totals: pd.Series,
    output_path: Path,
    dpi: int = DEFAULT_DPI,
    save_pdf: bool = False,
):
    _apply_theme()
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(16, 11),
        constrained_layout=False,
        sharex=False,
    )
    ax_d_counts, ax_d_pct, ax_c_counts, ax_c_pct = (
        axes[0, 0],
        axes[0, 1],
        axes[1, 0],
        axes[1, 1],
    )

    def _render_panel(
        ax: plt.Axes,
        pivot: pd.DataFrame,
        interval_totals: pd.Series | None,
        title: str,
        pct_mode: bool,
        ylabel: str | None,
    ):
        x = np.arange(len(pivot.index))
        bottoms = np.zeros(len(pivot.index), dtype=float)
        allowed_levels = {1, 2, 3, 4, 5}
        smmclid_levels = [
            lvl for lvl in _get_plot_order(pivot.columns) if lvl in allowed_levels
        ]
        for smmclid in smmclid_levels:
            vals = pivot[smmclid].to_numpy()
            ax.bar(
                x,
                vals,
                bottom=bottoms,
                color=SMMCLID_COLORS.get(smmclid, "#3b528b"),
                edgecolor="black",
                linewidth=0.6,
            )
            bottoms += vals
        ax.set_xticks(x)
        ax.set_xticklabels(pivot.index, fontsize=13 * TEXT_SCALE)
        ax.set_title(title, fontsize=15 * TEXT_SCALE)
        ax.grid(axis="y", linestyle="--", alpha=0.3)
        if pct_mode:
            # Normalized (fraction) panel
            ax.set_ylim(0, 1.0)
        else:
            y_max = bottoms.max() * 1.12 if bottoms.size else 1
            ax.set_ylim(0, y_max if y_max > 0 else 1)
            if interval_totals is not None:
                for xi, (idx, val) in zip(x, interval_totals.items()):
                    if np.isnan(val):
                        continue
                    ax.text(
                        xi,
                        val * 1.01 if val != 0 else 0.05,
                        _format_number_latex_style(val),
                        ha="center",
                        va="bottom",
                        fontsize=11 * TEXT_SCALE,
                        fontweight="bold",
                        color="#3b528b",
                    )
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=14 * TEXT_SCALE)

    _render_panel(
        ax=ax_d_counts,
        pivot=disjoint_counts,
        interval_totals=disjoint_interval_totals,
        title="Distinct Pair Counts",
        pct_mode=False,
        ylabel="Distinct Pair Counts",
    )
    _render_panel(
        ax=ax_d_pct,
        pivot=disjoint_norm,
        interval_totals=None,
        title="Distinct — normalized",
        pct_mode=True,
        ylabel="Proportion within interval",
    )
    _render_panel(
        ax=ax_c_counts,
        pivot=cumulative_counts,
        interval_totals=cumulative_interval_totals,
        title="Interval Pair Counts",
        pct_mode=False,
        ylabel="Pair Counts",
    )
    _render_panel(
        ax=ax_c_pct,
        pivot=cumulative_norm,
        interval_totals=None,
        title="Normalized",
        pct_mode=True,
        ylabel="Proportion within interval",
    )

    allowed_levels = {1, 2, 3, 4, 5}
    smmclid_levels = [5, 4, 3, 2, 1]
    # smmclid_levels = [lvl for lvl in _get_plot_order(disjoint_counts.columns) if lvl in allowed_levels]
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markersize=13 * TEXT_SCALE,
            markerfacecolor=SMMCLID_COLORS.get(s, "#3b528b"),
            markeredgecolor="black",
            label=f"Level {s}",
        )
        for s in smmclid_levels
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=min(len(legend_handles), 6),
        bbox_to_anchor=(0.5, 0.972),
        frameon=True,
        title="Small Molecule Confidence Level",
        columnspacing=1.4,
        handletextpad=0.7,
        fontsize=13 * TEXT_SCALE,
        title_fontsize=14 * TEXT_SCALE,
    )
    fig.suptitle(
        "Activity-Qualified Pair Identity Distribution (Counts & Normalized)",
        fontsize=20 * TEXT_SCALE,
        y=0.998,
    )
    fig.subplots_adjust(
        top=0.84,
        left=0.07,
        right=0.98,
        hspace=0.28,
        wspace=0.18,
        bottom=0.09,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return output_path


def create_pair_identity_distribution_activity_plots(
    input_csv: Path,
    output_dir: Path,
    output_prefix: str = "pair_identity_distribution_activity",
    counts: bool = True,
    percents: bool = True,  # 'percents' flag retained for CLI compatibility; now produces normalized (0..1) panels
    annotate: bool = True,
    dpi: int = DEFAULT_DPI,
    save_png: bool = True,
    save_pdf: bool = False,
) -> List[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    df = _read_input(input_csv)
    generated: List[Path] = []
    if df.empty:
        if counts or percents:
            fig, ax = plt.subplots(figsize=(4, 2))
            ax.text(0.5, 0.5, "No Data", ha="center", va="center")
            ax.axis("off")
            if counts and save_png:
                p = output_dir / f"{output_prefix}_counts.png"
                fig.savefig(p, dpi=dpi)
                generated.append(p)
            if percents and save_png:
                p = output_dir / f"{output_prefix}_percents.png"
                fig.savefig(p, dpi=dpi)
                generated.append(p)
            plt.close(fig)
        return generated
    pivot_counts, interval_totals = _pivot_counts(df)
    pivot_norm = _derive_percents_from_counts(pivot_counts)  # fractions 0..1
    if counts:
        if save_png:
            path_counts = output_dir / f"{output_prefix}_counts.png"
            _stacked_bar(
                pivot_counts,
                path_counts,
                interval_totals=interval_totals,
                title="Activity-Qualified Pair Identity Distribution (Counts)",
                ylabel="Distinct Pairs",
                annotate=annotate,
                pct_mode=False,
                dpi=dpi,
                save_pdf=save_pdf,
            )
            generated.append(path_counts)
        elif save_pdf:
            path_counts_pdf = output_dir / f"{output_prefix}_counts.pdf"
            _stacked_bar(
                pivot_counts,
                path_counts_pdf,
                interval_totals=interval_totals,
                title="Activity-Qualified Pair Identity Distribution (Counts)",
                ylabel="Distinct Pairs",
                annotate=annotate,
                pct_mode=False,
                dpi=dpi,
                save_pdf=True,
            )
            generated.append(path_counts_pdf)
    if percents:
        if save_png:
            path_norm = output_dir / f"{output_prefix}_percents.png"
            _stacked_bar(
                pivot_norm,
                path_norm,
                interval_totals=interval_totals,
                title="Activity-Qualified Pair Identity Distribution (Normalized)",
                ylabel="Proportion within interval",
                annotate=False,  # no top annotations for normalized panel
                pct_mode=True,
                dpi=dpi,
                save_pdf=save_pdf,
            )
            generated.append(path_norm)
        elif save_pdf:
            path_norm_pdf = output_dir / f"{output_prefix}_percents.pdf"
            _stacked_bar(
                pivot_norm,
                path_norm_pdf,
                interval_totals=interval_totals,
                title="Activity-Qualified Pair Identity Distribution (Normalized)",
                ylabel="Proportion within interval",
                annotate=False,
                pct_mode=True,
                dpi=dpi,
                save_pdf=True,
            )
            generated.append(path_norm_pdf)
    return generated


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot activity-qualified, coverage-agnostic ligand pair identity distribution."
    )
    p.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input CSV produced by PairIdentityDistributionActivity app.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for figure(s).",
    )
    p.add_argument(
        "--output-prefix",
        type=str,
        default="pair_identity_distribution_activity",
        help="Base filename prefix (default: pair_identity_distribution_activity).",
    )
    p.add_argument(
        "--no-counts",
        action="store_true",
        help="Do not generate counts stacked bar plot.",
    )
    p.add_argument(
        "--no-percents",
        action="store_true",
        help="Do not generate percent stacked bar plot.",
    )
    p.add_argument(
        "--no-annotate",
        action="store_true",
        help="Disable total/100% annotations.",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"DPI for raster outputs (default: {DEFAULT_DPI}).",
    )
    p.add_argument(
        "--no-png",
        action="store_true",
        help="Skip PNG output (only useful if --pdf used).",
    )
    p.add_argument(
        "--pdf",
        action="store_true",
        help="Also (or only, if --no-png) write PDF outputs.",
    )
    p.add_argument(
        "--combined",
        action="store_true",
        help="Generate a combined 2x2 panel comparing disjoint vs cumulative distributions (requires both CSVs).",
    )
    p.add_argument(
        "--cumulative-input",
        type=Path,
        required=False,
        help="Path to cumulative CSV (pair_identity_distribution_activity_cumulative.csv) when using --combined.",
    )
    return p


def main():
    args = _build_arg_parser().parse_args()
    create_pair_identity_distribution_activity_plots(
        input_csv=args.input,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        counts=not args.no_counts,
        percents=not args.no_percents,
        annotate=not args.no_annotate,
        dpi=args.dpi,
        save_png=not args.no_png,
        save_pdf=args.pdf,
    )
    if args.combined:
        if not args.cumulative_input:
            raise SystemExit("--combined requires --cumulative-input path.")
        cumulative_df = _read_input(args.cumulative_input)
        disjoint_df = _read_input(args.input)
        disjoint_counts_pivot, disjoint_interval_totals = _pivot_counts(disjoint_df)
        disjoint_norm_pivot = _derive_percents_from_counts(disjoint_counts_pivot)
        cumulative_counts_pivot, cumulative_interval_totals = _pivot_counts(
            cumulative_df
        )
        cumulative_norm_pivot = _derive_percents_from_counts(cumulative_counts_pivot)
        combined_path = args.output_dir / f"{args.output_prefix}_combined_panel.png"
        _combined_panel(
            disjoint_counts=disjoint_counts_pivot,
            disjoint_norm=disjoint_norm_pivot,
            disjoint_interval_totals=disjoint_interval_totals,
            cumulative_counts=cumulative_counts_pivot,
            cumulative_norm=cumulative_norm_pivot,
            cumulative_interval_totals=cumulative_interval_totals,
            output_path=combined_path,
            dpi=args.dpi,
            save_pdf=args.pdf,
        )


if __name__ == "__main__":
    main()
