import pathlib
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger
import sys

from activity_db_analysis.plotting.palettes import mono_color
from activity_db_analysis.apps.variant_mutation_plot.variant_mutation_sql import VARIANT_MUTATION_DATA_EXTRACTION_SQL


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


def generate_latex_table(
    df: pd.DataFrame,
    output_path: pathlib.Path,
    table_caption: str = "Distribution of Mutation Counts in ChEMBL Variant Sequences",
    table_label: str = "tab:mutation_count_distribution",
) -> None:
    """
    Generate a LaTeX table from the variant mutation data using regular table environment and siunitx.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns: num_mutations, num_sequences, percentage
    output_path : pathlib.Path
        Path where the LaTeX table file will be saved
    table_caption : str
        Caption for the LaTeX table
    table_label : str
        Label for the LaTeX table reference

    Note
    ----
    Requires \\usepackage{booktabs} and \\usepackage{siunitx} in LaTeX preamble.
    """
    latex_content = r"% Requires \usepackage{booktabs} and \usepackage{siunitx}" + "\n"
    latex_content += r"\begin{table}[!htbp]" + "\n"
    latex_content += r"    \centering" + "\n"
    latex_content += f"    \\caption{{{table_caption} This table is self-contained: columns summarize variant mutation categories, counts, and related statistics; all terminology corresponds directly to the processed mutation dataset so no external description is required to interpret the numbers.}}\n"
    latex_content += f"    \\label{{{table_label}}}\n"
    latex_content += r"    \begin{tabular}{rrr}" + "\n"
    latex_content += r"        \toprule" + "\n"
    latex_content += (
        r"        \textbf{Number of mutations} & \textbf{Number of sequences} & \textbf{Percentage (\%)} \\"
        + "\n"
    )
    latex_content += r"        \midrule" + "\n"

    # Show all individual values using siunitx formatting
    for _, row in df.iterrows():
        latex_content += f"        {int(row['num_mutations'])} & \\num{{{int(row['num_sequences'])}}} & {row['percentage']:.3f} \\\\\n"

    latex_content += r"        \bottomrule" + "\n"
    latex_content += r"    \end{tabular}" + "\n"
    latex_content += r"\end{table}"

    # Write to file
    with open(output_path, "w") as f:
        f.write(latex_content)

    logger.info(f"LaTeX table saved to {output_path}")


def plot_mutation_distribution(
    data_df: pd.DataFrame,
    output_dir: pathlib.Path,
    output_filename_base: str = "mutation_count_distribution",
) -> None:
    """
    Create bar plots for mutation count distribution.

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame with columns: num_mutations, num_sequences, percentage
    output_dir : pathlib.Path
        Directory where plot files will be saved
    output_filename_base : str
        Base filename for the output plots (without extension)
    """
    # Determine plotting strategy based on data size
    total_rows = len(data_df)

    if total_rows <= 50:
        # Show all data for reasonable-sized datasets
        plot_data = data_df.copy()
        x_labels = data_df["num_mutations"].astype(str).tolist()
        show_gap_info = False
        logger.info(f"Showing all {total_rows} data points")
    else:
        # Show first 30 and last 15 for very large datasets
        first_30 = data_df.head(30)
        last_15 = data_df.tail(15)

        # Create combined data
        plot_data = pd.concat([first_30, last_15], ignore_index=True)

        # Create labels with gap indicator
        first_labels = first_30["num_mutations"].astype(str).tolist()
        last_labels = last_15["num_mutations"].astype(str).tolist()
        x_labels = first_labels + ["..."] + last_labels

        show_gap_info = True
        skipped = total_rows - 45
        logger.info(
            f"Showing first 30 and last 15 data points, skipping {skipped} middle points"
        )

    # Create the plot
    plt.figure(figsize=(18, 10))

    # Create bar positions
    if show_gap_info:
        # Positions: 0-29 for first part, gap at 30, 31-45 for last part
        first_positions = list(range(30))
        last_positions = list(range(31, 46))
        all_positions = first_positions + [30] + last_positions

        # Values: first 30 values, 0 for gap, last 15 values
        first_values = plot_data["num_sequences"].head(30).tolist()
        last_values = plot_data["num_sequences"].tail(15).tolist()
        all_values = first_values + [0] + last_values

        # Create bars (viridis monochrome)
        bars = plt.bar(all_positions, all_values, width=0.7, color=mono_color())

        # Make gap bar invisible
        bars[30].set_alpha(0)

    else:
        # Simple plotting for manageable datasets
        positions = range(len(plot_data))
        bars = plt.bar(
            positions, plot_data["num_sequences"], width=0.7, color=mono_color()
        )
        all_positions = list(positions)

    # Format the plot
    plt.title(
        "Distribution of Mutation Counts in ChEMBL Variant Sequences",
        fontsize=32,
        pad=20,
    )
    plt.xlabel("Number of Mutations per Sequence", fontsize=24)
    plt.ylabel("Number of Sequences", fontsize=24)

    # Set x-axis labels
    plt.xticks(all_positions, x_labels, ha="center", fontsize=20)

    # Set logarithmic scale for y-axis
    plt.yscale("log")

    # Add grid
    plt.grid(True, axis="y", linestyle="--", alpha=0.7)
    # Increase tick label font sizes for publication readability
    plt.tick_params(axis="both", which="major", labelsize=20)

    # Add value labels on all bars
    bar_index = 0
    for i, pos in enumerate(all_positions):
        if show_gap_info and i == 30:  # Skip gap position
            continue

        if bar_index < len(plot_data):
            value = plot_data["num_sequences"].iloc[bar_index]
            plt.text(
                pos,
                value * 1.02,
                _format_number_latex_style(value),
                ha="center",
                va="bottom",
                fontsize=16,
            )
            bar_index += 1

    # Add statistics box
    total_sequences = data_df["num_sequences"].sum()
    max_mutations = data_df["num_mutations"].max()

    # Calculate mean and median from full dataset
    expanded = []
    for _, row in data_df.iterrows():
        expanded.extend([row["num_mutations"]] * int(row["num_sequences"]))
    mean_mutations = np.mean(expanded)
    median_mutations = np.median(expanded)

    stats_text = f"Total sequences: {_format_number_latex_style(total_sequences)}\n"
    if show_gap_info:
        stats_text += f"Showing: 45 of {total_rows} mutation counts\n"
        stats_text += f"(First 30 + Last 15)\n"
    stats_text += f"Max mutations: {int(max_mutations)}\n"
    stats_text += f"Mean mutations: {mean_mutations:.1f}\n"
    stats_text += f"Median mutations: {median_mutations:.1f}"

    plt.text(
        0.985,
        0.97,
        stats_text,
        transform=plt.gca().transAxes,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.9),
        verticalalignment="top",
        horizontalalignment="right",
        fontsize=20,
    )

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)

    # Save files
    figures_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    png_path = figures_dir / f"{output_filename_base}.png"
    pdf_path = figures_dir / f"{output_filename_base}.pdf"
    plt.savefig(png_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.savefig(pdf_path, bbox_inches="tight", facecolor="white")

    plt.close()  # Close the figure to free memory

    logger.info(f"Plot saved to {png_path}")
    logger.info(f"PDF plot saved to {pdf_path}")


def execute_mutation_analysis(
    db_interaction,
    output_dir: pathlib.Path,
    output_filename_base: str = "mutation_distribution",
    generate_latex: bool = True,
    create_plots: bool = True,
    sql_file_path: Optional[pathlib.Path] = None,
) -> pd.DataFrame:
    """
    Execute the mutation count distribution analysis: run SQL query, generate plots and LaTeX table.

    Parameters
    ----------
    db_interaction : UnifiedPostgresConnector
        Unified PostgreSQL connector instance configured with connection parameters
    output_dir : pathlib.Path
        Directory where all output files will be saved
    output_filename_base : str
        Base filename for output files
    generate_latex : bool
        Whether to generate LaTeX table
    create_plots : bool
        Whether to create histogram plots
    sql_file_path : pathlib.Path, optional
        Path to a custom SQL file. If None, uses the embedded VARIANT_MUTATION_DATA_EXTRACTION_SQL constant.

    Returns
    -------
    pd.DataFrame
        The resulting data DataFrame
    """

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get SQL query: use custom file if provided, otherwise use embedded constant
    if sql_file_path is not None:
        with open(sql_file_path, "r") as f:
            sql_query = f.read()
        logger.info(f"Executing SQL query from custom file: {sql_file_path}")
    else:
        sql_query = VARIANT_MUTATION_DATA_EXTRACTION_SQL
        logger.info("Executing embedded variant mutation data extraction SQL query")

    try:
        # Execute query using the unified connector interface
        results = db_interaction.execute(sql_query) or []

        # Convert results to DataFrame
        if results:
            # Extract column names - we know the expected columns from the SQL
            df = pd.DataFrame(
                results, columns=["num_mutations", "num_sequences", "percentage"]
            )

            # Convert data types
            df["num_mutations"] = df["num_mutations"].astype(int)
            df["num_sequences"] = df["num_sequences"].astype(int)
            df["percentage"] = df["percentage"].astype(float)
        else:
            # Empty results
            df = pd.DataFrame(columns=["num_mutations", "num_sequences", "percentage"])

        # Save the DataFrame as CSV
        tables_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        csv_path = tables_dir / f"{output_filename_base}_data.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"Query results saved to {csv_path}")

        # Generate plots if requested
        if create_plots:
            plot_mutation_distribution(
                data_df=df,
                output_dir=output_dir,
                output_filename_base=f"{output_filename_base}_plot",
            )

        # Generate LaTeX table if requested
        if generate_latex:
            latex_path = tables_dir / f"{output_filename_base}_table.tex"
            generate_latex_table(
                df=df,
                output_path=latex_path,
                table_caption="Distribution of Mutation Counts in ChEMBL Variant Sequences",
                table_label="tab:mutation_count_distribution_variant_sequences",
            )

        return df

    except Exception as e:
        logger.exception("Failed to execute SQL query: {}", e)
        sys.exit(1)
