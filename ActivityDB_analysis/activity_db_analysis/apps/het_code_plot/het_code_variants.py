import pathlib
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger
import sys

from activity_db_analysis.plotting.palettes import mono_color
from activity_db_analysis.apps.het_code_plot.het_code_sql import HET_CODE_DATA_EXTRACTION_SQL


def generate_latex_table(df: pd.DataFrame, output_path: pathlib.Path, table_caption: str = "Distribution of Ligands by Number of Variants", table_label: str = "tab:ligand_variant_distribution") -> None:
    """
    Generate a multipage LaTeX table from the HET code variants data using longtable and siunitx.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns: num_variants, num_ligands, percentage
    output_path : pathlib.Path
        Path where the LaTeX table file will be saved
    table_caption : str
        Caption for the LaTeX table
    table_label : str
        Label for the LaTeX table reference

    Note
    ----
    Requires \\usepackage{longtable} and \\usepackage{siunitx} in LaTeX preamble.
    """
    latex_content = r"% Requires \usepackage{longtable} and \usepackage{siunitx}" + "\n"
    latex_content += r"\begin{longtable}{rrr}" + "\n"
    latex_content += f"    \\caption{{{table_caption}}} \\label{{{table_label}}} \\\\\n"
    latex_content += r"    \toprule" + "\n"
    latex_content += r"    \textbf{Number of variants} & \textbf{Number of HET codes} & \textbf{Percentage (\%)} \\" + "\n"
    latex_content += r"    \midrule" + "\n"
    latex_content += r"    \endfirsthead" + "\n"
    latex_content += r"    \multicolumn{3}{c}{\textit{Continued from previous page}} \\" + "\n"
    latex_content += r"    \toprule" + "\n"
    latex_content += r"    \textbf{Number of variants} & \textbf{Number of HET codes} & \textbf{Percentage (\%)} \\" + "\n"
    latex_content += r"    \midrule" + "\n"
    latex_content += r"    \endhead" + "\n"
    latex_content += r"    \midrule" + "\n"
    latex_content += r"    \multicolumn{3}{c}{\textit{Continued on next page}} \\" + "\n"
    latex_content += r"    \endfoot" + "\n"
    latex_content += r"    \bottomrule" + "\n"
    latex_content += r"    \endlastfoot" + "\n"

    # Show all individual values using siunitx formatting
    for _, row in df.iterrows():
        latex_content += f"    {int(row['num_variants'])} & \\num{{{int(row['num_ligands'])}}} & {row['percentage']:.3f} \\\\\n"

    latex_content += r"\end{longtable}"

    # Write to file
    with open(output_path, 'w') as f:
        f.write(latex_content)

    logger.info(f"LaTeX table saved to {output_path}")


def plot_het_code_variants(
    data_df: pd.DataFrame,
    output_dir: pathlib.Path,
    output_filename_base: str = "variants_per_het_code_histogram"
) -> None:
    """
    Create histogram plots for HET code variants distribution with simplified approach.

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame with columns: num_variants, num_ligands, percentage
    output_dir : pathlib.Path
        Directory where plot files will be saved
    output_filename_base : str
        Base filename for the output plots (without extension)
    """
    import numpy as np

    # Determine what to plot
    total_rows = len(data_df)

    if total_rows <= 35:
        # Show all data for small datasets
        plot_data = data_df.copy()
        x_labels = data_df['num_variants'].astype(str).tolist()
        show_gap_info = False
        logger.info(f"Showing all {total_rows} data points")
    else:
        # Show first 25 and last 10 for large datasets
        first_25 = data_df.head(25)
        last_10 = data_df.tail(10)

        # Create combined data
        plot_data = pd.concat([first_25, last_10], ignore_index=True)

        # Create labels with gap indicator
        first_labels = first_25['num_variants'].astype(str).tolist()
        last_labels = last_10['num_variants'].astype(str).tolist()
        x_labels = first_labels + ['...'] + last_labels

        show_gap_info = True
        skipped = total_rows - 35
        logger.info(f"Showing first 25 and last 10 data points, skipping {skipped} middle points")

    # Create the plot
    plt.figure(figsize=(18, 10))

    # Create bar positions
    if show_gap_info:
        # Positions: 0-24 for first part, gap at 25, 26-35 for last part
        first_positions = list(range(25))
        last_positions = list(range(26, 36))
        all_positions = first_positions + [25] + last_positions

        # Values: first 25 values, 0 for gap, last 10 values
        first_values = plot_data['num_ligands'].head(25).tolist()
        last_values = plot_data['num_ligands'].tail(10).tolist()
        all_values = first_values + [0] + last_values

        # Create bars (viridis monochrome)
        bars = plt.bar(all_positions, all_values, width=0.7, color=mono_color())

        # Make gap bar invisible
        bars[25].set_alpha(0)

    else:
        # Simple plotting for small datasets
        positions = range(len(plot_data))
        bars = plt.bar(positions, plot_data['num_ligands'], width=0.7, color=mono_color())
        all_positions = list(positions)

    # Format the plot
    plt.title('Variants per HET Code Distribution', fontsize=32, pad=20)
    plt.xlabel('Number of Variants per HET Code', fontsize=24)
    plt.ylabel('Number of HET Codes', fontsize=24)

    # Set y-axis to log scale
    plt.yscale('log')

    # Set x-axis labels
    plt.xticks(all_positions, x_labels, rotation=45, ha='center', fontsize=20)

    # Add grid
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=20)

    # Add value labels on all bars
    bar_index = 0
    for i, pos in enumerate(all_positions):
        if show_gap_info and i == 25:  # Skip gap position
            continue

        if bar_index < len(plot_data):
            value = plot_data['num_ligands'].iloc[bar_index]
            plt.text(pos, value * 1.1, f'{int(value)}',
                    ha='center', va='bottom', fontsize=16)
            bar_index += 1

    # Add statistics box
    total_het_codes = data_df['num_ligands'].sum()
    max_variants = data_df['num_variants'].max()

    # Calculate median from full dataset
    expanded = []
    for _, row in data_df.iterrows():
        expanded.extend([row['num_variants']] * int(row['num_ligands']))
    median_variants = np.median(expanded)

    stats_text = f'Total HET codes: {total_het_codes:,}\n'
    if show_gap_info:
        stats_text += f'Showing: 35 of {total_rows} variant counts\n'
        stats_text += f'(First 25 + Last 10)\n'
    stats_text += f'Max variants: {int(max_variants)}\n'
    stats_text += f'Median variants: {median_variants:.1f}\n'
    stats_text += 'Basis: InChI-Key'

    plt.text(0.985, 0.97, stats_text, transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9),
             verticalalignment='top', horizontalalignment='right',
             fontsize=20)

    # (Removed separate variant definition annotation; info integrated into stats box.)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)

    # Save files
    figures_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    png_path = figures_dir / f'{output_filename_base}.png'
    pdf_path = figures_dir / f'{output_filename_base}.pdf'

    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')

    plt.close()  # Close the figure to free memory

    logger.info(f"Main plot saved to {png_path}")
    logger.info(f"PDF plot saved to {pdf_path}")


def execute_het_code_analysis(
    output_dir: pathlib.Path,
    output_filename_base: str = "variants_per_het_code",
    generate_latex: bool = True,
    create_plots: bool = True,
    sql_file_path: Optional[pathlib.Path] = None,
) -> pd.DataFrame:
    """
    Execute the HET code variants analysis: run SQL query, generate plots and LaTeX table.

    Database connection is obtained via AppConfig (appconfig.toml + environment variables).
    SSH tunneling is managed at the CLI layer (--manage-tunnel, --tunnel-mode).

    Parameters
    ----------
    output_dir : pathlib.Path
        Directory where all output files will be saved
    output_filename_base : str
        Base filename for output files
    generate_latex : bool
        Whether to generate LaTeX table
    create_plots : bool
        Whether to create histogram plots
    sql_file_path : pathlib.Path, optional
        Path to a custom SQL file. If None, uses the embedded HET_CODE_DATA_EXTRACTION_SQL constant.

    Returns
    -------
    pd.DataFrame
        The resulting data DataFrame
    """
    from activity_db_analysis.siteconfig.app_config import AppConfig

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get SQL query: use custom file if provided, otherwise use embedded constant
    if sql_file_path is not None:
        with open(sql_file_path, 'r') as f:
            sql_query = f.read()
        logger.info(f"Executing SQL query from custom file: {sql_file_path}")
    else:
        sql_query = HET_CODE_DATA_EXTRACTION_SQL
        logger.info("Executing embedded HET code data extraction SQL query")

    try:
        # Database connection via AppConfig + env variables
        with AppConfig.build_activity_sql_handler_from_config() as db:
            results = db.execute(sql_query) or []

        # Convert results to DataFrame
        # The first row should contain column headers
        if results:
            # Extract column names from the SQL query or use generic names
            # Since we know the expected columns, we can hardcode them
            df = pd.DataFrame(results, columns=['num_variants', 'num_ligands', 'percentage'])

            # Convert data types
            df['num_variants'] = df['num_variants'].astype(int)
            df['num_ligands'] = df['num_ligands'].astype(int)
            df['percentage'] = df['percentage'].astype(float)
        else:
            # Empty results
            df = pd.DataFrame(columns=['num_variants', 'num_ligands', 'percentage'])

        # Save the DataFrame as CSV to tables directory
        tables_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        csv_path = tables_dir / f'{output_filename_base}_data.csv'
        df.to_csv(csv_path, index=False)
        logger.info(f"Query results saved to {csv_path}")

        # Generate plots if requested
        if create_plots:
            plot_het_code_variants(
                data_df=df,
                output_dir=output_dir,
                output_filename_base=f'{output_filename_base}_histogram'
            )

        # Generate LaTeX table if requested
        if generate_latex:
            tables_dir = (output_dir if output_dir else pathlib.Path.cwd()) / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)
            latex_path = tables_dir / f'{output_filename_base}_table.tex'
            generate_latex_table(
                df=df,
                output_path=latex_path,
                table_caption="Distribution of HET Codes by Number of Variants (variants = distinct raw NAOMI InChI-Keys per HET code; based on standard_inchi_key)",
                table_label="tab:het_code_variant_distribution"
            )

        return df

    except Exception as e:
        logger.exception("Failed to execute SQL query: {}", e)
        sys.exit(1)
