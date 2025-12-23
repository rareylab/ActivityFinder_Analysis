import pathlib
import sys
from datetime import date
from typing import Optional

import typer
from loguru import logger

from activity_db_analysis.apps.ligand_match_level_counts import (
    run as run_ligand_match_level_counts,
)
from activity_db_analysis.db.cli_tunnel_wrapper import TunnelMode
from activity_db_analysis.db.ssh_tunnel import (
    SSHTunnelManager,
    TunnelConfig,
    get_tunnel_status,
)
from activity_db_analysis.siteconfig.app_config import (
    AppConfig,
    build_chembl_sql_handler_from_config,
    parse_app_config_toml,
)
from activity_db_analysis.siteconfig.logging import LoggingLevel

app = typer.Typer(
    name="Activity-DB-Analysis",
    help="Suites of tools to reproduce the analyses around the ActivityFinder publication.",
    epilog=f"(C) {date.today().year}",
    add_completion=False,
    invoke_without_command=True,
)

# ---------------------------------------------------------------------------
# SQL Artifacts Management Sub-Command Integration
# ---------------------------------------------------------------------------
try:
    # Import the sql artifacts management Typer app and mount it as a sub-command.
    # The mounted name 'sql-manager' enables: uv run actfind sql-manager <command>
    from activity_db_analysis.apps.sql_artifacts import (
        sql_artifacts_management as _sql_mgr_mod,  # type: ignore
    )

    app.add_typer(
        _sql_mgr_mod.app,
        name="sql-manager",
        help="Manage SQL artifacts (create, drop, refresh, list, validate, etc.)",
    )
except Exception as e:  # pragma: no cover
    # Fail soft so the core CLI still works even if optional module has an issue.
    logger.warning(f"Could not register sql-manager sub-command: {e}")

# ---------------------------------------------------------------------------
# Sequence Mapping Sub-Command Integration
# ---------------------------------------------------------------------------
try:
    from activity_db_analysis.apps.sequence_mapping import (
        sequence_mapping_cli as _seq_map_mod,  # type: ignore
    )

    app.add_typer(
        _seq_map_mod.app,
        name="sequence-mapper",
        help="Sequence mapping reports & plots (generate right-bound & interval metrics)",
    )
except Exception as e:  # pragma: no cover
    logger.warning(f"Could not register sequence-mapper sub-command: {e}")


# ---------------------------------------------------------------------------
# Global callback to load an app configuration TOML before executing commands.
# Provides a unified --app-config option available to all subcommands.
# Now also manages SSH tunnel for database connectivity.
# ---------------------------------------------------------------------------
@app.callback()
def cli_root_callback(
    ctx: typer.Context,
    app_config: Optional[pathlib.Path] = typer.Option(
        None,
        "--app-config",
        help="Optional path to configuration TOML. If omitted, will attempt ./appconfig.toml if present.",
        exists=True,
        dir_okay=False,
        readable=True,
    ),
    logging_level: LoggingLevel = typer.Option(
        LoggingLevel.INFO,
        "--logging-level",
        help="Logging verbosity (TRACE, DEBUG, INFO, SUCCESS, WARNING, ERROR, CRITICAL).",
    ),
    manage_tunnel: bool = typer.Option(
        True,
        "--manage-tunnel/--no-tunnel",
        help="Automatically manage SSH tunnel for database access. Use --no-tunnel for local development or when tunnel is externally managed.",
    ),
    tunnel_mode: TunnelMode = typer.Option(
        TunnelMode.AUTO,
        "--tunnel-mode",
        help="Tunnel management mode: auto (start if needed, don't stop external), required (fail if unavailable), optional (continue without), skip (no management).",
    ),
    show_tunnel_status: bool = typer.Option(
        False,
        "--tunnel-status",
        help="Show current tunnel status and exit.",
    ),
) -> None:
    logger.remove()
    logger.add(sys.stderr, level=logging_level.value)

    # Handle --tunnel-status flag
    if show_tunnel_status:
        status = get_tunnel_status()
        if status["is_active"]:
            logger.info(
                "✓ SSH tunnel ACTIVE: {} → {}",
                status["local_endpoint"],
                status["remote_endpoint"],
            )
            if status["owned_by_manager"]:
                logger.info("  (managed by this process)")
            else:
                logger.info("  (externally managed)")
        else:
            logger.warning(
                "✗ SSH tunnel INACTIVE: {} → {}",
                status["local_endpoint"],
                status["remote_endpoint"],
            )
            logger.info(
                "  Start with: ssh -f -N {}",
                status["ssh_host_alias"],
            )
        raise typer.Exit(0)

    # Load configuration
    cfg_path: Optional[pathlib.Path] = app_config
    if cfg_path is None:
        default_path = pathlib.Path("./appconfig.toml")
        if default_path.exists():
            cfg_path = default_path

    if cfg_path is None:
        logger.info(
            "No configuration file provided and default ./appconfig.toml not found."
        )
    else:
        config_dict = parse_app_config_toml(cfg_path)
        if config_dict:
            AppConfig.set_config_file_parameters(config_dict)
            logger.info(f"Loaded application configuration from {cfg_path}")
        else:
            logger.warning(f"Failed to load configuration from {cfg_path}")

    # Manage SSH tunnel based on options
    if manage_tunnel and tunnel_mode != TunnelMode.SKIP:
        # Build TunnelConfig from appconfig.toml [ssh] section if available
        tunnel_config = None
        if config_dict and "ssh" in config_dict:
            ssh_cfg = config_dict["ssh"]
            tunnel_config = TunnelConfig(
                ssh_host_alias=ssh_cfg.get("host_alias", "db_tunnel"),
                local_host=ssh_cfg.get("local_bind_host", "127.0.0.1"),
                local_port=ssh_cfg.get("local_bind_port", 5433),
                remote_host=ssh_cfg.get("remote_bind_host", "database.example.com"),
                remote_port=ssh_cfg.get("remote_bind_port", 5432),
            )
        tunnel_manager = SSHTunnelManager(config=tunnel_config)

        try:
            if tunnel_mode == TunnelMode.OPTIONAL:
                try:
                    tunnel_manager.start()
                    logger.debug("SSH tunnel started (optional mode)")
                except Exception as e:
                    logger.warning(
                        f"Could not establish SSH tunnel (continuing anyway): {e}"
                    )
            elif tunnel_mode in (TunnelMode.AUTO, TunnelMode.REQUIRED):
                tunnel_manager.start()
                if tunnel_manager.is_tunnel_active:
                    if tunnel_manager.owns_tunnel:
                        logger.info("SSH tunnel established (will be stopped on exit)")
                    else:
                        logger.info("Using existing SSH tunnel (externally managed)")
        except Exception as e:
            if tunnel_mode == TunnelMode.REQUIRED:
                logger.exception("SSH tunnel required but failed to start: {}", e)
                raise typer.Exit(1)
            elif tunnel_mode == TunnelMode.AUTO:
                logger.warning(f"Could not establish SSH tunnel: {e}")
                raise

        # Store tunnel manager in context for potential cleanup
        ctx.obj = ctx.obj or {}
        ctx.obj["tunnel_manager"] = tunnel_manager


@app.command(no_args_is_help=False)
def pair_identity_distribution_activity(
    output_dir: pathlib.Path = typer.Option(
        None,
        "--output-dir",
        help="Directory to write activity-qualified pair identity distribution outputs (CSV + optional LaTeX). Defaults to thesis tables directory.",
        resolve_path=True,
    ),
    latex: bool = typer.Option(
        True,
        "--latex/--no-latex",
        help="Generate LaTeX output (default: True).",
    ),
    cumulative: bool = typer.Option(
        False,
        "--cumulative",
        help="Count each pair in every identity interval it qualifies for (multi-bucket cumulative mode).",
    ),
):
    """
    Run the activity-qualified, coverage-agnostic ligand pair identity distribution analysis.

    NOTE: SSH tunnel is managed automatically by the CLI. Use --no-tunnel to disable.
    """
    from activity_db_analysis.apps.ligand_mapping.pair_identity_distribution_activity import (
        run as _run_pair_identity_distribution_activity,
    )

    _run_pair_identity_distribution_activity(
        output_dir=output_dir,
        no_tex=not latex,
        cumulative=cumulative,
    )


@app.command(no_args_is_help=False)
def plot_pair_identity_distribution_activity(
    input_csv: pathlib.Path = typer.Option(
        ...,
        "--input",
        exists=True,
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
        help="Path to disjoint pair_identity_distribution_activity.csv (max-identity assignment).",
    ),
    output_dir: pathlib.Path = typer.Option(
        ...,
        "--output-dir",
        resolve_path=True,
        help="Directory to write generated figure files.",
    ),
    output_prefix: str = typer.Option(
        "pair_identity_distribution_activity",
        "--output-prefix",
        help="Base filename prefix for output figures.",
    ),
    no_counts: bool = typer.Option(
        False,
        "--no-counts",
        help="Do not generate counts stacked bar plot.",
    ),
    no_percents: bool = typer.Option(
        False,
        "--no-percents",
        help="Do not generate normalized (within-interval composition) stacked bar plot.",
    ),
    no_annotate: bool = typer.Option(
        False,
        "--no-annotate",
        help="Disable total / 100% bar-top annotations.",
    ),
    dpi: int = typer.Option(
        250,
        "--dpi",
        help="DPI for raster (PNG) outputs.",
    ),
    no_png: bool = typer.Option(
        False,
        "--no-png",
        help="Skip PNG output (use with --pdf to produce only PDF).",
    ),
    pdf: bool = typer.Option(
        False,
        "--pdf",
        help="Also (or only if --no-png) write PDF output(s).",
    ),
    combined: bool = typer.Option(
        False,
        "--combined",
        help="Generate a combined 2x2 panel (disjoint vs cumulative) in addition to standard plots.",
    ),
    cumulative_input: Optional[pathlib.Path] = typer.Option(
        None,
        "--cumulative-input",
        exists=True,
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
        help="Path to cumulative pair_identity_distribution_activity_cumulative.csv (required if --combined).",
    ),
):
    """
    Plot the activity-qualified pair identity distribution (counts and/or normalized within-interval composition),
    stacked by structural confidence (smmclid).
    """
    from activity_db_analysis.apps.ligand_mapping.plot_pair_identity_distribution_activity import (
        _combined_panel,
        _derive_percents_from_counts,
        _pivot_counts,
        _read_input,
        create_pair_identity_distribution_activity_plots,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    created = create_pair_identity_distribution_activity_plots(
        input_csv=input_csv,
        output_dir=output_dir,
        output_prefix=output_prefix,
        counts=not no_counts,
        percents=not no_percents,
        annotate=not no_annotate,
        dpi=dpi,
        save_png=not no_png,
        save_pdf=pdf,
    )

    if combined:
        if cumulative_input is None:
            raise typer.BadParameter(
                "--combined requires --cumulative-input pointing to cumulative CSV."
            )

        disjoint_df = _read_input(input_csv)
        cumulative_df = _read_input(cumulative_input)

        d_counts, d_totals = _pivot_counts(disjoint_df)
        d_norm = _derive_percents_from_counts(d_counts)

        c_counts, c_totals = _pivot_counts(cumulative_df)
        c_norm = _derive_percents_from_counts(c_counts)

        combined_path = output_dir / f"{output_prefix}_combined_panel.png"
        _combined_panel(
            disjoint_counts=d_counts,
            disjoint_norm=d_norm,
            disjoint_interval_totals=d_totals,
            cumulative_counts=c_counts,
            cumulative_norm=c_norm,
            cumulative_interval_totals=c_totals,
            output_path=combined_path,
            dpi=dpi,
            save_pdf=pdf,
        )
        created.append(combined_path)
        logger.info(
            "Generated combined disjoint vs cumulative panel: {}", combined_path
        )

    if created:
        logger.info(
            f"Generated {len(created)} activity-qualified pair identity distribution figure file(s)."
        )
    else:
        logger.warning("No figures generated (input may be empty).")


@app.command(no_args_is_help=False)
def plot_het_code_variants(
    sql_file: Optional[pathlib.Path] = typer.Option(
        None,
        "--sql-file",
        help="Path to custom SQL file. If not provided, uses embedded SQL query.",
        resolve_path=True,
    ),
    output_dir: pathlib.Path = typer.Option(
        pathlib.Path("./het_code_analysis"),
        help="Directory where output files will be saved",
        resolve_path=True,
    ),
    output_prefix: str = typer.Option(
        "variants_per_het_code",
        "--output-prefix",
        help="Base filename prefix for output files (plots and data)",
    ),
    latex: bool = typer.Option(
        True,
        "--latex/--no-latex",
        help="Generate LaTeX table from the data (default: True)",
    ),
    create_plots: bool = typer.Option(
        True,
        help="Create histogram plots from the data",
    ),
    dpi: int = typer.Option(
        250,
        "--dpi",
        help="DPI for raster (PNG) outputs.",
    ),
    pdf: bool = typer.Option(
        False,
        "--pdf",
        help="Also write PDF output(s).",
    ),
    show: bool = typer.Option(
        False,
        "--show/--no-show",
        help="Display the plot interactively.",
    ),
) -> None:
    """
    Generate histogram plots and LaTeX table for HET code variants distribution.

    By default, uses an embedded SQL query. Use --sql-file to provide a custom query.
    """
    from activity_db_analysis.apps.het_code_plot.het_code_variants import (
        execute_het_code_analysis,
    )

    # Validate config is loaded (execute_het_code_analysis will use AppConfig internally)
    AppConfig.set_activity_db_database_connection_parameters(cli_options={})

    try:
        df = execute_het_code_analysis(
            output_dir=output_dir,
            output_filename_base=output_prefix,
            generate_latex=latex,
            create_plots=create_plots,
            sql_file_path=sql_file,
        )

        logger.info(
            f"Analysis completed successfully. Found {len(df)} unique variant counts."
        )
        logger.info(f"Output files saved in: {output_dir}")

    except Exception as e:
        logger.exception("Failed to execute HET code analysis: {}", e)
        raise typer.Exit(1)


@app.command(no_args_is_help=False)
def plot_variant_mutations(
    sql_file: Optional[pathlib.Path] = typer.Option(
        None,
        "--sql-file",
        help="Path to custom SQL file. If not provided, uses embedded SQL query.",
        resolve_path=True,
    ),
    output_dir: pathlib.Path = typer.Option(
        pathlib.Path("./variant_mutation_analysis"),
        help="Directory where output files will be saved",
        resolve_path=True,
    ),
    output_prefix: str = typer.Option(
        "mutation_distribution",
        "--output-prefix",
        help="Base filename prefix for output files (plots and data)",
    ),
    latex: bool = typer.Option(
        True,
        "--latex/--no-latex",
        help="Generate LaTeX table from the data (default: True)",
    ),
    create_plots: bool = typer.Option(
        True,
        help="Create bar plots from the data",
    ),
    dpi: int = typer.Option(
        250,
        "--dpi",
        help="DPI for raster (PNG) outputs.",
    ),
    pdf: bool = typer.Option(
        False,
        "--pdf",
        help="Also write PDF output(s).",
    ),
    show: bool = typer.Option(
        False,
        "--show/--no-show",
        help="Display the plot interactively.",
    ),
) -> None:
    """
    Generate bar plots and LaTeX table for variant mutation count distribution.

    By default, uses an embedded SQL query. Use --sql-file to provide a custom query.
    """
    from activity_db_analysis.apps.variant_mutation_plot.variant_mutation_analysis import (
        execute_mutation_analysis,
    )

    chembl_db_params: AppConfig.ChemblDatabaseConnectionParameters = (
        AppConfig.set_chembl_db_database_connection_parameters(cli_options={})
    )

    required_params = ["name", "user", "password", "host", "port", "db_schema"]
    missing_params = [
        k for k in required_params if not getattr(chembl_db_params, k, None)
    ]
    if missing_params:
        logger.error(
            f"Missing required ChEMBL database parameters in config: {', '.join(missing_params)}"
        )
        logger.info(
            "Please ensure these parameters are set in the [chembl_database] section of appconfig.toml"
        )
        raise typer.Exit(1)

    chembl_db_interaction = build_chembl_sql_handler_from_config(
        cli_overrides={
            "name": chembl_db_params.name,
            "user": chembl_db_params.user,
            "password": chembl_db_params.password,
            "host": chembl_db_params.host,
            "port": chembl_db_params.port,
            "schema": chembl_db_params.db_schema,
        }
    )

    try:
        df = execute_mutation_analysis(
            db_interaction=chembl_db_interaction,
            output_dir=output_dir,
            output_filename_base=output_prefix,
            generate_latex=latex,
            create_plots=create_plots,
            sql_file_path=sql_file,
        )

        logger.info(
            f"Analysis completed successfully. Found {len(df)} unique mutation counts."
        )
        logger.info(f"Output files saved in: {output_dir}")

    except Exception as e:
        logger.exception("Failed to execute variant mutation analysis: {}", e)
        raise typer.Exit(1)


@app.command(no_args_is_help=False)
def plot_amino_acid_mapping(
    input_csv: pathlib.Path = typer.Option(
        ...,
        "--input",
        help="Path to CSV file containing three-letter to one-letter amino acid mappings.",
        resolve_path=True,
    ),
    output_dir: pathlib.Path = typer.Option(
        pathlib.Path("./csv_plot_results"),
        help="Directory where output files (PNG/PDF/CSV) will be saved",
        resolve_path=True,
    ),
    output_prefix: str = typer.Option(
        "amino_acid_distribution",
        "--output-prefix",
        help="Base filename prefix for output files (plots and summary CSV)",
    ),
    delimiter: str = typer.Option(
        ";",
        help="Delimiter used in the input CSV file",
    ),
    show: bool = typer.Option(
        False,
        "--show/--no-show",
        help="Display the plot interactively",
    ),
    write_csv_summary: bool = typer.Option(
        True,
        help="Write a CSV summary of counts by one-letter code",
    ),
    latex: bool = typer.Option(
        True,
        "--latex/--no-latex",
        help="Generate LaTeX longtable from the mapping CSV (default: True)",
    ),
    latex_only: bool = typer.Option(
        False,
        "--latex-only",
        help="Only generate LaTeX table, skip plot generation",
    ),
    latex_output: Optional[pathlib.Path] = typer.Option(
        None,
        "--latex-output",
        help="Custom path for LaTeX output file (.tex). Defaults to <output_dir>/tables/<output_prefix>.tex",
        resolve_path=True,
    ),
    dpi: int = typer.Option(
        250,
        "--dpi",
        help="DPI for raster (PNG) outputs.",
    ),
    pdf: bool = typer.Option(
        False,
        "--pdf",
        help="Also write PDF output(s).",
    ),
) -> None:
    """
    Create a colored bar chart (canonical vs non-canonical) showing how many three-letter
    identifiers map to each single-letter amino acid code. Optionally writes a summary CSV
    and/or a LaTeX longtable.
    """
    try:
        from activity_db_analysis.apps.amino_acid_mapping import (
            plot_amino_acid_distribution,
        )
        from activity_db_analysis.apps.amino_acid_mapping.amino_acid_mapping import (
            generate_latex_table,
        )

        written: list[pathlib.Path] = []

        # Generate plot unless --latex-only
        if not latex_only:
            plot_outputs = plot_amino_acid_distribution(
                input_csv=input_csv,
                output_dir=output_dir,
                output_filename_base=output_prefix,
                delimiter=delimiter,
                write_csv_summary=write_csv_summary,
            )
            written.extend(plot_outputs)

        # Generate LaTeX table if requested
        if latex or latex_only:
            if latex_output is not None:
                tex_path = latex_output
            else:
                tables_dir = output_dir / "tables"
                tables_dir.mkdir(parents=True, exist_ok=True)
                tex_path = tables_dir / f"{output_prefix}.tex"

            latex_path = generate_latex_table(
                input_csv=input_csv,
                output_path=tex_path,
                delimiter=delimiter,
            )
            written.append(latex_path)

        for p in written:
            logger.info(f"Wrote: {p}")
        logger.info(f"Amino acid mapping completed. Outputs in: {output_dir}")

    except Exception as e:
        logger.exception("Failed to generate amino acid mapping: {}", e)
        raise typer.Exit(1)


@app.command(no_args_is_help=False)
def ligand_match_level_counts(
    output_dir: pathlib.Path = typer.Option(
        None,
        "--output-dir",
        help="Directory for CSV/TeX ligand match level outputs (defaults to thesis tables directory).",
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
    ),
    latex: bool = typer.Option(
        True,
        "--latex/--no-latex",
        help="Generate LaTeX table (default: True).",
    ),
    output_prefix: str = typer.Option(
        "ligand_match_levels.csv",
        "--output-prefix",
        help="Output CSV filename prefix (TeX will share the same stem).",
    ),
) -> None:
    """
    Generate counts of PDB–ChEMBL ligand mapping confidence levels (already best-per-pair).
    Writes one CSV (and optional LaTeX table) to the specified output directory.

    NOTE: SSH tunnel is managed automatically by the CLI. Use --no-tunnel to disable.
    """
    run_ligand_match_level_counts(
        output_dir=output_dir,
        raw_basename=output_prefix,
        write_tex=latex,
    )


@app.command(name="tunnel-status", no_args_is_help=False)
def tunnel_status_command() -> None:
    """
    Display the current SSH tunnel status.

    Shows whether the tunnel is active, the configured endpoints, and whether
    the tunnel is managed by this process or externally.
    """
    status = get_tunnel_status()

    typer.echo("\n=== SSH Tunnel Status ===\n")
    typer.echo(f"SSH Host Alias:    {status['ssh_host_alias']}")
    typer.echo(f"Local Endpoint:    {status['local_endpoint']}")
    typer.echo(f"Remote Endpoint:   {status['remote_endpoint']}")
    typer.echo("")

    if status["is_active"]:
        typer.secho("Status: ✓ ACTIVE", fg=typer.colors.GREEN, bold=True)
        if status["owned_by_manager"]:
            typer.echo("        (managed by this process)")
        else:
            typer.echo("        (externally managed)")
    else:
        typer.secho("Status: ✗ INACTIVE", fg=typer.colors.RED, bold=True)
        typer.echo(f"\nTo start tunnel manually:")
        typer.echo(f"  ssh -f -N {status['ssh_host_alias']}")

    typer.echo("")


if __name__ == "__main__":
    app()
