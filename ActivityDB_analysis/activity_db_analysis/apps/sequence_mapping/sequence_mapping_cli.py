"""
Sequence Mapping Typer Sub-Application
======================================

Provides a unified Typer-based CLI interface for generating sequence mapping
reports (right-bound + interval + simplified publication tables) leveraging
the existing business logic in `sequence_mapping.generate_target_matching_reports`.

Scope (initial skeleton):
- Implements `generate` command (single coverage mode; sequence or chain-level).
- Normalizes coverage mode tokens (hyphen / underscore variants).
- Adds basic empty-output detection (optional non-zero exit).
- Adds chain-level MV existence guard (optional fallback or fail).
- Prepares helper utilities for forthcoming commands (`generate-all`, `plot`,
  `plot-batch`, `summary`, `validate`) not yet implemented.

NOT included intentionally at this stage (per instructions):
- Parity check integration
- Plotting / validation / summary commands
- Multi-mode orchestration (`generate-all`)
- JSON output formats

Integration:
------------
Register this sub-app in the main CLI (e.g. in `cli.py`):

    from activity_db_analysis.apps.sequence_mapping import sequence_mapping_cli as _seq_cli
    app.add_typer(_seq_cli.app, name="sequence-mapper", help="Sequence mapping reports & plots")

Exit Codes:
-----------
0  Success
1  Runtime / validation error (e.g. unsupported coverage mode, chain MV missing without fallback)
2  Empty result set (MV produced no rows) and user did not allow empty outputs

Future TODO (subsequent iterations):
- Add multi-mode generation command
- Add plotting (single + batch) commands
- Add summary & validate commands (with JSON)
- Add tests (Typer CliRunner) once full surface is implemented

Author: Auto-generated skeleton based on integration plan checklist.
"""

from __future__ import annotations

import pathlib
from typing import Dict, Iterable, Optional, Tuple

import pandas as pd
import typer
from loguru import logger

from activity_db_analysis.apps.sequence_mapping.sequence_mapping import (
    generate_target_matching_reports,
)
from activity_db_analysis.apps.sequence_mapping.sequence_mapping_parity_check import (
    run_parity_check,
    VALID_COVERAGE_MODES,
)

try:
    # Local import inside try to keep skeleton robust even when DB layer not fully configured
    from activity_db_analysis.db.connector import (
        UnifiedPostgresConnector,
        build_unified_connector,
    )
except Exception:  # pragma: no cover
    UnifiedPostgresConnector = object  # type: ignore
    build_unified_connector = None  # type: ignore


# ---------------------------------------------------------------------------
# Typer Application
# ---------------------------------------------------------------------------
app = typer.Typer(
    name="sequence-mapper",
    help="Generate sequence / chain-level target mapping reports (identity thresholds & intervals).",
    add_completion=False,
    no_args_is_help=True,
)


# ---------------------------------------------------------------------------
# Constants & Helpers
# ---------------------------------------------------------------------------

VALID_COVERAGE_MODES_HYPHEN = {"full", "no-hit", "identity-only"}
VALID_COVERAGE_MODES_UNDERSCORE = {"full", "no_hit", "identity_only"}


def normalize_coverage_mode(token: str) -> str:
    """
    Normalize user-provided coverage mode to canonical hyphenated form.

    Accepts: full, no-hit, no_hit, identity-only, identity_only (case-insensitive).

    Returns canonical: full | no-hit | identity-only
    """
    t = token.strip().lower()
    if t in VALID_COVERAGE_MODES_HYPHEN:
        return t
    if t in VALID_COVERAGE_MODES_UNDERSCORE:
        return t.replace("_", "-")
    raise typer.BadParameter(
        f"Unsupported coverage mode '{token}'. Expected one of: full | no-hit | identity-only."
    )


def coverage_mode_filename_fragment(coverage_mode: str) -> str:
    """
    Convert canonical hyphenated coverage mode to file-name-safe underscore form.
    e.g. 'identity-only' -> 'identity_only'
    """
    return coverage_mode.replace("-", "_")


def build_output_file_stems(
    output_dir: pathlib.Path,
    coverage_mode: str,
    chain_level: bool,
) -> Dict[str, pathlib.Path]:
    """
    Build expected output base paths (without sidecar TeX variants).
    Keys:
      right_bound_csv
      interval_csv
      simplified_counts_tex
      simplified_props_tex
    """
    out_dir = output_dir
    cm_fragment = coverage_mode_filename_fragment(coverage_mode)
    prefix = "chain_level_" if chain_level else ""
    stems = {
        "right_bound_csv": out_dir
        / f"{prefix}target_matches_by_sequences_right_bound_{cm_fragment}.csv",
        "interval_csv": out_dir
        / f"{prefix}target_matches_by_sequences_sequence_identity_ranges_{cm_fragment}.csv",
        "simplified_counts_tex": out_dir
        / f"{prefix}target_matching_interval_counts_{cm_fragment}{'_chain_level' if chain_level else ''}.tex",
        "simplified_props_tex": out_dir
        / f"{prefix}target_matching_interval_proportions_{cm_fragment}{'_chain_level' if chain_level else ''}.tex",
    }
    return stems


def detect_empty_interval_file(path: pathlib.Path) -> bool:
    """
    Heuristic emptiness check:
      - File missing => treat as empty
      - DataFrame loads but contains no interval rows starting with '['
    """
    if not path.exists():
        return True
    try:
        df = pd.read_csv(path, index_col=0)
        interval_like = [
            idx for idx in df.index if isinstance(idx, str) and idx.startswith("[")
        ]
        return len(interval_like) == 0
    except Exception as e:
        # If unreadable treat as non-empty to avoid false positives; caller logs error separately
        logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)
        return False


def chain_level_mv_exists(db: UnifiedPostgresConnector) -> bool:
    """
    Lightweight existence check for chain-level MV.
    """
    try:
        sql = "SELECT 1 FROM pg_matviews WHERE matviewname = 'mv_alignment_sequence_metrics_chain' LIMIT 1"
        rows = db.execute(sql) or []
        return bool(rows)
    except Exception as e:  # pragma: no cover
        logger.warning(
            f"Failed to verify chain-level MV existence (assuming missing): {e}"
        )
        return False


# ---------------------------------------------------------------------------
# Shared DB Connection Helper
# ---------------------------------------------------------------------------


def _build_connector_with_cli_overrides(
    name: Optional[str],
    user: Optional[str],
    password: Optional[str],
    host: Optional[str],
    port: Optional[int],
    schema: Optional[str],
) -> UnifiedPostgresConnector:
    """
    Build a unified connector from CLI overrides.

    If any of (db_name, user, password, host, port) are missing,
    attempt to hydrate them from the already-loaded AppConfig activity_db section. This lets users
    rely purely on --app-config without passing duplicate flags.

    Priority (per field):
      1. Explicit CLI option (non-None)
      2. AppConfig activity_db value (including password from env injection)
      3. Omit (only for schema which is optional)

    NOTE: SSH tunneling is managed at the CLI layer. Configure host=127.0.0.1 and the tunnel port
          in appconfig.toml when using the CLI-managed tunnel.

    Raises:
        RuntimeError if connector factory unavailable.
    """
    if build_unified_connector is None:  # pragma: no cover
        raise RuntimeError(
            "Database connector factory not available; ensure dependencies are installed."
        )

    # Import lazily to avoid hard dependency if siteconfig not present during certain packaging steps.
    try:
        from activity_db_analysis.siteconfig.app_config import AppConfig  # type: ignore
    except Exception:
        AppConfig = None  # type: ignore

    cfg_params = {}
    if AppConfig is not None:
        try:
            # Request parameters without overriding (empty cli_options) so AppConfig
            # returns its internally stored values (populated by root callback).
            act_params = AppConfig.set_activity_db_database_connection_parameters(
                cli_options={}
            )
            # act_params provides: name, user, password (from env), host, port
            cfg_params = {
                "db_name": getattr(act_params, "name", None),
                "user": getattr(act_params, "user", None),
                "password": getattr(act_params, "password", None),
                "host": getattr(act_params, "host", None),
                "port": getattr(act_params, "port", None),
            }
        except Exception as e:
            # Silent fallback; explicit CLI or environment-only usage will still work.
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)
            cfg_params = {}

    # Merge precedence: CLI > AppConfig > None
    merged = {
        "db_name": name if name is not None else cfg_params.get("db_name"),
        "user": user if user is not None else cfg_params.get("user"),
        "password": password if password is not None else cfg_params.get("password"),
        "host": host if host is not None else cfg_params.get("host"),
        "port": port if port is not None else cfg_params.get("port"),
        "schema": schema,  # schema not typically in AppConfig activity_db; leave as given
        "auto_connect": True,
    }

    # Remove keys still None so connector can apply its own fallbacks (if any)
    clean = {k: v for k, v in merged.items() if v is not None}
    return build_unified_connector(**clean)  # type: ignore


# ---------------------------------------------------------------------------
# Command: generate
# ---------------------------------------------------------------------------


@app.command("generate")
def cmd_generate(
    output_dir: pathlib.Path = typer.Option(
        ...,
        "--output-dir",
        help="Directory to write report artifacts (created if missing).",
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    coverage_mode: str = typer.Option(
        "full",
        "--coverage-mode",
        help="Coverage filtering mode: full | no-hit | identity-only (hyphen or underscore accepted).",
    ),
    chain_level: bool = typer.Option(
        False,
        "--chain-level",
        help="Use chain-level (intersection) MV (requires mv_alignment_sequence_metrics_chain).",
    ),
    allow_chain_fallback: bool = typer.Option(
        False,
        "--allow-chain-fallback",
        help="If chain-level requested but MV missing, fall back to sequence-level instead of exiting.",
    ),
    min_identity_floor: int = typer.Option(
        80,
        "--min-identity-floor",
        help="Reserved for future use (MV fixed at >=80). Warning if not 80.",
    ),
    allow_empty: bool = typer.Option(
        False,
        "--allow-empty",
        help="Do not treat empty MV snapshot (no interval rows) as an error (exit code 2).",
    ),
    # DB override options (optional; if omitted rely on previously loaded app config)
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="Database name override.", rich_help_panel="DB Overrides"
    ),
    user: Optional[str] = typer.Option(
        None, help="Database user override.", rich_help_panel="DB Overrides"
    ),
    password: Optional[str] = typer.Option(
        None, help="Database password override.", rich_help_panel="DB Overrides"
    ),
    host: Optional[str] = typer.Option(
        None, help="Database host override.", rich_help_panel="DB Overrides"
    ),
    port: Optional[int] = typer.Option(
        None, help="Database port override.", rich_help_panel="DB Overrides"
    ),
    schema: Optional[str] = typer.Option(
        None, help="Search path / schema override.", rich_help_panel="DB Overrides"
    ),
) -> None:
    """
    Generate sequence (or chain-level) target mapping reports for a single coverage mode.

    Artifacts written (prefix 'chain_level_' when chain-level):
      - target_matches_by_sequences_right_bound_<mode>.csv (+ .tex)
      - target_matches_by_sequences_sequence_identity_ranges_<mode>.csv (+ .tex)
      - target_matching_interval_counts_<mode>[ _chain_level].tex
      - target_matching_interval_proportions_<mode>[ _chain_level].tex

    Empty detection:
      - If the interval CSV lacks any interval rows (labels beginning with '[') it is treated as empty.
      - Exit code 2 unless --allow-empty provided.

    Chain-level guard:
      - Verifies presence of mv_alignment_sequence_metrics_chain.
      - If missing: exit (code 1) unless --allow-chain-fallback.
    """
    try:
        cov_mode = normalize_coverage_mode(coverage_mode)
    except typer.BadParameter as e:
        raise typer.Exit(code=1) from e

    if min_identity_floor != 80:
        logger.warning(
            "min_identity_floor=%s ignored (MVs fixed at >=80). Future versions may remove this option.",
            min_identity_floor,
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Sequence Mapping Report Generation")
    logger.info("Coverage mode .......... {}", cov_mode)
    logger.info("Chain-level ............ {}", "YES" if chain_level else "NO")
    logger.info("Output directory ....... {}", output_dir)

    # Establish DB connection
    try:
        db = _build_connector_with_cli_overrides(
            name=db_name,
            user=user,
            password=password,
            host=host,
            port=port,
            schema=schema,
        )
    except Exception as e:
        logger.exception("Failed to establish database connection: {}", e)
        raise typer.Exit(code=1)

    # Chain-level MV presence check
    if chain_level:
        exists = chain_level_mv_exists(db)
        if not exists:
            msg = "Chain-level MV 'mv_alignment_sequence_metrics_chain' not found."
            if allow_chain_fallback:
                logger.warning(f"{msg} Falling back to sequence-level snapshot.")
                chain_level_local = False
            else:
                logger.error(
                    f"{msg} (pass --allow-chain-fallback to continue as sequence-level)."
                )
                try:
                    db.close()
                except Exception as e:
                    logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)
                raise typer.Exit(code=1)
        else:
            chain_level_local = True
    else:
        chain_level_local = False

    # Perform generation
    try:
        generate_target_matching_reports(
            db_handler=db,
            output_file_directory=output_dir,
            chain_level=chain_level_local,
            coverage_mode=cov_mode,
            min_identity_floor=min_identity_floor,
        )
    except Exception as e:
        logger.exception("Report generation failed: {}", e)
        try:
            db.close()
        except Exception as e:
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)
        raise typer.Exit(code=1)
    finally:
        try:
            db.close()
        except Exception as e:
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)

    # Post-generation emptiness check
    stems = build_output_file_stems(output_dir, cov_mode, chain_level_local)
    interval_csv = stems["interval_csv"]
    is_empty = detect_empty_interval_file(interval_csv)
    if is_empty:
        if allow_empty:
            logger.warning(
                "Interval CSV appears empty (no interval rows) but --allow-empty set. Exit 0."
            )
        else:
            logger.error(
                "Interval CSV appears empty (no interval rows). Treating as error (use --allow-empty to override)."
            )
            raise typer.Exit(code=2)

    logger.success("Sequence mapping reports generated successfully.")


# ---------------------------------------------------------------------------
# Additional Commands: generate-all, plot, plot-batch, summary, validate
# ---------------------------------------------------------------------------


@app.command("generate-all")
def cmd_generate_all(
    output_dir: pathlib.Path = typer.Option(
        ...,
        "--output-dir",
        help="Directory to write all report artifacts.",
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    coverage_modes: list[str] = typer.Option(
        ["full", "no-hit", "identity-only"],
        "--coverage-modes",
        help="Coverage modes to generate (hyphen or underscore forms accepted).",
    ),
    include_chain_level: bool = typer.Option(
        False,
        "--include-chain-level",
        help="Also generate chain-level reports for each coverage mode.",
    ),
    allow_chain_fallback: bool = typer.Option(
        False,
        "--allow-chain-fallback",
        help="If chain-level MV missing, fall back to sequence-level for that mode.",
    ),
    allow_empty: bool = typer.Option(
        False,
        "--allow-empty",
        help="Permit empty outputs without non-zero exit.",
    ),
    fail_on_missing: bool = typer.Option(
        False,
        "--fail-on-missing",
        help="Exit with code 2 if ANY generated interval CSV is empty.",
    ),
    min_identity_floor: int = typer.Option(
        80,
        "--min-identity-floor",
        help="Retained for compatibility (ignored if !=80).",
    ),
    # DB overrides
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="Database name override.", rich_help_panel="DB Overrides"
    ),
    user: Optional[str] = typer.Option(
        None, help="Database user override.", rich_help_panel="DB Overrides"
    ),
    password: Optional[str] = typer.Option(
        None, help="Database password override.", rich_help_panel="DB Overrides"
    ),
    host: Optional[str] = typer.Option(
        None, help="Database host override.", rich_help_panel="DB Overrides"
    ),
    port: Optional[int] = typer.Option(
        None, help="Database port override.", rich_help_panel="DB Overrides"
    ),
    schema: Optional[str] = typer.Option(
        None, help="Search path / schema override.", rich_help_panel="DB Overrides"
    ),
) -> None:
    """
    Generate reports for multiple coverage modes (sequence-level and optionally chain-level).
    """
    if min_identity_floor != 80:
        logger.warning(f"min_identity_floor={min_identity_floor} ignored (fixed >=80).")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build one DB connection reused for all modes (simpler: open/close per mode if chain fallback differs)
    try:
        db = _build_connector_with_cli_overrides(
            name=db_name,
            user=user,
            password=password,
            host=host,
            port=port,
            schema=schema,
        )
    except Exception as e:
        logger.exception("Failed establishing DB connection: {}", e)
        raise typer.Exit(code=1)

    overall_empty = False
    try:
        chain_mv_available = chain_level_mv_exists(db) if include_chain_level else False
        for mode in coverage_modes:
            try:
                cov_mode = normalize_coverage_mode(mode)
            except typer.BadParameter as e:
                logger.error(str(e))
                raise typer.Exit(code=1)
            for cl_flag in [False, True] if include_chain_level else [False]:
                cl_actual = cl_flag
                if cl_flag and not chain_mv_available:
                    if allow_chain_fallback:
                        logger.warning(
                            f"Chain-level MV missing; falling back to sequence-level for coverage={cov_mode}"
                        )
                        cl_actual = False
                    else:
                        logger.error(
                            f"Chain-level MV missing for coverage={cov_mode} (use --allow-chain-fallback)."
                        )
                        raise typer.Exit(code=1)

                logger.info(
                    f"Generating coverage_mode={cov_mode} chain_level={'YES' if cl_actual else 'NO'}"
                )
                try:
                    generate_target_matching_reports(
                        db_handler=db,
                        output_file_directory=output_dir,
                        chain_level=cl_actual,
                        coverage_mode=cov_mode,
                        min_identity_floor=min_identity_floor,
                    )
                except Exception as e:
                    logger.exception(
                        f"Failed generating coverage_mode={cov_mode} chain_level={cl_actual}: {e}"
                    )
                    raise typer.Exit(code=1)

                stems = build_output_file_stems(output_dir, cov_mode, cl_actual)
                if detect_empty_interval_file(stems["interval_csv"]):
                    logger.warning(
                        f"Empty interval output for coverage_mode={cov_mode} chain_level={cl_actual}"
                    )
                    overall_empty = True
                    if not allow_empty and fail_on_missing:
                        # Record and continue to finish all modes for completeness; final exit later.
                        pass
        if overall_empty and fail_on_missing and not allow_empty:
            logger.error(
                "At least one coverage mode produced empty output. Exiting with code 2."
            )
            raise typer.Exit(code=2)
    finally:
        try:
            db.close()
        except Exception as e:
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)

    logger.success("Multi-mode generation complete.")


# ---------------------------------------------------------------------------
# Plotting Commands (wrap existing plotting function)
# ---------------------------------------------------------------------------


def _infer_coverage_mode_from_filename(path: pathlib.Path) -> Optional[str]:
    name = path.name
    for token in ("full", "no_hit", "no-hit", "identity_only", "identity-only"):
        if token in name:
            # Normalize
            token_norm = token.replace("_", "-")
            return token_norm
    return None


@app.command("plot")
def cmd_plot(
    input_csv: pathlib.Path = typer.Option(
        ...,
        "--input",
        exists=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Interval CSV produced by generate command.",
    ),
    output_file: pathlib.Path = typer.Option(
        ...,
        "--output",
        dir_okay=False,
        writable=True,
        resolve_path=True,
        help="Output figure path (.png or .pdf).",
    ),
    title: Optional[str] = typer.Option(
        None, "--title", help="Optional custom figure title."
    ),
    coverage_mode: Optional[str] = typer.Option(
        None,
        "--coverage-mode",
        help="Override / explicit coverage mode (for title inference).",
    ),
) -> None:
    """
    Create a stacked bar figure from one interval CSV.
    """
    try:
        from activity_db_analysis.apps.sequence_mapping.plot_sequence_matching_interval_stacked_bars import (
            plot_target_matching_interval_stacked_bars,
        )
    except Exception as e:
        logger.exception("Failed to import plotting function: {}", e)
        raise typer.Exit(code=1)

    if coverage_mode:
        try:
            cov_mode = normalize_coverage_mode(coverage_mode)
        except typer.BadParameter as e:
            logger.error(str(e))
            raise typer.Exit(code=1)
    else:
        inferred = _infer_coverage_mode_from_filename(input_csv)
        cov_mode = inferred if inferred else "full"

    if title is None:
        title = f"Unique PDB–ChEMBL sequence pairs by identity intervals"

    output_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        plot_target_matching_interval_stacked_bars(
            csv_path=input_csv,
            out_path=output_file,
            figure_title=title,
        )
    except Exception as e:
        logger.exception("Plotting failed: {}", e)
        raise typer.Exit(code=1)
    logger.success("Wrote figure: %s", output_file)


@app.command("plot-batch")
def cmd_plot_batch(
    input_dir: pathlib.Path = typer.Option(
        ..., "--input-dir", exists=True, file_okay=False, resolve_path=True
    ),
    output_dir: pathlib.Path = typer.Option(
        ..., "--output-dir", file_okay=False, resolve_path=True
    ),
    coverage_modes: list[str] = typer.Option(
        ["full", "no-hit", "identity-only"],
        "--coverage-modes",
        help="Subset of coverage modes to plot.",
    ),
    include_chain_level: bool = typer.Option(
        True,
        "--include-chain-level/--no-chain-level",
        help="Include chain-level interval CSVs when present.",
    ),
    pdf: bool = typer.Option(False, "--pdf", help="Also produce PDF alongside PNG."),
    overwrite: bool = typer.Option(
        True, "--overwrite/--no-overwrite", help="Overwrite existing figures."
    ),
) -> None:
    """
    Batch plot all (sequence + optionally chain-level) interval CSVs for selected coverage modes.
    """
    try:
        from activity_db_analysis.apps.sequence_mapping.plot_sequence_matching_interval_stacked_bars import (
            plot_target_matching_interval_stacked_bars,
        )
    except Exception as e:
        logger.exception("Failed to import plotting function: {}", e)
        raise typer.Exit(code=1)

    from typing import List

    output_dir.mkdir(parents=True, exist_ok=True)
    requested = [normalize_coverage_mode(c) for c in coverage_modes]

    # Added detailed context logging
    logger.info(f"Plot-batch input directory: {input_dir}")
    logger.info(f"Plot-batch output directory: {output_dir}")
    logger.info(f"Requested coverage modes: {', '.join(requested)}")
    logger.info(
        f"Include chain-level: {include_chain_level} | PDF: {pdf} | Overwrite: {overwrite}"
    )

    # Tracking
    produced_sets = 0
    created: List[pathlib.Path] = []
    skipped_existing: List[pathlib.Path] = []
    missing_inputs: List[pathlib.Path] = []
    failed: List[tuple[pathlib.Path, str]] = []

    for cov_mode in requested:
        fragment = coverage_mode_filename_fragment(cov_mode)
        seq_csv = (
            input_dir
            / f"target_matches_by_sequences_sequence_identity_ranges_{fragment}.csv"
        )
        for is_chain, csv_path in (
            (False, seq_csv),
            (
                True,
                input_dir
                / f"chain_level_target_matches_by_sequences_sequence_identity_ranges_{fragment}.csv",
            ),
        ):
            if is_chain and not include_chain_level:
                logger.debug(
                    f"Skipping chain-level for {cov_mode} (chain-level disabled)"
                )
                continue
            if not csv_path.exists():
                logger.warning(f"Missing interval CSV (skipping): {csv_path}")
                missing_inputs.append(csv_path)
                continue

            suffix_chain = "_chain_level" if is_chain else ""
            png_out = (
                output_dir
                / f"target_matching_interval_stacked_bars{suffix_chain}_{fragment}.png"
            )
            pdf_out = (
                output_dir
                / f"target_matching_interval_stacked_bars{suffix_chain}_{fragment}.pdf"
            )

            if png_out.exists() and not overwrite:
                logger.info(f"Skipping existing figure (no-overwrite): {png_out}")
                skipped_existing.append(png_out)
                if pdf and pdf_out.exists():
                    skipped_existing.append(pdf_out)
                continue

            title = f"{'Chain-level unique' if is_chain else 'Unique'} PDB–ChEMBL sequence pairs"
            try:
                plot_target_matching_interval_stacked_bars(
                    csv_path=csv_path,
                    out_path=png_out,
                    figure_title=title,
                    additional_output_paths=[pdf_out] if pdf else None,
                )
                produced_sets += 1
                created.append(png_out)
                if pdf:
                    created.append(pdf_out)
                logger.success(
                    f"Created figure set ({png_out.name}{' + PDF' if pdf else ''}) from {csv_path.name}"
                )
            except Exception as e:
                failed.append((csv_path, str(e)))
                logger.error(f"Failed plotting {csv_path}: {e}")

    # Summary
    if missing_inputs:
        logger.warning(
            f"Missing {len(missing_inputs)} required interval CSV(s). First few: {', '.join(p.name for p in missing_inputs[:5])}"
        )
        # (converted to f-string above)
        # (converted to f-string above))
    if skipped_existing:
        logger.info(
            f"Skipped {len(skipped_existing)} existing figure file(s) due to --no-overwrite or default keep: {', '.join(p.name for p in skipped_existing[:5])}"
        )
        # (converted to f-string above)
        # (converted to f-string above))
    if failed:
        logger.error(
            f"Failed to produce {len(failed)} figure set(s). Examples: {', '.join(f[0].name for f in failed[:5])}"
        )
        # (converted to f-string above)
        # (converted to f-string above))

    if produced_sets == 0:
        logger.warning("No new figure sets produced.")
    else:
        logger.success(f"Produced {produced_sets} new figure set(s).")


# ---------------------------------------------------------------------------
# Summary & Validation Helpers
# ---------------------------------------------------------------------------


def _list_interval_csvs(directory: pathlib.Path) -> Dict[str, Dict[str, pathlib.Path]]:
    """
    Return mapping:
      {
        coverage_mode (hyphen): {
           'sequence': Path | None,
           'chain': Path | None
        },
        ...
      }
    """
    mapping: Dict[str, Dict[str, pathlib.Path]] = {}
    for p in directory.glob(
        "target_matches_by_sequences_sequence_identity_ranges_*.csv"
    ):
        name = p.name
        # Extract coverage fragment after last underscore group
        frag = name.replace(
            "target_matches_by_sequences_sequence_identity_ranges_", ""
        ).replace(".csv", "")
        is_chain = False
        cov = frag
        if name.startswith("chain_level_"):
            continue  # guard: chain-level prefixed naming uses chain_level_* pattern handled below
        # Already sequence-level
        cov_mode = frag
        cov_mode = cov_mode.replace("_", "-")
        mapping.setdefault(cov_mode, {})["sequence"] = p
    for p in directory.glob(
        "chain_level_target_matches_by_sequences_sequence_identity_ranges_*.csv"
    ):
        frag = p.name.replace(
            "chain_level_target_matches_by_sequences_sequence_identity_ranges_", ""
        ).replace(".csv", "")
        cov_mode = frag.replace("_", "-")
        mapping.setdefault(cov_mode, {})["chain"] = p
    return mapping


@app.command("summary")
def cmd_summary(
    directory: pathlib.Path = typer.Option(
        ..., "--directory", exists=True, file_okay=False, resolve_path=True
    ),
) -> None:
    """
    Summarize available interval + right-bound + simplified artifacts in a directory.
    """
    mapping = _list_interval_csvs(directory)
    if not mapping:
        logger.warning(f"No interval CSVs found in {directory}")
        raise typer.Exit(code=0)

    for cov_mode, sides in sorted(mapping.items()):
        fragment = coverage_mode_filename_fragment(cov_mode)
        rb_seq = directory / f"target_matches_by_sequences_right_bound_{fragment}.csv"
        rb_chain = (
            directory
            / f"chain_level_target_matches_by_sequences_right_bound_{fragment}.csv"
        )
        simp_counts_seq = directory / f"target_matching_interval_counts_{fragment}.tex"
        simp_props_seq = (
            directory / f"target_matching_interval_proportions_{fragment}.tex"
        )
        simp_counts_chain = (
            directory / f"target_matching_interval_counts_{fragment}_chain_level.tex"
        )
        simp_props_chain = (
            directory
            / f"target_matching_interval_proportions_{fragment}_chain_level.tex"
        )

        logger.info(f"Coverage: {cov_mode}")
        seq_csv = sides.get("sequence")
        chain_csv = sides.get("chain")
        logger.info(
            f"  Sequence interval CSV ......... {seq_csv if seq_csv and seq_csv.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Chain interval CSV ............ {chain_csv if chain_csv and chain_csv.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Right-bound (sequence) ........ {'OK' if rb_seq.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Right-bound (chain) ........... {'OK' if rb_chain.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Simplified counts (seq) ....... {'OK' if simp_counts_seq.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Simplified props (seq) ........ {'OK' if simp_props_seq.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Simplified counts (chain) ...... {'OK' if simp_counts_chain.exists() else 'MISSING'}"
        )
        logger.info(
            f"  Simplified props (chain) ....... {'OK' if simp_props_chain.exists() else 'MISSING'}"
        )


def _load_csv_indexed(path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def _validate_interval_vs_right_bound(
    interval_df: pd.DataFrame, right_df: pd.DataFrame
) -> list[str]:
    issues: list[str] = []
    # right_df index: ">=100", ">=95", ... ; interval_df index: "[100]", "[95,100)", ...
    if interval_df.empty or right_df.empty:
        return issues
    # For each side/seqtype RAW columns: COMPONENT_ATOMSEQ_RAW etc (interval) vs >=80 row right-bound
    if ">=80" not in right_df.index:
        issues.append("Right-bound table missing '>=80' row.")
        return issues
    right_row = right_df.loc[">=80"]
    interval_rows = [r for r in interval_df.index if r.startswith("[")]
    sum_row = interval_df.loc[interval_rows].sum(numeric_only=True)
    # Columns to compare: *_RAW columns only
    for col in [c for c in sum_row.index if c.endswith("_RAW")]:
        if col not in right_row.index:
            # right-bound RAW columns have no COMPONENT_ prefix until after we added; in combined we prefixed earlier.
            continue
        lhs = int(sum_row[col] or 0)
        rhs = int(right_row[col] or 0)
        if lhs != rhs:
            issues.append(
                f"Interval sum mismatch for {col}: intervals={lhs} vs >=80={rhs}"
            )
    return issues


def _validate_monotonic(right_df: pd.DataFrame) -> list[str]:
    issues: list[str] = []
    order = [">=100", ">=95", ">=90", ">=85", ">=80"]
    raw_cols = [c for c in right_df.columns if c.endswith("_RAW")]
    for col in raw_cols:
        last_val = None
        for thr in order:
            if thr not in right_df.index:
                continue
            v = int(right_df.loc[thr, col] or 0)
            if last_val is not None and v < last_val:
                issues.append(
                    f"Monotonicity violation {col} at {thr}: {v} < {last_val}"
                )
            last_val = v
    return issues


def _validate_chain_not_exceed(
    sequence_interval: Optional[pd.DataFrame], chain_interval: Optional[pd.DataFrame]
) -> list[str]:
    issues: list[str] = []
    if (
        sequence_interval is None
        or chain_interval is None
        or sequence_interval.empty
        or chain_interval.empty
    ):
        return issues
    interval_rows = [r for r in chain_interval.index if r.startswith("[")]
    for col in [c for c in chain_interval.columns if c.endswith("_RAW")]:
        chain_sum = int(chain_interval.loc[interval_rows, col].sum())
        if col in sequence_interval.columns:
            seq_sum = int(
                sequence_interval.loc[
                    [r for r in sequence_interval.index if r.startswith("[")], col
                ].sum()
            )
            if chain_sum > seq_sum:
                issues.append(
                    f"Chain-level exceeds sequence-level total for {col}: chain={chain_sum} seq={seq_sum}"
                )
    return issues


@app.command("validate")
def cmd_validate(
    directory: pathlib.Path = typer.Option(
        ..., "--directory", exists=True, file_okay=False, resolve_path=True
    ),
    coverage_modes: Optional[list[str]] = typer.Option(
        None,
        "--coverage-modes",
        help="Limit validation to these coverage modes (hyphen/underscore accepted). If omitted, auto-discover.",
    ),
    include_chain_checks: bool = typer.Option(
        True,
        "--include-chain-checks/--no-chain-checks",
        help="Validate chain ≤ sequence totals where both exist.",
    ),
) -> None:
    """
    Validate internal consistency of generated artifacts.
    """
    mapping = _list_interval_csvs(directory)
    if not mapping:
        logger.error(f"No interval CSVs found in {directory}")
        raise typer.Exit(code=1)

    selected_covs = set()
    if coverage_modes:
        for m in coverage_modes:
            selected_covs.add(normalize_coverage_mode(m))
    else:
        selected_covs = {normalize_coverage_mode(k) for k in mapping.keys()}

    failures: list[str] = []
    for cov_mode in sorted(selected_covs):
        fragment = coverage_mode_filename_fragment(cov_mode)
        seq_interval = mapping.get(cov_mode, {}).get("sequence")
        chain_interval = mapping.get(cov_mode, {}).get("chain")
        seq_interval_df = (
            _load_csv_indexed(seq_interval)
            if seq_interval and seq_interval.exists()
            else None
        )
        chain_interval_df = (
            _load_csv_indexed(chain_interval)
            if chain_interval and chain_interval.exists()
            else None
        )

        # Right-bound files
        seq_right = (
            directory / f"target_matches_by_sequences_right_bound_{fragment}.csv"
        )
        chain_right = (
            directory
            / f"chain_level_target_matches_by_sequences_right_bound_{fragment}.csv"
        )
        seq_right_df = _load_csv_indexed(seq_right) if seq_right.exists() else None
        chain_right_df = (
            _load_csv_indexed(chain_right) if chain_right.exists() else None
        )

        logger.info(f"Validating coverage={cov_mode}")
        if seq_interval_df is not None and seq_right_df is not None:
            failures.extend(
                [
                    f"[{cov_mode}] {msg}"
                    for msg in _validate_interval_vs_right_bound(
                        seq_interval_df, seq_right_df
                    )
                ]
            )
            failures.extend(
                [f"[{cov_mode}] {msg}" for msg in _validate_monotonic(seq_right_df)]
            )
        if chain_interval_df is not None and chain_right_df is not None:
            failures.extend(
                [
                    f"[{cov_mode} chain] {msg}"
                    for msg in _validate_interval_vs_right_bound(
                        chain_interval_df, chain_right_df
                    )
                ]
            )
            failures.extend(
                [
                    f"[{cov_mode} chain] {msg}"
                    for msg in _validate_monotonic(chain_right_df)
                ]
            )
        if include_chain_checks:
            failures.extend(
                [
                    f"[{cov_mode}] {msg}"
                    for msg in _validate_chain_not_exceed(
                        seq_interval_df, chain_interval_df
                    )
                ]
            )

    if failures:
        logger.error(f"Validation FAILED with {len(failures)} issue(s):")
        for f in failures:
            logger.error(f"  {f}")
        raise typer.Exit(code=1)

    logger.success("Validation PASSED.")


# ---------------------------------------------------------------------------
# Module Exports
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Parity Check Command
# ---------------------------------------------------------------------------


@app.command("parity-check")
def cmd_parity_check(
    coverage_modes: list[str] = typer.Option(
        ["identity-only", "no-hit", "full"],
        "--coverage-modes",
        help=f"Coverage modes to validate (subset of {VALID_COVERAGE_MODES}).",
    ),
    allow_differences: bool = typer.Option(
        False,
        "--allow-differences",
        help="Do not exit non-zero on differences (informational run).",
    ),
    min_identity_floor: int = typer.Option(
        80,
        "--min-identity-floor",
        help="Minimum identity floor for interval bin inclusion.",
    ),
    # DB override options
    db_name: Optional[str] = typer.Option(
        None, "--db-name", help="Database name override.", rich_help_panel="DB Overrides"
    ),
    user: Optional[str] = typer.Option(
        None, help="Database user override.", rich_help_panel="DB Overrides"
    ),
    password: Optional[str] = typer.Option(
        None, help="Database password override.", rich_help_panel="DB Overrides"
    ),
    host: Optional[str] = typer.Option(
        None, help="Database host override.", rich_help_panel="DB Overrides"
    ),
    port: Optional[int] = typer.Option(
        None, help="Database port override.", rich_help_panel="DB Overrides"
    ),
    schema: Optional[str] = typer.Option(
        None, help="Search path / schema override.", rich_help_panel="DB Overrides"
    ),
) -> None:
    """
    Run parity check comparing legacy (inline) vs unified tier-view logic.

    Validates both right-bound thresholds (>=100,95,90,85,80) and interval
    buckets ([100], [95,100), etc.) across dimensions: chembl_type, pdb_stid,
    pair_count_raw, distinct_sequences.

    Exit Codes:
      0 = Perfect parity (or differences allowed)
      1 = Differences detected (and not allowed)
      2 = Runtime / connection / SQL error
    """
    try:
        db = _build_connector_with_cli_overrides(
            name=db_name,
            user=user,
            password=password,
            host=host,
            port=port,
            schema=schema,
        )
    except Exception as e:
        logger.error(f"Failed to establish database connection: {e}")
        raise typer.Exit(code=2)

    try:
        rc = run_parity_check(
            db=db,
            coverage_modes=coverage_modes,
            allow_differences=allow_differences,
            min_identity_floor=min_identity_floor,
        )
        raise typer.Exit(code=rc)
    except typer.Exit:
        raise
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        raise typer.Exit(code=2)
    except Exception as e:
        logger.exception(f"Unexpected error during parity check: {e}")
        raise typer.Exit(code=2)
    finally:
        try:
            db.close()
        except Exception as e:
            logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)


# ---------------------------------------------------------------------------
# Module Exports
# ---------------------------------------------------------------------------

__all__ = [
    "app",
    "normalize_coverage_mode",
    "coverage_mode_filename_fragment",
    "build_output_file_stems",
    "cmd_generate",
    "cmd_generate_all",
    "cmd_plot",
    "cmd_plot_batch",
    "cmd_summary",
    "cmd_validate",
    "cmd_parity_check",
]
