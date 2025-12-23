"""
CLI Tunnel Wrapper for database commands.

This module provides integration between the SSH tunnel manager and typer CLI
commands, ensuring database connectivity is available when needed.

Key Features
------------
- Typer callback that checks/establishes tunnel before database commands
- Decorator for wrapping individual CLI commands with tunnel management
- Integration with existing AppConfig for configuration
- Clear error messages when tunnel cannot be established

Usage
-----
```python
import typer
from activity_db_analysis.db.cli_tunnel_wrapper import (
    with_tunnel,
    tunnel_callback,
    TunnelMode,
)

app = typer.Typer()

# Option 1: Use as a typer callback for all commands in an app
@app.callback()
def main(ctx: typer.Context):
    tunnel_callback(ctx)

# Option 2: Wrap individual commands
@app.command()
@with_tunnel
def my_command():
    # Tunnel is guaranteed active
    ...

# Option 3: Explicit control with mode selection
@app.command()
@with_tunnel(mode=TunnelMode.REQUIRED)
def strict_command():
    # Will fail if tunnel cannot be established
    ...
```
"""

from __future__ import annotations

import functools
from enum import Enum
from typing import Callable, Optional, ParamSpec, TypeVar

import typer
from loguru import logger

from activity_db_analysis.db.ssh_tunnel import (
    SSHTunnelManager,
    TunnelConfig,
    TunnelError,
    TunnelStartError,
    get_tunnel_status,
)

# Type variables for decorators
P = ParamSpec("P")
R = TypeVar("R")


class TunnelMode(str, Enum):
    """Tunnel handling modes for CLI commands."""

    AUTO = "auto"
    """Automatically start tunnel if not running, stop when done if started."""

    REQUIRED = "required"
    """Require tunnel to be active; fail if cannot establish."""

    OPTIONAL = "optional"
    """Try to establish tunnel but continue without if it fails."""

    SKIP = "skip"
    """Skip tunnel management entirely (for local development)."""


def with_tunnel(
    func: Optional[Callable[P, R]] = None,
    *,
    mode: TunnelMode = TunnelMode.AUTO,
    config: Optional[TunnelConfig] = None,
) -> Callable[P, R]:
    """
    Decorator that wraps a CLI command with SSH tunnel management.

    Parameters
    ----------
    func : Callable, optional
        The function to wrap. If None, returns a partial decorator.
    mode : TunnelMode
        How to handle tunnel management:
        - AUTO: Start if needed, stop when done (default)
        - REQUIRED: Fail if tunnel cannot be established
        - OPTIONAL: Try to start but continue without if fails
        - SKIP: Don't manage tunnel at all
    config : TunnelConfig, optional
        Custom tunnel configuration. Uses defaults if not provided.

    Returns
    -------
    Callable
        Wrapped function with tunnel management.

    Examples
    --------
    ```python
    @with_tunnel
    def fetch_data():
        # Tunnel is active here
        ...

    @with_tunnel(mode=TunnelMode.REQUIRED)
    def critical_fetch():
        # Will raise if tunnel fails
        ...
    ```
    """

    def decorator(fn: Callable[P, R]) -> Callable[P, R]:
        @functools.wraps(fn)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            if mode == TunnelMode.SKIP:
                logger.debug("Tunnel management skipped for {}", fn.__name__)
                return fn(*args, **kwargs)

            tunnel_manager = SSHTunnelManager(config)

            try:
                if mode == TunnelMode.OPTIONAL:
                    try:
                        tunnel_manager.start()
                    except TunnelError as e:
                        logger.warning(
                            "Could not establish SSH tunnel (continuing anyway): {}",
                            e,
                        )
                else:
                    # AUTO or REQUIRED
                    tunnel_manager.start()

                return fn(*args, **kwargs)

            except TunnelStartError as e:
                if mode == TunnelMode.REQUIRED:
                    logger.exception("SSH tunnel required but failed to start: {}", e)
                    raise typer.Exit(code=1) from e
                raise

            finally:
                # Only stop if we're in AUTO mode and we own the tunnel
                if mode == TunnelMode.AUTO:
                    tunnel_manager.stop()

        return wrapper

    # Handle both @with_tunnel and @with_tunnel() syntax
    if func is not None:
        return decorator(func)
    return decorator


def tunnel_callback(
    ctx: typer.Context,
    manage_tunnel: bool = True,
    tunnel_mode: TunnelMode = TunnelMode.AUTO,
) -> None:
    """
    Typer callback that manages tunnel state for all commands in an app.

    This should be used as or called from a typer app callback to ensure
    tunnel is available for all database-accessing commands.

    Parameters
    ----------
    ctx : typer.Context
        Typer context object.
    manage_tunnel : bool
        Whether to manage the tunnel. Set False to disable.
    tunnel_mode : TunnelMode
        How to handle tunnel management.

    Examples
    --------
    ```python
    app = typer.Typer()

    @app.callback()
    def main(ctx: typer.Context):
        tunnel_callback(ctx)

    @app.command()
    def my_command():
        # Tunnel is available
        ...
    ```
    """
    if not manage_tunnel or tunnel_mode == TunnelMode.SKIP:
        return

    # Store tunnel manager in context for cleanup
    tunnel_manager = SSHTunnelManager()

    try:
        if tunnel_mode == TunnelMode.OPTIONAL:
            try:
                tunnel_manager.start()
                ctx.obj = ctx.obj or {}
                ctx.obj["tunnel_manager"] = tunnel_manager
            except TunnelError as e:
                logger.warning("Could not establish SSH tunnel: {}", e)
        else:
            tunnel_manager.start()
            ctx.obj = ctx.obj or {}
            ctx.obj["tunnel_manager"] = tunnel_manager

    except TunnelStartError as e:
        if tunnel_mode == TunnelMode.REQUIRED:
            logger.exception("SSH tunnel required but failed: {}", e)
            raise typer.Exit(code=1) from e
        raise


def add_tunnel_option(
    manage_tunnel: bool = typer.Option(
        True,
        "--manage-tunnel/--no-tunnel",
        help="Automatically manage SSH tunnel for database access.",
    ),
    tunnel_host: Optional[str] = typer.Option(
        None,
        "--tunnel-host",
        help="Override SSH host alias for tunnel (default: from config).",
    ),
    tunnel_local_port: Optional[int] = typer.Option(
        None,
        "--tunnel-local-port",
        help="Override local port for tunnel (default: 5433).",
    ),
) -> tuple[bool, Optional[str], Optional[int]]:
    """
    Standard tunnel options to add to CLI commands.

    Use this as a dependency in typer commands to provide consistent
    tunnel configuration options.

    Returns
    -------
    tuple
        (manage_tunnel, tunnel_host, tunnel_local_port)

    Examples
    --------
    ```python
    @app.command()
    def my_command(
        tunnel_opts: tuple = typer.Depends(add_tunnel_option),
    ):
        manage_tunnel, host, port = tunnel_opts
        ...
    ```
    """
    return (manage_tunnel, tunnel_host, tunnel_local_port)


def print_tunnel_status() -> None:
    """Print current tunnel status to console."""
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


def create_tunnel_aware_callback(
    original_callback: Optional[Callable] = None,
    mode: TunnelMode = TunnelMode.AUTO,
) -> Callable:
    """
    Create a typer callback that wraps an existing callback with tunnel management.

    Parameters
    ----------
    original_callback : Callable, optional
        Existing callback to wrap. If None, creates a new callback.
    mode : TunnelMode
        Tunnel management mode.

    Returns
    -------
    Callable
        New callback with tunnel management.

    Examples
    --------
    ```python
    @app.callback()
    @create_tunnel_aware_callback
    def main(ctx: typer.Context):
        # Original callback logic
        ...
    ```
    """

    def callback_wrapper(
        ctx: typer.Context,
        *args,
        manage_tunnel: bool = typer.Option(
            True,
            "--manage-tunnel/--no-tunnel",
            help="Automatically manage SSH tunnel.",
        ),
        **kwargs,
    ):
        # Manage tunnel
        if manage_tunnel:
            tunnel_callback(ctx, manage_tunnel=True, tunnel_mode=mode)

        # Call original callback if provided
        if original_callback is not None:
            return original_callback(ctx, *args, **kwargs)

    if original_callback is not None:
        # Preserve function metadata
        callback_wrapper.__name__ = original_callback.__name__
        callback_wrapper.__doc__ = original_callback.__doc__

    return callback_wrapper


class TunnelAwareCommand:
    """
    Mixin or base class for commands that require database access.

    Provides standardized tunnel management for command classes.

    Examples
    --------
    ```python
    class MyDatabaseCommand(TunnelAwareCommand):
        def execute(self):
            with self.tunnel():
                # Database access here
                ...
    ```
    """

    def __init__(
        self,
        tunnel_mode: TunnelMode = TunnelMode.AUTO,
        tunnel_config: Optional[TunnelConfig] = None,
    ):
        self._tunnel_mode = tunnel_mode
        self._tunnel_config = tunnel_config
        self._tunnel_manager: Optional[SSHTunnelManager] = None

    def tunnel(self) -> SSHTunnelManager:
        """
        Get a tunnel manager configured for this command.

        Returns a context manager that can be used with `with` statement.
        """
        if self._tunnel_manager is None:
            self._tunnel_manager = SSHTunnelManager(self._tunnel_config)
        return self._tunnel_manager

    def ensure_tunnel(self) -> None:
        """Ensure tunnel is active, starting it if necessary."""
        manager = self.tunnel()
        if self._tunnel_mode != TunnelMode.SKIP:
            manager.start()

    def release_tunnel(self) -> None:
        """Release tunnel if we own it."""
        if self._tunnel_manager is not None:
            self._tunnel_manager.stop()
