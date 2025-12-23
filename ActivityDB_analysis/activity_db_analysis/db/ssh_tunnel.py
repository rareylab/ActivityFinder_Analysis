"""
SSH Tunnel Manager using native SSH subprocess.

This module provides a robust SSH tunnel implementation that leverages the native
`ssh` command via subprocess, enabling full support for ~/.ssh/config features
including ProxyCommand, ProxyJump, and complex multi-hop configurations.

Key Features
------------
- Uses native `ssh -f -N` for backgrounded tunnel with LocalForward
- Context manager pattern for automatic setup/teardown
- Checks for existing tunnels to avoid conflicts
- Supports configurable SSH host alias from ~/.ssh/config
- Thread-safe singleton pattern for shared tunnel access
- Graceful cleanup on exit

Usage
-----
```python
from activity_db_analysis.db.ssh_tunnel import SSHTunnelManager, ensure_tunnel

# Context manager pattern (recommended for scripts)
with SSHTunnelManager() as tunnel:
    # tunnel is active, connect to localhost:<local_port>
    conn = psycopg2.connect(host="127.0.0.1", port=5433, ...)

# Decorator pattern for CLI commands
@ensure_tunnel
def my_database_command():
    # tunnel is guaranteed to be active
    ...

# Manual control
tunnel = SSHTunnelManager()
tunnel.start()
try:
    # do work
finally:
    tunnel.stop()
```
"""

from __future__ import annotations

import atexit
import os
import signal
import socket
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import Callable, Optional, ParamSpec, TypeVar

from loguru import logger

# Type variables for decorator
P = ParamSpec("P")
R = TypeVar("R")


@dataclass
class TunnelConfig:
    """Configuration for SSH tunnel."""

    ssh_host_alias: str = "db_tunnel"  # Configure in ~/.ssh/config
    local_host: str = "127.0.0.1"
    local_port: int = 5433
    remote_host: str = "database.example.com"  # Set via appconfig.toml
    remote_port: int = 5432
    ssh_config_path: Path = field(
        default_factory=lambda: Path("~/.ssh/config").expanduser()
    )
    connection_timeout: float = 30.0
    health_check_interval: float = 5.0


class TunnelError(Exception):
    """Base exception for tunnel-related errors."""

    pass


class TunnelStartError(TunnelError):
    """Raised when tunnel fails to start."""

    pass


class TunnelHealthCheckError(TunnelError):
    """Raised when tunnel health check fails."""

    pass


class SSHTunnelManager:
    """
    Manages an SSH tunnel using native ssh subprocess.

    This implementation uses the native `ssh` command to establish tunnels,
    which provides full support for ~/.ssh/config features like ProxyCommand,
    ProxyJump, IdentityFile, and other advanced configurations.

    The manager tracks whether it created the tunnel or if one was already
    running, and only tears down tunnels it created.
    """

    _instance: Optional["SSHTunnelManager"] = None
    _lock: Lock = Lock()

    def __new__(cls, config: Optional[TunnelConfig] = None) -> "SSHTunnelManager":
        """Singleton pattern to ensure only one tunnel manager exists."""
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
                cls._instance._initialized = False
            return cls._instance

    def __init__(self, config: Optional[TunnelConfig] = None) -> None:
        if self._initialized:
            return
        self._initialized = True
        self._config = config or TunnelConfig()
        self._process: Optional[subprocess.Popen] = None
        self._owned_tunnel: bool = False
        self._started: bool = False

        # Register cleanup on exit
        atexit.register(self._cleanup_on_exit)

    @property
    def config(self) -> TunnelConfig:
        """Return the tunnel configuration."""
        return self._config

    @property
    def is_tunnel_active(self) -> bool:
        """Check if a tunnel is currently active on the configured port."""
        return self._is_port_listening(self._config.local_host, self._config.local_port)

    @property
    def owns_tunnel(self) -> bool:
        """Return True if this manager created the currently active tunnel."""
        return self._owned_tunnel and self._process is not None

    def start(self) -> "SSHTunnelManager":
        """
        Start the SSH tunnel if not already running.

        Returns self for method chaining.

        Raises
        ------
        TunnelStartError
            If tunnel fails to start within the timeout period.
        """
        if self._started:
            logger.debug("Tunnel manager already started, skipping")
            return self

        # Check if tunnel is already running (externally managed)
        if self.is_tunnel_active:
            logger.info(
                "SSH tunnel already active on {}:{} (externally managed)",
                self._config.local_host,
                self._config.local_port,
            )
            self._owned_tunnel = False
            self._started = True
            return self

        # Start new tunnel
        logger.info(
            "Starting SSH tunnel: {} -> {}:{} via {}",
            f"{self._config.local_host}:{self._config.local_port}",
            self._config.remote_host,
            self._config.remote_port,
            self._config.ssh_host_alias,
        )

        self._start_ssh_tunnel()
        self._owned_tunnel = True
        self._started = True

        return self

    def stop(self) -> None:
        """
        Stop the SSH tunnel if this manager owns it.

        If the tunnel was already running when start() was called (externally
        managed), this method does nothing.
        """
        if not self._started:
            return

        if not self._owned_tunnel:
            logger.debug(
                "Tunnel was externally managed, not stopping (port {})",
                self._config.local_port,
            )
            self._started = False
            return

        if self._process is not None:
            logger.info("Stopping SSH tunnel (PID: {})", self._process.pid)
            self._terminate_process()
            self._process = None

        self._owned_tunnel = False
        self._started = False

    def health_check(self) -> bool:
        """
        Perform a health check on the tunnel.

        Returns True if the tunnel port is listening and accepting connections.
        """
        return self._is_port_listening(self._config.local_host, self._config.local_port)

    def __enter__(self) -> "SSHTunnelManager":
        """Context manager entry - start tunnel."""
        return self.start()

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - stop tunnel if owned."""
        self.stop()

    def _start_ssh_tunnel(self) -> None:
        """Start the SSH tunnel subprocess."""
        cfg = self._config

        # Verify SSH config exists
        if not cfg.ssh_config_path.exists():
            raise TunnelStartError(f"SSH config file not found: {cfg.ssh_config_path}")

        # Build SSH command
        # -N: don't execute remote command (tunnel only)
        # -o ExitOnForwardFailure=yes: exit if port forwarding fails
        # -o ServerAliveInterval=60: send keepalive every 60s
        # -o ServerAliveCountMax=3: disconnect after 3 missed keepalives
        # -o BatchMode=yes: don't prompt for passwords, rely on ssh-agent
        # Note: We do NOT use -f because we want to manage the process ourselves
        cmd = [
            "ssh",
            "-N",
            "-o",
            "ExitOnForwardFailure=yes",
            "-o",
            "ServerAliveInterval=60",
            "-o",
            "ServerAliveCountMax=3",
            "-o",
            "BatchMode=yes",
            "-L",
            f"{cfg.local_host}:{cfg.local_port}:{cfg.remote_host}:{cfg.remote_port}",
            cfg.ssh_host_alias,
        ]

        logger.debug("Executing SSH command: {}", " ".join(cmd))

        try:
            # Start SSH process in background (no -f flag, we manage it ourselves)
            # BatchMode=yes means it will fail fast if key isn't in agent
            self._process = subprocess.Popen(
                cmd,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                start_new_session=True,  # Detach from terminal
            )

            # Give SSH a moment to start and potentially fail
            time.sleep(0.5)

            # Check if process died immediately (e.g., auth failure)
            poll_result = self._process.poll()
            if poll_result is not None:
                # Process exited - get error message
                stderr_output = ""
                if self._process.stderr:
                    try:
                        stderr_output = self._process.stderr.read().decode().strip()
                    except Exception as e:
                        logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)

                if "Permission denied" in stderr_output or poll_result == 255:
                    # Use env var for SSH key name (allows users to configure their key)
                    ssh_key_name = os.getenv("ACTFIND_SSH_KEY_NAME", "id_ed25519")
                    raise TunnelStartError(
                        f"SSH authentication failed. Ensure your SSH key is loaded in ssh-agent. "
                        f"Run: ssh-add ~/.ssh/{ssh_key_name}\n"
                        f"stderr: {stderr_output}"
                    )
                raise TunnelStartError(
                    f"SSH tunnel failed to start (exit code {poll_result}): {stderr_output}"
                )

        except FileNotFoundError:
            raise TunnelStartError("ssh command not found in PATH")

        # Wait for tunnel to become active
        if not self._wait_for_port(cfg.local_host, cfg.local_port, timeout=15.0):
            # Process may have failed after initial check
            if self._process and self._process.poll() is not None:
                stderr_output = ""
                if self._process.stderr:
                    try:
                        stderr_output = self._process.stderr.read().decode().strip()
                    except Exception as e:
                        logger.debug("Cleanup/close operation failed (safe to ignore): {}", e)
                raise TunnelStartError(
                    f"SSH tunnel process exited unexpectedly: {stderr_output}"
                )
            raise TunnelStartError(
                f"Tunnel port {cfg.local_host}:{cfg.local_port} not listening after start. "
                "The SSH connection may have succeeded but port forwarding failed."
            )

        # Wrap process in PIDWrapper for consistent interface
        if self._process:
            pid = self._process.pid
            self._process = _PIDWrapper(pid)
            logger.info(
                "SSH tunnel established (PID: {}, local: {}:{})",
                pid,
                cfg.local_host,
                cfg.local_port,
            )

    def _find_tunnel_process(self) -> Optional[subprocess.Popen]:
        """
        Find the SSH tunnel process by checking listening processes on the port.

        Returns a Popen-like object wrapping the PID, or None if not found.
        """
        cfg = self._config

        try:
            # Use lsof to find the process listening on our port
            result = subprocess.run(
                ["lsof", "-ti", f"TCP:{cfg.local_port}", "-sTCP:LISTEN"],
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and result.stdout.strip():
                pids = result.stdout.strip().split("\n")
                if pids:
                    pid = int(pids[0])
                    # Create a minimal object to track the PID
                    return _PIDWrapper(pid)

        except (subprocess.SubprocessError, ValueError) as e:
            logger.debug("Could not find tunnel process: {}", e)

        return None

    def _terminate_process(self) -> None:
        """Terminate the SSH tunnel process gracefully."""
        if self._process is None:
            return

        pid = self._process.pid

        try:
            # Try graceful termination first
            os.kill(pid, signal.SIGTERM)
            logger.debug("Sent SIGTERM to SSH tunnel (PID: {})", pid)

            # Wait briefly for process to exit
            for _ in range(10):
                time.sleep(0.1)
                try:
                    os.kill(pid, 0)  # Check if process exists
                except OSError:
                    logger.debug("SSH tunnel process terminated (PID: {})", pid)
                    return

            # Force kill if still running
            logger.warning("SSH tunnel didn't respond to SIGTERM, sending SIGKILL")
            os.kill(pid, signal.SIGKILL)

        except OSError as e:
            # Process already gone
            logger.debug("Process {} already terminated: {}", pid, e)

    def _is_port_listening(self, host: str, port: int) -> bool:
        """Check if a port is listening and accepting connections."""
        try:
            with socket.create_connection((host, port), timeout=2.0):
                return True
        except (socket.error, socket.timeout, OSError):
            return False

    def _wait_for_port(self, host: str, port: int, timeout: float = 10.0) -> bool:
        """Wait for a port to become available."""
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self._is_port_listening(host, port):
                return True
            time.sleep(0.2)
        return False

    def _cleanup_on_exit(self) -> None:
        """Cleanup handler called on program exit."""
        if self._owned_tunnel and self._process is not None:
            logger.debug("Cleaning up SSH tunnel on exit")
            self._terminate_process()


class _PIDWrapper:
    """Minimal wrapper around a PID to provide a consistent interface."""

    def __init__(self, pid: int) -> None:
        self.pid = pid

    def poll(self) -> Optional[int]:
        """Check if process is still running. Returns None if running."""
        try:
            os.kill(self.pid, 0)
            return None
        except OSError:
            return -1


def ensure_tunnel(func: Callable[P, R]) -> Callable[P, R]:
    """
    Decorator that ensures SSH tunnel is active before calling the function.

    The tunnel will be started if not already active, and stopped after the
    function completes (only if this decorator started it).

    Usage
    -----
    ```python
    @ensure_tunnel
    def my_database_function():
        # Tunnel is guaranteed active here
        conn = psycopg2.connect(host="127.0.0.1", port=5433, ...)
    ```
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        with SSHTunnelManager():
            return func(*args, **kwargs)

    return wrapper


def get_tunnel_status() -> dict:
    """
    Get the current status of the SSH tunnel.

    Returns a dictionary with tunnel status information.
    """
    manager = SSHTunnelManager()
    cfg = manager.config

    status = {
        "configured": True,
        "ssh_host_alias": cfg.ssh_host_alias,
        "local_endpoint": f"{cfg.local_host}:{cfg.local_port}",
        "remote_endpoint": f"{cfg.remote_host}:{cfg.remote_port}",
        "is_active": manager.is_tunnel_active,
        "owned_by_manager": manager.owns_tunnel,
    }

    return status


def require_tunnel_or_fail() -> None:
    """
    Assert that the tunnel is active or raise an error.

    Use this at the start of functions that require database access to provide
    a clear error message if the tunnel isn't running.

    Raises
    ------
    TunnelError
        If the tunnel is not active.
    """
    manager = SSHTunnelManager()
    if not manager.is_tunnel_active:
        raise TunnelError(
            f"SSH tunnel is not active. Expected tunnel on "
            f"{manager.config.local_host}:{manager.config.local_port}. "
            f"Start the tunnel with: ssh -f -N {manager.config.ssh_host_alias}"
        )
