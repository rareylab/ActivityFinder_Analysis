"""
Database connectivity module for ActivityDB analysis.

This module provides:
- UnifiedPostgresConnector: PostgreSQL connection layer
- SSHTunnelManager: Native SSH tunnel management using subprocess
- CLI integration helpers for tunnel management

Architecture
------------
SSH tunnel management is centralized at the CLI layer:

1. CLI callback (cli.py) starts the tunnel via SSHTunnelManager before commands run
2. Connectors connect to localhost via the configured tunnel port
3. CLI callback stops the tunnel after commands complete (if it started it)

Key Components
--------------
Connector (connector.py):
    - PostgresConnectionParams: Connection parameter container
    - UnifiedPostgresConnector: Main database connector class
    - build_unified_connector: Factory function

SSH Tunnel (ssh_tunnel.py):
    - SSHTunnelManager: Context manager for SSH tunnel lifecycle
    - TunnelConfig: Tunnel configuration dataclass
    - ensure_tunnel: Decorator for functions requiring tunnel
    - get_tunnel_status: Query current tunnel state

CLI Integration (cli_tunnel_wrapper.py):
    - with_tunnel: Decorator for CLI commands
    - TunnelMode: Enum for tunnel handling modes
    - tunnel_callback: Typer callback for tunnel management
"""

from activity_db_analysis.db.cli_tunnel_wrapper import (
    TunnelAwareCommand,
    TunnelMode,
    print_tunnel_status,
    tunnel_callback,
    with_tunnel,
)
from activity_db_analysis.db.connector import (
    PostgresConnectionParams,
    UnifiedPostgresConnector,
    build_unified_connector,
)
from activity_db_analysis.db.ssh_tunnel import (
    SSHTunnelManager,
    TunnelConfig,
    TunnelError,
    TunnelHealthCheckError,
    TunnelStartError,
    ensure_tunnel,
    get_tunnel_status,
    require_tunnel_or_fail,
)

__all__ = [
    # Connector
    "PostgresConnectionParams",
    "UnifiedPostgresConnector",
    "build_unified_connector",
    # SSH Tunnel
    "SSHTunnelManager",
    "TunnelConfig",
    "TunnelError",
    "TunnelStartError",
    "TunnelHealthCheckError",
    "ensure_tunnel",
    "get_tunnel_status",
    "require_tunnel_or_fail",
    # CLI Integration
    "TunnelMode",
    "with_tunnel",
    "tunnel_callback",
    "print_tunnel_status",
    "TunnelAwareCommand",
]
