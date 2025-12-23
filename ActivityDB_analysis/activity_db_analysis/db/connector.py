"""
Unified PostgreSQL connector with context manager support.

Key Features
------------
- Single implementation for both ActivityDB and ChEMBL connections.
- Deterministic lifecycle: context manager (__enter__/__exit__) + explicit close().
- Consistent INFO logging (db/user, host:port).
- Convenience execution helpers: execute(...), run_scalar(...).
- Optional schema attribute (search_path set if provided).

This module is the authoritative database connectivity layer for the project.
It intentionally remains lightweight (no ORM, no pooling yet).

NOTE: SSH tunnel management is handled at the CLI layer (see db/ssh_tunnel.py).
      This connector expects to connect directly to the provided host:port,
      which should be localhost:<tunnel_port> when using the CLI-managed tunnel.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Sequence

import psycopg2
from loguru import logger


@dataclass(frozen=True)
class PostgresConnectionParams:
    """
    Immutable container for connection parameter values.
    """

    db_name: str
    user: str
    password: str
    host: str
    port: int = 5432
    schema: Optional[str] = None
    connect_timeout_sec: int = 10


class UnifiedPostgresConnector:
    """
    Lightweight PostgreSQL access layer with:
      - Deterministic context-managed lifecycle
      - Simple execution helpers (execute / run_scalar)
      - Optional schema search_path initialization

    Design goals: explicit, minimal, testable. No ORM, pooling, or retry logic
    (those can be layered on later if required).

    NOTE: SSH tunneling is managed at the CLI layer. This connector connects
          directly to the configured host:port (typically localhost via the
          configured tunnel port when using the CLI-managed tunnel).
    """

    def __init__(self, params: PostgresConnectionParams) -> None:
        self._params = params
        self._connection: Optional["psycopg2.extensions.connection"] = None
        self._closed: bool = False

    # --------------------------------------------------------------------- #
    # Properties
    # --------------------------------------------------------------------- #
    @property
    def params(self) -> PostgresConnectionParams:
        return self._params

    @property
    def host(self) -> str:
        """Host used for the psycopg2 connection."""
        return self._params.host

    @property
    def port(self) -> int:
        """Port used for the psycopg2 connection."""
        return self._params.port

    @property
    def raw_connection(self):
        """
        Direct psycopg2 connection object (may be None if not connected yet).
        """
        return self._connection

    @property
    def schema(self) -> Optional[str]:
        """
        Convenience accessor for the configured database schema (may be None).
        Mirrors params.schema to avoid reaching into ._params at call sites.
        """
        return self._params.schema

    # --------------------------------------------------------------------- #
    # Public API
    # --------------------------------------------------------------------- #
    def connect(self) -> "UnifiedPostgresConnector":
        """
        Establish the database connection.

        Safe to call multiple times; subsequent calls are no-ops once connected.
        """
        if self._connection is not None:
            return self  # Already connected

        p = self._params

        # Log summary before establishing connection
        logger.info(
            "Establishing PostgreSQL connection db='{}' user='{}' host='{}:{}' schema='{}'",
            p.db_name,
            p.user,
            p.host,
            p.port,
            p.schema or "(default)",
        )

        # Build DSN
        dsn = (
            f"dbname={p.db_name} user={p.user} password={p.password} "
            f"host={p.host} port={p.port} "
            f"connect_timeout={p.connect_timeout_sec}"
        )

        try:
            self._connection = psycopg2.connect(dsn)
        except Exception as e:
            self._abort_with_error(
                f"Failed to establish psycopg2 connection (host={p.host} port={p.port}): {e}"
            )

        logger.info(
            "PostgreSQL connection established (host='{}' port='{}')",
            p.host,
            p.port,
        )

        # Optional: set schema if requested
        if p.schema:
            try:
                with self._connection.cursor() as cur:
                    cur.execute(f"SET search_path TO {p.schema};")
                self._connection.commit()
                logger.debug("Set search_path to schema '{}'", p.schema)
            except Exception as e:
                logger.warning(
                    "Could not set search_path to schema '{}': {}", p.schema, e
                )

        return self

    def close(self) -> None:
        """
        Idempotent shutdown of DB connection.
        """
        if self._closed:
            return
        self._closed = True

        # Close DB connection
        if self._connection is not None:
            try:
                self._connection.close()
                logger.debug("Closed psycopg2 database connection")
            except Exception as e:
                logger.debug("Error closing database connection: {}", e)
            finally:
                self._connection = None

    # Context manager protocol
    def __enter__(self) -> "UnifiedPostgresConnector":
        return self.connect()

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    # Provide nicer debug representation
    def __repr__(self) -> str:  # pragma: no cover - trivial
        status = "connected" if self._connection else "disconnected"
        return (
            f"<UnifiedPostgresConnector db={self._params.db_name} user={self._params.user} "
            f"host={self._params.host}:{self._params.port} status={status}>"
        )

    # --------------------------------------------------------------------- #
    # Query helpers
    # --------------------------------------------------------------------- #
    def execute(
        self,
        query: str,
        params: Optional[Sequence[Any]] = None,
        fetch: bool = True,
    ) -> Optional[list[tuple[Any, ...]]]:
        """
        Execute a SQL statement.

        Parameters
        ----------
        query : str
            SQL text (param placeholders must align with psycopg2 style: %s).
        params : sequence | None
            Optional param values.
        fetch : bool
            If True, fetchall() is returned; otherwise None.

        Returns
        -------
        list[tuple] | None
        """
        self._assert_connected()
        assert self._connection is not None  # for type checkers

        logger.trace("Executing query (fetch={}):\n{}", fetch, query)
        try:
            with self._connection.cursor() as cur:
                cur.execute(query, params)
                if fetch:
                    rows = cur.fetchall()
                    logger.trace("Fetched {} row(s)", len(rows))
                    return rows
                else:
                    self._connection.commit()
                    return None
        except Exception as e:
            logger.error("Query execution failed: {}", e)
            raise

    def run_scalar(self, query: str, params: Optional[Sequence[Any]] = None) -> Any:
        """
        Execute a query expected to return a single row, single column.
        Raises ValueError if the shape is not as expected.
        """
        rows = self.execute(query, params=params, fetch=True)
        if not rows:
            raise ValueError("No rows returned for scalar query")
        if len(rows[0]) != 1:
            raise ValueError(
                f"Expected one column in scalar result, got {len(rows[0])}"
            )
        return rows[0][0]

    # --------------------------------------------------------------------- #
    # Internal helpers
    # --------------------------------------------------------------------- #
    def _abort_with_error(self, message: str) -> None:
        logger.error(message)
        # Attempt partial cleanup if something failed mid-setup
        try:
            self.close()
        except Exception:  # pragma: no cover - defensive
            pass
        raise RuntimeError(message)

    def _assert_connected(self) -> None:
        if self._connection is None:
            raise RuntimeError(
                "No active database connection. Call .connect() or use context manager."
            )


# --------------------------------------------------------------------- #
# Convenience factory
# --------------------------------------------------------------------- #
def build_unified_connector(
    *,
    db_name: str,
    user: str,
    password: str,
    host: str,
    port: int = 5432,
    schema: Optional[str] = None,
    connect_timeout_sec: int = 10,
    auto_connect: bool = True,
) -> UnifiedPostgresConnector:
    """
    Factory function to construct (and optionally connect) a UnifiedPostgresConnector.

    Parameters mirror PostgresConnectionParams plus auto_connect toggle.

    NOTE: SSH tunneling is managed at the CLI layer. When using the CLI,
          the tunnel is established before commands run, and you should
          configure host=127.0.0.1 and the tunnel port to connect through it.
    """
    params = PostgresConnectionParams(
        db_name=db_name,
        user=user,
        password=password,
        host=host,
        port=port,
        schema=schema,
        connect_timeout_sec=connect_timeout_sec,
    )
    connector = UnifiedPostgresConnector(params)
    if auto_connect:
        connector.connect()
    return connector
