import enum
import os
import pathlib
from typing import Union

import rtoml
from loguru import logger
from pydantic import BaseModel, HttpUrl

from activity_db_analysis.db.connector import (
    UnifiedPostgresConnector,
    build_unified_connector,
)


# ---------------------------------------------------------------------------
# Backwards Compatibility Wrappers
# ---------------------------------------------------------------------------
# Several modules still import `build_chembl_sql_handler_from_config` directly
# from this module (legacy pattern). We now provide a thin wrapper that delegates
# to the new AppConfig staticmethod returning a UnifiedPostgresConnector.
# This preserves import compatibility while the codebase is migrated.
def build_chembl_sql_handler_from_config(cli_overrides: dict | None = None):
    return AppConfig.build_chembl_sql_handler_from_config(cli_overrides=cli_overrides)


#
# Additional deprecated alias for legacy Activity DB builder name.
# Some modules may still attempt to import build_activity_db_sql_handler_from_config.
# Emit a warning and delegate to the new unified builder.
def build_activity_db_sql_handler_from_config(cli_overrides: dict | None = None):
    logger.warning(
        "Deprecated alias used: build_activity_db_sql_handler_from_config -> build_activity_sql_handler_from_config"
    )
    return AppConfig.build_activity_sql_handler_from_config(cli_overrides=cli_overrides)


class ConfigFileDictIsNone(Exception):
    def __str__(self):
        return (
            "Dictionary of config file parameters cannot be None. Cannot run software without configuration "
            "parameters."
        )


class FileFormatEnum(str, enum.Enum):
    pdf = "PDF"
    png = "PNG"


class AppConfig:
    config_file_parameters: dict = {}

    @staticmethod
    def set_config_file_parameters(parameters_from_config_file: dict) -> None:
        if parameters_from_config_file is None:
            raise ConfigFileDictIsNone
        AppConfig.config_file_parameters = parameters_from_config_file

    @staticmethod
    def get_config_dict() -> dict:
        return AppConfig.config_file_parameters

    class ChemblReleasesOptions(BaseModel):
        chembl_ftp: HttpUrl
        cache_file: str
        use_cache: bool

    chembl_releases_options: ChemblReleasesOptions = {}

    @staticmethod
    def set_chembl_releases_options(cli_options: dict) -> ChemblReleasesOptions:
        AppConfig.chembl_releases_options: AppConfig.ChemblReleasesOptions = (
            AppConfig.ChemblReleasesOptions(
                chembl_ftp=cli_options.get(
                    "ftp_url",
                    AppConfig.config_file_parameters.get("chembl", {}).get("ftp_url"),
                ),
                cache_file=cli_options.get(
                    "cache_file_name",
                    AppConfig.config_file_parameters.get("cache", {}).get(
                        "chembl_release_dates"
                    ),
                ),
                use_cache=cli_options.get(
                    "use_cache",
                    AppConfig.config_file_parameters.get("cache", {}).get(
                        "use_cache", False
                    ),
                ),
            )
        )

        return AppConfig.chembl_releases_options

    class ChemblDatabaseConnectionParameters(BaseModel):
        name: str
        user: str
        password: str
        host: str
        port: int
        db_schema: str

    chembl_db_connection_parameters: ChemblDatabaseConnectionParameters = {}

    @staticmethod
    def set_chembl_db_database_connection_parameters(
        cli_options: dict,
    ) -> ChemblDatabaseConnectionParameters:
        if not AppConfig.config_file_parameters:
            raise RuntimeError("No configuration loaded. Provide --app-config.")
        cfg = AppConfig.config_file_parameters
        chembl_section = cfg.get("chembl_db", {})
        # Enforce env-only password policy (ignore any provided CLI password override)
        if "password" in cli_options:
            logger.warning(
                "[AppConfig] Ignoring ChEMBL password provided via CLI override; environment variable CHEMBL_PWD is required."
            )
        resolved_chembl_pwd = os.environ.get("CHEMBL_PWD")
        if resolved_chembl_pwd is None:
            raise RuntimeError(
                "Missing required environment variable CHEMBL_PWD for ChEMBL database password."
            )
        logger.info(
            "[AppConfig] Using CHEMBL_PWD for ChEMBL user '{}'",
            cli_options.get("user", chembl_section.get("user")),
        )
        AppConfig.chembl_db_connection_parameters = (
            AppConfig.ChemblDatabaseConnectionParameters(
                name=cli_options.get("name", chembl_section.get("name")),
                user=cli_options.get("user", chembl_section.get("user")),
                password=resolved_chembl_pwd,
                host=cli_options.get("host", chembl_section.get("host")),
                port=cli_options.get("port", chembl_section.get("port")),
                db_schema=cli_options.get(
                    "schema", chembl_section.get("schema", "public")
                ),
            )
        )

        return AppConfig.chembl_db_connection_parameters

    class ActivityDatabaseConnectionParameters(BaseModel):
        name: str
        user: str
        password: str
        host: str
        port: int

    activity_db_connection_parameters: ActivityDatabaseConnectionParameters = {}

    @staticmethod
    def set_activity_db_database_connection_parameters(
        cli_options: dict,
    ) -> ActivityDatabaseConnectionParameters:
        if not AppConfig.config_file_parameters:
            raise RuntimeError("No configuration loaded. Provide --app-config.")
        cfg = AppConfig.config_file_parameters
        activity_section = cfg.get("activity_db", {})
        user_value = cli_options.get("user", activity_section.get("user"))
        # Derive expected password environment variable based on user suffix to ensure correct pairing.
        expected_pwd_env = None
        if isinstance(user_value, str):
            lv = user_value.lower()
            if lv.endswith("_dev"):
                expected_pwd_env = "ACT_DB_DEV_PWD"
            elif lv.endswith("_test"):
                expected_pwd_env = "ACT_DB_TEST_PWD"
            elif lv.endswith("_production") or lv.endswith("_prod"):
                expected_pwd_env = "ACT_DB_PROD_PWD"
        password_from_env = (
            os.environ.get(expected_pwd_env) if expected_pwd_env else None
        )
        # Enforce recognized user suffix and corresponding password env
        if expected_pwd_env is None:
            raise RuntimeError(
                f"Cannot resolve Activity DB password: user '{user_value}' must end with one of "
                "'_dev', '_test', '_production' (or '_prod')."
            )
        if password_from_env is None:
            raise RuntimeError(
                f"Expected environment variable {expected_pwd_env} for Activity DB user '{user_value}' is not set."
            )
        logger.info(
            f"[AppConfig] Using {expected_pwd_env} for Activity DB user '{user_value}'"
        )
        if "password" in cli_options:
            logger.warning(
                "[AppConfig] Ignoring Activity DB password provided via CLI override; environment variable {} is enforced.".format(
                    expected_pwd_env
                )
            )
        AppConfig.activity_db_connection_parameters = (
            AppConfig.ActivityDatabaseConnectionParameters(
                name=cli_options.get("name", activity_section.get("name")),
                user=user_value,
                password=password_from_env,
                host=cli_options.get("host", activity_section.get("host")),
                port=cli_options.get("port", activity_section.get("port")),
            )
        )

        return AppConfig.activity_db_connection_parameters

    @staticmethod
    def build_chembl_sql_handler_from_config(
        cli_overrides: dict | None = None,
    ) -> UnifiedPostgresConnector:
        """
        Build a unified PostgreSQL connector for the ChEMBL database using the resolved configuration.

        NOTE: SSH tunneling is managed at the CLI layer. Configure host=127.0.0.1 and the tunnel port
              in appconfig.toml when using the CLI-managed tunnel.
        """
        params = AppConfig.set_chembl_db_database_connection_parameters(
            cli_options=cli_overrides or {}
        )
        return build_unified_connector(
            db_name=params.name,
            user=params.user,
            password=params.password,
            host=params.host,
            port=params.port,
            schema=params.db_schema,
            auto_connect=True,
        )

    @staticmethod
    def build_activity_sql_handler_from_config(
        cli_overrides: dict | None = None,
    ) -> UnifiedPostgresConnector:
        """
        Build a unified PostgreSQL connector for the Activity database using the resolved configuration.

        NOTE: SSH tunneling is managed at the CLI layer. Configure host=127.0.0.1 and the tunnel port
              in appconfig.toml when using the CLI-managed tunnel.
        """
        params = AppConfig.set_activity_db_database_connection_parameters(
            cli_options=cli_overrides or {}
        )
        return build_unified_connector(
            db_name=params.name,
            user=params.user,
            password=params.password,
            host=params.host,
            port=params.port,
            schema=None,  # Activity DB does not currently expose schema differentiation
            auto_connect=True,
        )


def parse_app_config_toml(config_path: Union[str, pathlib.Path]) -> Union[dict, None]:
    """
    Safely parse a TOML configuration file and return its dictionary.

    Convenience:
        Accepts either a pathlib.Path or a str. A str will be coerced to Path early.

    Implementation notes:
    - Use text mode + explicit UTF-8 decoding.
    - Use rtoml.loads on the full string (rtoml.load expects already-parsed content or a file-like
      supporting .read returning str; passing a raw binary handle caused a failure).

    Side-effect:
        On successful parse, populates AppConfig.config_file_parameters so downstream
        builder methods (e.g. set_activity_db_database_connection_parameters) can proceed.
    """
    # Coerce string inputs
    if isinstance(config_path, str):
        config_path = pathlib.Path(config_path)

    if not config_path.exists():
        logger.error(f"Config file does not exist at {config_path}")
        return None
    if not config_path.is_file():
        logger.error(f"Provided path does not point to a file {config_path}")
        return None

    try:
        logger.trace(f"Parsing config file from {config_path}")
        with config_path.open("r", encoding="utf-8") as fh:
            toml_text = fh.read()
        data = rtoml.loads(toml_text)
        if not isinstance(data, dict):
            logger.error(f"Parsed config file did not yield a dictionary: {type(data)}")
            return None
        logger.debug(f"Loaded configuration keys: {list(data.keys())}")
        # Persist into AppConfig singleton state
        AppConfig.set_config_file_parameters(data)
        return data
    except Exception as e:
        logger.error(f"Failed to parse TOML config {config_path}: {e}")
        return None
