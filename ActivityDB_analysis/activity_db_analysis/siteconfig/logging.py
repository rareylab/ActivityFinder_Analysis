from enum import Enum


class LoggingLevel(str, Enum):
    """
    Define allowed logging levels for the application.

    Inherits from str for Typer CLI integration (enables autocompletion and validation).
    """

    TRACE = "TRACE"
    DEBUG = "DEBUG"
    INFO = "INFO"
    SUCCESS = "SUCCESS"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"
