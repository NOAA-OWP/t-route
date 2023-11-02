from pydantic import BaseModel

from typing_extensions import Literal

LogLevel = Literal["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]


class LoggingParameters(BaseModel, extra='forbid'):
    log_level: LogLevel = "DEBUG"
    """
    Python logging level. Can either be a string or an integer from the list below optional,
    defaults to DEBUG (10). All logging statements at or above the level specified will be
    displayed.

    To check the contents of the YAML configuration file, set the logging level at or lower than
    DEBUG For the most verbose output, set log_level to DEBUG
    """
    showtiming: bool = False
    """
    logical, it True a timing summary is provided that reports the total time required for each
    major process in the simulation sequence.  optional, defaults to None and no timing summary is
    reported
    """
