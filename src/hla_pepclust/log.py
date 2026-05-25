"""Rich logging for the pipeline (modernized from the original cli/logger.py).

A single configured logger drives every stage; output goes through a
``RichHandler`` sharing the recording ``CONSOLE`` so it can be exported to a
log file. Use ``get_logger()`` to obtain it and ``save_console_log()`` to dump
the full coloured session to disk.
"""

from __future__ import annotations

import logging

from rich.logging import RichHandler

from hla_pepclust.console import CONSOLE

LOG_LEVELS = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}

_LOGGER_NAME = "hla_pepclust"


def configure_logging(
    level: str = "info",
    log_to_file: bool = False,
    file_name: str = "cluster_search_pipeline.log",
) -> logging.Logger:
    """Configure the package logger with a RichHandler (+ optional file handler)."""
    log_level = LOG_LEVELS.get(level.lower(), logging.INFO)
    handlers: list[logging.Handler] = [
        RichHandler(
            console=CONSOLE,
            rich_tracebacks=True,
            tracebacks_show_locals=True,
            markup=True,
            show_path=False,
            log_time_format="[%X]",
        )
    ]
    if log_to_file:
        fh = logging.FileHandler(file_name)
        fh.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "[%X]")
        )
        handlers.append(fh)

    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=handlers,
        force=True,
    )
    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(log_level)
    return logger


def get_logger() -> logging.Logger:
    """Return the package logger (configure_logging() should have run once)."""
    return logging.getLogger(_LOGGER_NAME)


def save_console_log(file_name: str = "cluster_search.log") -> None:
    """Write the full recorded console session (with styling stripped) to a file."""
    with open(file_name, "w") as fh:
        fh.write(CONSOLE.export_text())
