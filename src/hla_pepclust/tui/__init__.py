"""Terminal UI: Rich console, banner, results table, logging, error messages."""

from hla_pepclust.tui.console import CONSOLE, banner, results_table
from hla_pepclust.tui.errors import (
    file_not_found,
    png_not_found,
    reference_not_found,
)
from hla_pepclust.tui.log import (
    LOG_LEVELS,
    configure_logging,
    get_logger,
    save_console_log,
)

__all__ = [
    "CONSOLE",
    "banner",
    "results_table",
    "configure_logging",
    "get_logger",
    "save_console_log",
    "LOG_LEVELS",
    "file_not_found",
    "png_not_found",
    "reference_not_found",
]
