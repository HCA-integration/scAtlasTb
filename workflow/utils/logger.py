#!/usr/bin/env python

"""
Logging configuration for the project.

This module sets up the logging format and levels for the project.
It uses the standard logging library to create a logger that outputs formatted
log messages to the console. It is recommended to use this module at the
beginning of your scripts to ensure consistent logging across the analysis.
"""

import os, sys, re, warnings
from pathlib import Path
import logging

## Config variables ## ---------------------------------------------------------

LOGGER_FORMAT = "[%(asctime)s] %(levelname)-8s %(name)s [%(filename)s:%(funcName)s:%(lineno)d] %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

LEVEL_COLORS = {
    logging.DEBUG:    "\033[1;34m",  # blue bold
    logging.INFO:     "\033[0;32m",  # green
    logging.WARNING:  "\033[1;33m",  # yellow bold
    logging.ERROR:    "\033[0;31m",  # red
    logging.CRITICAL: "\033[1;31m",  # red bold
}
RESET = "\033[0m"

## Functions ## ----------------------------------------------------------------

def _repo_root() -> Path | None:
    import subprocess
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ).stdout.strip()
        p = Path(out)
        return p if p.exists() else None
    except Exception:
        return None

def _project_name(repo_root: Path | None, exclude: str = "_|-") -> str:
    """Return styled project name from repo root or cwd."""
    name = (repo_root or Path.cwd()).name or "App"
    if re.findall(exclude, name): name = re.sub(exclude, " ", name).title()
    return name

def _hide_base_path(repo_root: Path | None, ignore: list[str] | None = None) -> str:
    """Return a cwd string with some superfluous base paths removed."""
    cwd = str(repo_root or Path.cwd().resolve())
    # Build default ignore list if none provided
    if ignore is None:
        user = os.environ.get('USER', os.environ.get('USERNAME', ''))
        ignore = [rf".*{re.escape(user)}"]
        ignore.extend([r"\.os\.py", r".*mamba", r".*conda", r".*projects"])
    # Apply filters sequentially
    for pat in ignore: cwd = re.sub(pat, "", cwd)
    # Fallback: collapse to last 2 parts if result is empty
    if not cwd.strip(): cwd = str(Path(*Path.cwd().parts[-2:]))
    return cwd

def _is_nested(logger: logging.Logger) -> bool:
    handlers = []
    for h in logger.handlers:
        temp = isinstance(h, logging.StreamHandler)
        handlers.append(temp and getattr(h, "_from_setup", False))
    return any(handlers)

class _ColorFormatter(logging.Formatter):
    def __init__(self, fmt: str, datefmt: str, use_color: bool) -> None:
        super().__init__(fmt=fmt, datefmt=datefmt)
        self.use_color = use_color

    def format(self, record: logging.LogRecord) -> str:
        if self.use_color:
            color = LEVEL_COLORS.get(record.levelno, "")
            if color:
                record.levelname = f"{color}{record.levelname}{RESET}"
        return super().format(record)

def _stream_supports_color(stream: object) -> bool:
    """Respect NO_COLOR; require TTY"""
    if os.environ.get("NO_COLOR"): return False
    return hasattr(stream, "isatty") and stream.isatty()

## Configure the logger ## -----------------------------------------------------

def setup_logger(level: int = logging.INFO, stream = sys.stdout) -> logging.Logger:
    logger_name = _project_name(_repo_root())
    # Root/logger setup without global side effects
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    # Avoid duplicate handlers if setup_logger() is called multiple times
    if not _is_nested(logger):
        handler = logging.StreamHandler(stream)
        handler._from_setup = True  # sentinel
        handler.setLevel(level)
        handler.setFormatter(
            _ColorFormatter(
                LOGGER_FORMAT, DATE_FORMAT, _stream_supports_color(stream)
            )
        )
        logger.addHandler(handler)
    # Show where the script is running from
    logger.info("Working at %s", _hide_base_path(_repo_root()))
    return logger

## logger and shorthands ## ----------------------------------------------------

logger = setup_logger(level=logging.INFO)
info = logger.info
warning = logger.warning
debug = logger.debug
error = logger.error
critical = logger.critical

## Suppress specific warnings ## -----------------------------------------------
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)
