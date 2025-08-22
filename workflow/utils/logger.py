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


## Configure the logger ## -----------------------------------------------------
logging.basicConfig(
    format=LOGGER_FORMAT,
    datefmt=DATE_FORMAT,
    level=logging.INFO,
    stream=sys.stdout, # send to output (no red background)
)
for level, format in LEVEL_COLORS.items():
    logging.addLevelName(level, format + logging.getLevelName(level) + RESET)

# Log the current working directory, removing user-specific and irrelevant paths
logging.info("Working at %s", _hide_base_path(_repo_root()))

# Add project name to logger
logger_name = _project_name(_repo_root())
logger = logging.getLogger(logger_name)

## logger and shorthands ## ----------------------------------------------------
logger.setLevel(logging.INFO)
# Define functions for different logging levels
info = logger.info
warning = logger.warning
debug = logger.debug
error = logger.error
critical = logger.critical

## Suppress specific warnings ## -----------------------------------------------
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)
