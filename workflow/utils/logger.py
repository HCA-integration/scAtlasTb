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

## Config ## -------------------------------------------------------------------

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

# Configure the root logger
logging.basicConfig(
    format=LOGGER_FORMAT,
    datefmt=DATE_FORMAT,
    level=logging.INFO,
    stream=sys.stdout, # send to output (no red background)
)
for level, format in LEVEL_COLORS.items():
    logging.addLevelName(level, format + logging.getLevelName(level) + RESET)

# Log the current working directory, removing user-specific and irrelevant paths
base_paths = ".*" + os.environ.get('USER', os.environ.get('USERNAME'))
base_paths = base_paths + "|.os.py|.*mamba|.*conda|.*projects"
temp = re.sub(base_paths, "", Path().cwd().__str__())
logging.info("Working at %s", temp)

# Suppress specific warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)

# Add project name to logger
project_path = os.popen("git rev-parse --show-toplevel 2>&1").read().rstrip()
logger_name = Path(project_path).name if not 'fatal' in project_path else "App"
if re.findall("_|-", logger_name):
    logger_name = re.sub("_|-", " ", logger_name).title()
logger = logging.getLogger(logger_name)

## logger and shorthands ## ----------------------------------------------------
logger.setLevel(logging.INFO)
# Define functions for different logging levels
info = logger.info
warning = logger.warning
debug = logger.debug
error = logger.error
critical = logger.critical