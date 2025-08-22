#!/usr/bin/env python

"""
Logging configuration for the project.

This module sets up the logging format and level for the project.
It uses the standard logging library to create a logger that outputs formatted
log messages to the console. It is recommended to use this module at the
beginning of your scripts to ensure consistent logging across the analysis.
"""

import os, sys, re, warnings

import logging

logger_format = "[%(asctime)s] %(name)s %(levelname)-2s [%(filename)s:%(funcName)s:%(lineno)s] %(message)s"
FORMATS = {
    logging.DEBUG: "\033[1;34m", # blue bold
    logging.INFO: "\033[0;32m", # green
    logging.WARNING: "\033[1;33m", # yellow bold
    logging.ERROR: "\033[0;31m", # red
    logging.CRITICAL: "\033[1;31m", # red bold
}

# temp = ["__vsc_ipynb_file__"] # when working on vscode
# __file__ = [globals()[i] for i in temp if i in globals().keys()][0]
# os.chdir(os.path.dirname(os.path.dirname(__file__)))
logging.basicConfig(
    format=logger_format,
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    stream=sys.stdout, # send to output (no red background)
)
for level, format in FORMATS.items():
    logging.addLevelName(level, format + logging.getLevelName(level) + "\033[0m")

# Name logger with the project name
project_path = os.popen("git rev-parse --show-toplevel 2>&1").read().rstrip()
logger = logging.getLogger(os.path.basename(project_path))

# Set up logging levels
logger.setLevel(logging.INFO)

info = logger.info
warning = logger.warning
debug = logger.debug
error = logger.error
critical = logger.critical