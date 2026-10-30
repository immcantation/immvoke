"""
Commandline interface helpers
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import logging
import sys
from argparse import ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter


class CommonHelpFormatter(RawDescriptionHelpFormatter, ArgumentDefaultsHelpFormatter):
    """
    Custom argparse.HelpFormatter preserving epilog layout and showing defaults.

    Matches the formatter used across Change-O and pRESTO so that ``immvoke``
    help output reads like the rest of Immcantation.
    """
    pass


def setupLogging(verbose=False, quiet=False):
    """
    Configure package logging.

    Progress bars go to stderr via tqdm; this configures everything else. Status
    output uses the logging module rather than print so that it can be captured,
    redirected and silenced.

    Arguments:
      verbose (bool): if True set the level to DEBUG.
      quiet (bool): if True set the level to ERROR. Takes precedence over verbose.

    Returns:
      logging.Logger: the configured package logger.
    """
    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(StatusFormatter('%(message)s'))
    logging.basicConfig(level=level, handlers=[handler])

    return logging.getLogger('immvoke')


class StatusFormatter(logging.Formatter):
    """
    Write status (INFO) lines as plain text, and name the level of every
    other line: 'WARNING ...', 'ERROR ...', 'DEBUG ...'.
    """

    def format(self, record):
        text = super().format(record)
        if record.levelno == logging.INFO:
            return text

        return '%s %s' % (record.levelname, text)
