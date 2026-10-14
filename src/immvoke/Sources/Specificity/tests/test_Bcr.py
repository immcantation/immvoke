"""
Unit tests for the `immvoke specificity bcr` shortcut (BcrCli.py)
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import io
import unittest
from pathlib import Path
from unittest import mock

# Immvoke imports
from immvoke.Cli import getArgParser
from immvoke.Sources.Specificity.BcrCli import DEFAULT_OUTDIR, handleBcr


class TestBcrParser(unittest.TestCase):
    """
    Tests for `bcr`'s parser wiring under `specificity`
    """

    def test_bcr_is_a_specificity_subcommand(self):
        parser = getArgParser()
        with mock.patch('sys.stdout', new_callable=io.StringIO) as out:
            with self.assertRaises(SystemExit):
                parser.parse_args(['specificity', 'bcr', '--help'])

        text = out.getvalue()
        self.assertIn('--outdir', text)
        # IEDB's own columns are always kept; there is nothing to strip.
        self.assertNotIn('--strict-airr', text)
        # No 'download' verb, no 'iedb' source name -- naming the table is enough.
        self.assertNotIn('download', text.split('\n')[0].lower())

    def test_outdir_is_optional(self):
        parser = getArgParser()
        args = parser.parse_args(['specificity', 'bcr'])

        self.assertEqual(args.outdir, DEFAULT_OUTDIR)


class TestHandleBcr(unittest.TestCase):
    """
    Tests for the handler, which delegates to handleSpecificityDownload
    """

    def test_delegates_to_handle_specificity_download_with_iedb_bcr(self):
        from argparse import Namespace

        seen = {}

        def _fake(args):
            seen.update(vars(args))
            return 0

        with mock.patch('immvoke.Cli.handleSpecificityDownload', side_effect=_fake):
            status = handleBcr(Namespace(outdir=Path('somewhere')))

        self.assertEqual(status, 0)
        self.assertEqual(seen['db'], 'iedb')
        self.assertEqual(seen['table'], 'bcr')
        self.assertEqual(seen['outdir'], Path('somewhere'))
        self.assertFalse(seen['strict_airr'])
        self.assertFalse(seen['dry_run'])
