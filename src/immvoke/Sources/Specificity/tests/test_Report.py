"""
Unit tests for what `immvoke specificity <db> download` reports as it runs
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
import unittest
from pathlib import Path

import pandas

# Immvoke imports
from immvoke.Commandline import StatusFormatter
from immvoke.Sources.Specificity import Report
from immvoke.Sources.Specificity.Iedb import IedbSource


class TestReport(unittest.TestCase):
    """
    Tests for the download's status lines
    """

    def test_durations_read_well(self):
        self.assertEqual(Report.humanDuration(8.4), '8 s')
        self.assertEqual(Report.humanDuration(125), '2 min 5 s')
        self.assertEqual(Report.humanDuration(3 * 3600 + 610), '3 h 10 min')

    def test_tally_counts_distinct_values_as_chunks_pass(self):
        tally = Report.Tally(IedbSource, 'bcr')
        chunks = [pandas.DataFrame({'iedb_receptor_group': ['g1', 'g1'],
                                    'reactivity_ref': ['e1', '']}),
                  pandas.DataFrame({'iedb_receptor_group': ['g2'],
                                    'reactivity_ref': ['e1']})]

        passed = list(tally.wrap(chunks))

        self.assertEqual(len(passed), 2)
        self.assertEqual(tally.distinct(), [('receptor groups', 2), ('epitopes', 1)])

    def test_completed_says_what_was_written(self):
        tally = Report.Tally(IedbSource, 'bcr')
        list(tally.wrap([pandas.DataFrame({'iedb_receptor_group': ['g1'],
                                           'reactivity_ref': ['e1']})]))
        validation = {'header_valid': True, 'rows_checked': 1, 'rows_invalid': 0}

        with self.assertLogs(Report.log, 'INFO') as logs:
            Report.completed(Path('databases/iedb/iedb_bcr.tsv'), 15.2, 3, 1, tally,
                             validation)
        text = '\n'.join(record.getMessage() for record in logs.records)

        self.assertIn('Download completed in 15 s: databases/iedb/iedb_bcr.tsv', text)
        self.assertIn('rows             1 of 3 kept', text)
        self.assertIn('every row valid (iedb_bcr.airr_validation.log)', text)

    def test_status_lines_carry_no_level_but_warnings_do(self):
        formatter = StatusFormatter('%(message)s')

        def line(level):
            return formatter.format(logging.LogRecord('immvoke', level, __file__, 1,
                                                      'hello', None, None))

        self.assertEqual(line(logging.INFO), 'hello')
        self.assertEqual(line(logging.WARNING), 'WARNING hello')


if __name__ == '__main__':
    unittest.main()
