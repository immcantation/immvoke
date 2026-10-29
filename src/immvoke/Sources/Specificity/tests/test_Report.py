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
        self.assertEqual(Report.humanDuration(0.42), '0.4 s')
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

    def test_an_annotation_says_what_it_compares_and_what_it_found(self):
        """The run's own options, in words, rather than the raw values."""
        query = {'path': 'mine.tsv', 'rows': 43225}
        reference = {'path': 'databases/iedb/iedb_bcr.tsv', 'rows': 18650}
        options = {'sequence_region': 'CDR3H', 'partition_by': 'VJ',
                   'paired': False, 'identity_threshold': 0.7,
                   'distance_metric': 'hamming'}

        with self.assertLogs(Report.log, 'INFO') as logs:
            Report.annotationStarted(query, reference, options)
            Report.annotationCompleted(
                56.2, {'query_rows_compared': 21412, 'query_rows_with_hits': 17,
                       'hits': 22}, paired=False)
        text = '\n'.join(record.getMessage() for record in logs.records)

        self.assertIn('Annotation started: 43,225 rows of mine.tsv', text)
        self.assertIn('against 18,650 rows of databases/iedb/iedb_bcr.tsv', text)
        self.assertIn('comparing: the heavy chain CDR3 (cdr3_aa)', text)
        self.assertIn('a hit is: 70% identity or better, hamming', text)
        self.assertIn('only rows with the same V and J gene are compared', text)
        self.assertIn('Annotation completed in 56 s', text)
        self.assertIn('query rows matched         17 of 21,412', text)

    def test_an_exact_paired_annotation_says_so(self):
        options = {'sequence_region': 'CDR3H', 'partition_by': None,
                   'paired': False, 'identity_threshold': None,
                   'distance_metric': 'hamming'}
        with self.assertLogs(Report.log, 'INFO') as logs:
            Report.annotationStarted({'path': None, 'rows': 2},
                                     {'path': None, 'rows': 3}, options)
        text = '\n'.join(x.getMessage() for x in logs.records)

        self.assertIn('a table passed from Python', text)
        self.assertIn('a hit is: an exact match', text)
        self.assertIn('any pair of rows is compared', text)

        with self.assertLogs(Report.log, 'INFO') as logs:
            Report.annotationStarted({'path': None, 'rows': 2},
                                     {'path': None, 'rows': 3},
                                     dict(options, paired=True))
            Report.annotationCompleted(
                1.5, {'query_antibodies': 2, 'classification': {'both': 1, 'neither': 1}},
                paired=True)
        text = '\n'.join(x.getMessage() for x in logs.records)

        self.assertIn('heavy and light chain of the same antibody', text)
        self.assertIn('query antibodies', text)
        self.assertIn('both', text)

    def test_status_lines_carry_no_level_but_warnings_do(self):
        formatter = StatusFormatter('%(message)s')

        def line(level):
            return formatter.format(logging.LogRecord('immvoke', level, __file__, 1,
                                                      'hello', None, None))

        self.assertEqual(line(logging.INFO), 'hello')
        self.assertEqual(line(logging.WARNING), 'WARNING hello')


if __name__ == '__main__':
    unittest.main()
