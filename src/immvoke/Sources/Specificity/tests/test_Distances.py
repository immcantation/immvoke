"""
Unit tests for every gated pair's distance (Distances.py)
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import unittest
from unittest import mock

import pandas

# Immvoke imports
from immvoke.Sources.Specificity import Distances
from immvoke.Sources.Specificity.Distances import allPairDistances
from immvoke.Sources.Specificity.Symscan import runSymscanAnnotation
from immvoke.Sources.Specificity.tests.BruteForce import bruteForceAllPairs
from immvoke.Sources.Specificity.tests.test_Symscan import randomTables


def frame(rows):
    return pandas.DataFrame(rows)


def allPairs(*args, **kwargs):
    """allPairDistances' tables, as one."""
    return pandas.concat(list(allPairDistances(*args, **kwargs)), ignore_index=True)


class TestAllPairDistances(unittest.TestCase):
    """
    Tests for scoring every gated pair, whatever the threshold
    """

    def test_every_pair_is_scored_beyond_the_threshold(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'}, {'id': 'r2', 'cdr3': 'CQQQYW'},
                           {'id': 'r3', 'cdr3': 'CARD'}])

        pairs = allPairs(queries, reference, ['cdr3'], 1, distance_type='levenshtein',
                         id_col='id')

        self.assertEqual(pairs[['reference_id', 'mismatches', 'within_threshold']]
                         .values.tolist(),
                         [['r1', 0, True], ['r3', 2, False], ['r2', 3, False]])
        self.assertEqual(pairs['identity'].round(4).tolist(), [1.0, 0.6667, 0.5])

    def test_hamming_scores_same_length_pairs_only(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CQQQYW'}, {'id': 'r2', 'cdr3': 'CARD'},
                           {'id': 'r3', 'cdr3': 'XARDYF'}])

        pairs = allPairs(queries, reference, ['cdr3'], 0, id_col='id')

        # An X never counts under Hamming; CARD is another length.
        self.assertEqual(pairs[['reference_id', 'mismatches', 'mismatch_positions']]
                         .values.tolist(), [['r3', 1, '5'], ['r1', 3, '1,2,3']])

    def test_gates_still_apply(self):
        queries = frame([{'cdr3': 'CARDYW', 'v_call': 'IGHV1-2*02'},
                         {'cdr3': '', 'v_call': 'IGHV1-2*02'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CQQQQQ', 'v_call': 'IGHV1-2*04'},
                           {'id': 'r2', 'cdr3': 'CARDYW', 'v_call': 'IGHV3-23*01'},
                           {'id': 'r3', 'cdr3': 'CARDYW', 'v_call': ''}])

        pairs = allPairs(queries, reference, ['cdr3'], 0, id_col='id',
                         gene_cols=['v_call'])

        # Only the same gene, and nothing for the blank query.
        self.assertEqual(pairs[['query_row', 'reference_id', 'mismatches']]
                         .values.tolist(), [[0, 'r1', 5]])

    def test_no_pair_yields_one_empty_table(self):
        tables = list(allPairDistances(frame([{'cdr3': 'CARDYW'}]),
                                       frame([{'id': 'r1', 'cdr3': 'CA'}]), ['cdr3'], 0,
                                       id_col='id'))

        self.assertEqual(len(tables), 1)
        self.assertEqual(len(tables[0]), 0)
        self.assertEqual(tables[0].columns.tolist(),
                         ['query_row', 'reference_id', 'mismatches', 'identity',
                          'within_threshold', 'mismatch_positions'])


class TestAllPairDistancesMatchBruteForce(unittest.TestCase):
    """
    allPairDistances must return exactly the pairs a brute-force comparison
    of every pair returns, and its within_threshold pairs exactly symscan's
    -- on the seeded random tables test_Symscan checks symscan with.
    """

    CASES = [(['cdr3'], 1, {}),
             (['cdr3'], 2, {'exact_cols': ['j_call']}),
             (['cdr1', 'cdr3'], 2, {'gene_cols': ['v_call']}),
             (['cdr3'], 0.34, {'threshold_is_fraction': True,
                               'gene_cols': ['v_call', 'j_call']}),
             (['cdr3'], 1, {'distance_type': 'levenshtein'}),
             (['cdr3'], 2, {'distance_type': 'levenshtein',
                            'gene_cols': ['v_call', 'j_call']})]

    def _check(self, queries, reference, cols, tolerance, **kwargs):
        pairs = allPairs(queries, reference, cols, tolerance, id_col='id', **kwargs)
        expected = bruteForceAllPairs(queries, reference, cols, tolerance, id_col='id',
                                      **kwargs)
        pandas.testing.assert_frame_equal(pairs, expected, check_dtype=False)

        _, hits = runSymscanAnnotation(queries, reference, cols, tolerance, id_col='id',
                                       return_pairs=True, **kwargs)
        within = pairs[pairs['within_threshold']].drop(columns='within_threshold')
        pandas.testing.assert_frame_equal(within.reset_index(drop=True), hits,
                                          check_dtype=False)

    def test_random_tables(self):
        for seed in range(5):
            queries, reference = randomTables(seed)
            for cols, tolerance, kwargs in self.CASES:
                with self.subTest(seed=seed, cols=cols, tolerance=tolerance, **kwargs):
                    self._check(queries, reference, cols, tolerance, **kwargs)

    def test_small_blocks_give_the_same_pairs(self):
        """Query blocks and vectorized steps small enough to split every
        group many times over."""
        queries, reference = randomTables(0)
        with mock.patch.object(Distances, 'QUERY_BLOCK', 7), \
             mock.patch.object(Distances, 'RESIDUE_BUDGET', 50):
            for cols, tolerance, kwargs in self.CASES:
                with self.subTest(cols=cols, tolerance=tolerance, **kwargs):
                    self._check(queries, reference, cols, tolerance, **kwargs)


if __name__ == '__main__':
    unittest.main()
