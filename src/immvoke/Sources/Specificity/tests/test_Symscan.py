"""
Unit tests for symscan-accelerated annotation (Symscan.py)
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import random
import unittest
from unittest import mock

import pandas

# Immvoke imports
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity import Symscan
from immvoke.Sources.Specificity.Annotate import (
    runFuzzyAnnotation,
    runLevenshteinAnnotation,
)
from immvoke.Sources.Specificity.Symscan import (
    MAX_DISTANCE,
    runSymscanAnnotation,
)


def frame(rows):
    return pandas.DataFrame(rows)


class TestRunSymscanAnnotation(unittest.TestCase):
    """
    Tests for runSymscanAnnotation's criteria and result table
    """

    def test_hamming_within_budget_hits(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARQYW'}])  # 1 substitution

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])
        self.assertEqual(result['mismatch_positions'].tolist(), ['3'])

    def test_hamming_over_budget_does_not_hit(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CQRQYW'}])  # 2 substitutions

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(result['n_hits_total'].tolist(), [0])
        self.assertTrue(pandas.isna(result['min_mismatches'].iloc[0]))
        self.assertEqual(result['hit_ids'].tolist(), [''])

    def test_exact_cols_gate_the_pairs(self):
        queries = frame([{'cdr3': 'CARDYW', 'v_call': 'IGHV1-2'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW', 'v_call': 'IGHV3-9'}])

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id',
                                      exact_cols=['v_call'])

        self.assertEqual(result['n_hits_total'].tolist(), [0])

    def test_pooled_columns_must_line_up_column_by_column(self):
        """'AB'+'CDE' and 'ABC'+'DE' concatenate identically but are not a
        column-by-column match."""
        queries = frame([{'cdr1': 'AB', 'cdr3': 'CDE'}])
        reference = frame([{'id': 'r1', 'cdr1': 'ABC', 'cdr3': 'DE'},
                           {'id': 'r2', 'cdr1': 'AB', 'cdr3': 'CDF'}])

        result = runSymscanAnnotation(queries, reference, ['cdr1', 'cdr3'], 1,
                                      id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r2'])
        self.assertEqual(result['mismatch_positions'].tolist(), ['1:2'])

    def test_threshold_is_fraction_scales_with_query_length(self):
        queries = frame([{'cdr3': 'AAAAAAAAAA'}])  # length 10 -> 1 allowed at 0.1
        reference = frame([{'id': 'one', 'cdr3': 'AAAAAAAAAB'},
                           {'id': 'two', 'cdr3': 'AAAAAAAABB'}])

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 0.1, id_col='id',
                                      threshold_is_fraction=True)

        self.assertEqual(result['hit_ids'].tolist(), ['one'])

    def test_levenshtein_tolerates_an_indel(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDGYW'}])  # 1 insertion

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 1,
                                      distance_type='levenshtein', id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])
        self.assertNotIn('mismatch_positions', result.columns)

    def test_best_hit_per_id_is_kept(self):
        queries = frame([{'cdr3': 'AAAA'}])
        reference = frame([{'id': 'r1', 'cdr3': 'AABB'},   # distance 2
                           {'id': 'r1', 'cdr3': 'AAAB'}])  # distance 1, same id

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 2, id_col='id')

        self.assertEqual(result['n_hits_total'].tolist(), [1])
        self.assertEqual(result['min_mismatches'].tolist(), [1.0])
        self.assertEqual(result['mismatch_positions'].tolist(), ['3'])

    def test_wildcard_rows_are_compared_in_python_and_the_rest_by_symscan(self):
        """An 'X' matches anything, and only the rows holding one leave
        symscan: the clean query row still goes through it."""
        queries = frame([{'cdr3': 'CARXYW'}, {'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'},
                           {'id': 'r2', 'cdr3': 'CXRDYW'}])

        real = Symscan.symscan.get_neighbors_across
        with mock.patch.object(Symscan.symscan, 'get_neighbors_across',
                               side_effect=real) as spy:
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 0,
                                          id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1;r2', 'r1;r2'])
        self.assertEqual(spy.call_count, 1)
        self.assertEqual(spy.call_args.args, (['CARDYW'], ['CARDYW']))

    def test_levenshtein_compares_a_single_column(self):
        queries = frame([{'a': 'AA', 'b': 'BB'}])
        reference = frame([{'id': 'r1', 'a': 'AA', 'b': 'BB'}])

        with self.assertRaises(AnnotationError):
            runSymscanAnnotation(queries, reference, ['a', 'b'], 1,
                                 distance_type='levenshtein', id_col='id')

    def test_unknown_distance_type_raises(self):
        queries = frame([{'cdr3': 'AA'}])
        reference = frame([{'id': 'r1', 'cdr3': 'AA'}])

        with self.assertRaises(AnnotationError):
            runSymscanAnnotation(queries, reference, ['cdr3'], 1,
                                 distance_type='blosum', id_col='id')

    def _spy(self):
        return mock.patch.object(Symscan.symscan, 'get_neighbors_across',
                                 side_effect=Symscan.symscan.get_neighbors_across)

    def test_a_distance_symscan_would_crash_on_is_compared_in_python(self):
        """symscan panics (killing the process) above MAX_DISTANCE, so it
        must never be called with one."""
        queries = frame([{'cdr3': 'A' * 40}])
        reference = frame([{'id': 'r1', 'cdr3': 'B' * 15 + 'A' * 25}])

        with self._spy() as spy:
            result = runSymscanAnnotation(queries, reference, ['cdr3'],
                                          MAX_DISTANCE + 1, id_col='id')

        spy.assert_not_called()
        self.assertEqual(result['hit_ids'].tolist(), ['r1'])
        self.assertEqual(result['min_mismatches'].tolist(), [15.0])

    def test_a_reference_over_the_variant_budget_is_compared_in_python(self):
        for distance_type in ('hamming', 'levenshtein'):
            with self.subTest(distance_type), self._spy() as spy, \
                 mock.patch.object(Symscan, 'MAX_VARIANTS', 1):
                queries = frame([{'cdr3': 'CARDYW'}])
                reference = frame([{'id': 'r1', 'cdr3': 'CARQYW'}])

                result = runSymscanAnnotation(queries, reference, ['cdr3'], 1,
                                              distance_type=distance_type,
                                              id_col='id')

                spy.assert_not_called()
                self.assertEqual(result['hit_ids'].tolist(), ['r1'])

    def test_queries_are_chunked_to_fit_the_variant_budget(self):
        """At distance 1 a length-4 string has 5 variants: the reference
        (5) plus two queries (10) fit a budget of 15, so 4 queries take two
        calls."""
        queries = frame([{'cdr3': 'AAAA'}, {'cdr3': 'AAAC'}, {'cdr3': 'AACC'},
                         {'cdr3': 'CCCC'}])
        reference = frame([{'id': 'r1', 'cdr3': 'AAAA'}])

        with self._spy() as spy, mock.patch.object(Symscan, 'MAX_VARIANTS', 15):
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(spy.call_count, 2)
        self.assertEqual(result['hit_ids'].tolist(), ['r1', 'r1', '', ''])

    def test_hamming_only_indexes_same_length_references(self):
        queries = frame([{'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARQYW'},
                           {'id': 'r2', 'cdr3': 'CARDYWW'}])

        with self._spy() as spy:
            runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(spy.call_args.args, (['CARDYW'], ['CARQYW']))


class TestSymscanExactMatch(unittest.TestCase):
    """
    Tests for an exact match: runSymscanAnnotation at 0 substitutions
    """

    def test_exact_cols_gate_an_exact_match(self):
        queries = frame([{'cdr3': 'CARDYW', 'v_call': 'V1', 'j_call': 'J1'},
                         {'cdr3': 'CARDYW', 'v_call': 'V2', 'j_call': 'J1'},
                         {'cdr3': 'CARDY', 'v_call': 'V1', 'j_call': 'J1'}])
        reference = frame([{'id': 'r2', 'cdr3': 'CARDYW', 'v_call': 'V1', 'j_call': 'J1'},
                           {'id': 'r1', 'cdr3': 'CARDYW', 'v_call': 'V1', 'j_call': 'J1'}])

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 0, id_col='id',
                                      exact_cols=['v_call', 'j_call'])

        self.assertEqual(result['hit_ids'].tolist(), ['r1;r2', '', ''])
        self.assertEqual(result['is_unique_hit'].tolist(), [False, True, True])

    def test_blank_matches_blank(self):
        queries = frame([{'cdr3': ''}])
        reference = frame([{'id': 'r1', 'cdr3': ''}, {'id': 'r2', 'cdr3': 'A'}])

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 0, id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['r1'])


class TestSymscanMatchesPython(unittest.TestCase):
    """
    runSymscanAnnotation must return exactly the table its pure-Python
    counterpart returns -- checked on seeded random tables built to hit the
    edge cases (wildcards, blanks, duplicate ids, gene gates, pooled
    columns, fraction thresholds).
    """

    def _tables(self, seed):
        rng = random.Random(seed)

        def residues(n):
            return ''.join(rng.choice('ACDX' if rng.random() < 0.15 else 'ACD')
                           for _ in range(n))

        def row(i=None):
            values = {'cdr1': residues(rng.randint(0, 3)),
                      'cdr3': residues(rng.randint(3, 6)),
                      'v_call': rng.choice(['V1', 'V2']),
                      'j_call': rng.choice(['J1', 'J2'])}
            if i is not None:
                values['id'] = 'r%d' % (i % 50)  # some ids span several rows
            return values

        return (frame([row() for _ in range(40)]),
                frame([row(i) for i in range(70)]))

    def _budgets(self):
        """The real budget, then smaller ones that (on these tables) chunk
        the queries, send some groups to Python and keep others on symscan,
        and send most groups to Python."""
        return (Symscan.MAX_VARIANTS, 1500, 400, 150)

    def test_hamming(self):
        cases = [(['cdr3'], 1, None, False), (['cdr3'], 2, ['v_call', 'j_call'], False),
                 (['cdr1', 'cdr3'], 2, ['v_call'], False),
                 (['cdr3'], 0.34, ['v_call', 'j_call'], True),
                 (['cdr1', 'cdr3'], 0.25, None, True)]
        for seed in range(5):
            queries, reference = self._tables(seed)
            for cols, tolerance, exact_cols, fraction in cases:
                with self.subTest(seed=seed, cols=cols, tolerance=tolerance,
                                  exact_cols=exact_cols):
                    expected = runFuzzyAnnotation(queries, reference, cols, tolerance,
                                                  id_col='id', exact_cols=exact_cols,
                                                  threshold_is_fraction=fraction)
                    for budget in self._budgets():
                        with mock.patch.object(Symscan, 'MAX_VARIANTS', budget):
                            result = runSymscanAnnotation(
                                queries, reference, cols, tolerance, id_col='id',
                                exact_cols=exact_cols, threshold_is_fraction=fraction)
                        pandas.testing.assert_frame_equal(result, expected,
                                                          check_dtype=False)

    def test_levenshtein(self):
        for seed in range(5):
            queries, reference = self._tables(seed)
            for max_edits, exact_cols in ((1, None), (2, ['v_call', 'j_call'])):
                with self.subTest(seed=seed, max_edits=max_edits):
                    expected = runLevenshteinAnnotation(queries, reference, 'cdr3',
                                                        max_edits, id_col='id',
                                                        exact_cols=exact_cols)
                    for budget in self._budgets():
                        with mock.patch.object(Symscan, 'MAX_VARIANTS', budget):
                            result = runSymscanAnnotation(
                                queries, reference, ['cdr3'], max_edits,
                                'levenshtein', id_col='id', exact_cols=exact_cols)
                        pandas.testing.assert_frame_equal(result, expected,
                                                          check_dtype=False)

    def test_exact(self):
        """An exact match is zero wildcard-aware substitutions, gated on
        whichever genes are given."""
        for seed in range(5):
            queries, reference = self._tables(seed)
            for exact_cols in (None, ['v_call'], ['v_call', 'j_call']):
                with self.subTest(seed=seed, exact_cols=exact_cols):
                    expected = runFuzzyAnnotation(queries, reference, ['cdr3'], 0,
                                                  id_col='id', exact_cols=exact_cols)
                    result = runSymscanAnnotation(queries, reference, ['cdr3'], 0,
                                                  id_col='id', exact_cols=exact_cols)
                    pandas.testing.assert_frame_equal(result, expected,
                                                      check_dtype=False)


if __name__ == '__main__':
    unittest.main()
