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
from immvoke.Sources.Specificity.Symscan import (
    MAX_DISTANCE,
    runSymscanAnnotation,
)
from immvoke.Sources.Specificity.tests.BruteForce import bruteForceAnnotation


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

    def test_x_matches_any_residue_through_symscan_alone(self):
        """No pair is compared in Python: X rows go to symscan too."""
        queries = frame([{'cdr3': 'CARXYW'}, {'cdr3': 'CARDYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'},
                           {'id': 'r2', 'cdr3': 'CXRDYW'}])

        with self._spy() as spy:
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 0,
                                          id_col='id')

        self.assertTrue(spy.called)
        self.assertEqual(result['hit_ids'].tolist(), ['r1;r2', 'r1;r2'])
        self.assertEqual(result['min_mismatches'].tolist(), [0.0, 0.0])

    def test_x_slack_is_rescored_away(self):
        """CARXYW vs CAQDYW is 2 letters apart but 1 real difference: within
        the slack symscan is searched with, beyond an exact match."""
        queries = frame([{'cdr3': 'CARXYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CAQDYW'}])

        exact = runSymscanAnnotation(queries, reference, ['cdr3'], 0, id_col='id')
        one = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(exact['hit_ids'].tolist(), [''])
        self.assertEqual(one['hit_ids'].tolist(), ['r1'])
        self.assertEqual(one['min_mismatches'].tolist(), [1.0])
        self.assertEqual(one['mismatch_positions'].tolist(), ['2'])

    def test_x_is_an_ordinary_letter_under_levenshtein(self):
        queries = frame([{'cdr3': 'CARXYW'}])
        reference = frame([{'id': 'r1', 'cdr3': 'CARDYW'}])

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 0,
                                      distance_type='levenshtein', id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), [''])

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

    def _assertSkipped(self, queries, reference, *args, **kwargs):
        """symscan is never called; the query row gets no hit, with a
        warning naming it."""
        with self._spy() as spy, \
             self.assertLogs('immvoke.Sources.Specificity.Symscan', 'WARNING') as logs:
            result = runSymscanAnnotation(queries, reference, *args, id_col='id',
                                          **kwargs)

        spy.assert_not_called()
        self.assertEqual(result['hit_ids'].tolist(), [''])
        self.assertIn('query_row 0', logs.output[0])

    def test_a_distance_symscan_would_crash_on_is_skipped(self):
        """symscan panics (killing the process) above MAX_DISTANCE, so it
        must never be called with one."""
        self._assertSkipped(frame([{'cdr3': 'A' * 40}]),
                            frame([{'id': 'r1', 'cdr3': 'B' * 15 + 'A' * 25}]),
                            ['cdr3'], MAX_DISTANCE + 1)

    def test_reference_xs_are_cut_out_not_searched_around(self):
        """A reference's X positions are removed before symscan sees them:
        15 X's cost nothing, and the search stays exact."""
        queries = frame([{'cdr3': 'A' * 20}, {'cdr3': 'B' + 'A' * 19}])
        reference = frame([{'id': 'r1', 'cdr3': 'X' * 15 + 'A' * 5}])

        with self._spy() as spy:
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 0,
                                          id_col='id')

        self.assertEqual(spy.call_args.args, (['AAAAA'] * 2, ['AAAAA']))
        self.assertEqual(result['hit_ids'].tolist(), ['r1', 'r1'])
        self.assertEqual(result['mismatch_positions'].tolist(), ['', ''])

    def test_query_x_slack_beyond_symscan_is_skipped(self):
        """A query's own X's are covered by slack: 15 of them need 15."""
        self._assertSkipped(frame([{'cdr3': 'X' * 15 + 'A' * 5}]),
                            frame([{'id': 'r1', 'cdr3': 'A' * 20}]),
                            ['cdr3'], 0)

    def test_a_sequence_over_the_variant_budget_is_skipped(self):
        for distance_type in ('hamming', 'levenshtein'):
            with self.subTest(distance_type), \
                 mock.patch.object(Symscan, 'MAX_VARIANTS', 1):
                self._assertSkipped(frame([{'cdr3': 'CARDYW'}]),
                                    frame([{'id': 'r1', 'cdr3': 'CARQYW'}]),
                                    ['cdr3'], 1, distance_type=distance_type)

    def test_only_the_rows_that_do_not_fit_are_skipped(self):
        """At distance 1 a length-4 string has 5 variants and a length-12
        one 13: with half a budget of 20 the long reference cannot be
        indexed, but the short pair is still compared."""
        queries = frame([{'cdr3': 'AAAA'}])
        reference = frame([{'id': 'short', 'cdr3': 'AAAC'},
                           {'id': 'long', 'cdr3': 'AAAAAAAAAAAA'}])

        with mock.patch.object(Symscan, 'MAX_VARIANTS', 20), \
             self.assertLogs('immvoke.Sources.Specificity.Symscan', 'WARNING'):
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 1,
                                          distance_type='levenshtein', id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['short'])

    def test_queries_are_chunked_to_fit_the_variant_budget(self):
        """At distance 1 a length-4 string has 5 variants: each side gets
        half a budget of 20, so two queries per call -- 4 queries take two
        calls."""
        queries = frame([{'cdr3': 'AAAA'}, {'cdr3': 'AAAC'}, {'cdr3': 'AACC'},
                         {'cdr3': 'CCCC'}])
        reference = frame([{'id': 'r1', 'cdr3': 'AAAA'}])

        with self._spy() as spy, mock.patch.object(Symscan, 'MAX_VARIANTS', 20):
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(spy.call_count, 2)
        self.assertEqual(result['hit_ids'].tolist(), ['r1', 'r1', '', ''])

    def test_references_are_chunked_too(self):
        """A reference side larger than half the budget is split across
        calls rather than skipped."""
        queries = frame([{'cdr3': 'AAAA'}])
        reference = frame([{'id': 'r%d' % i, 'cdr3': cdr3}
                           for i, cdr3 in enumerate(['AAAA', 'AAAC', 'AACA', 'CCCC'])])

        with self._spy() as spy, mock.patch.object(Symscan, 'MAX_VARIANTS', 20):
            result = runSymscanAnnotation(queries, reference, ['cdr3'], 1, id_col='id')

        self.assertEqual(spy.call_count, 2)
        self.assertEqual(result['hit_ids'].tolist(), ['r0;r1;r2'])

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

    def test_a_blank_value_never_matches(self):
        """Not even another blank: there is nothing to compare."""
        queries = frame([{'cdr3': ''}, {'cdr3': 'A'}])
        reference = frame([{'id': 'r1', 'cdr3': ''}, {'id': 'r2', 'cdr3': 'A'}])

        result = runSymscanAnnotation(queries, reference, ['cdr3'], 0, id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), ['', 'r2'])

    def test_a_blank_in_any_pooled_column_never_matches(self):
        queries = frame([{'cdr1': '', 'cdr3': 'CAR'}])
        reference = frame([{'id': 'r1', 'cdr1': '', 'cdr3': 'CAR'}])

        result = runSymscanAnnotation(queries, reference, ['cdr1', 'cdr3'], 1,
                                      id_col='id')

        self.assertEqual(result['hit_ids'].tolist(), [''])


class TestSymscanMatchesBruteForce(unittest.TestCase):
    """
    runSymscanAnnotation must return exactly the table a brute-force
    comparison of every pair returns -- checked on seeded random tables
    built to hit the edge cases (X residues, blanks, duplicate ids, exact
    and gene gates, pooled columns, fraction thresholds).
    """

    def _tables(self, seed):
        rng = random.Random(seed)

        def residues(n):
            value = ''.join(rng.choice('ACD') for _ in range(n))
            if value and rng.random() < 0.2:  # at most one X per column
                pos = rng.randrange(len(value))
                value = value[:pos] + 'X' + value[pos + 1:]
            return value

        def row(i=None):
            values = {'cdr1': residues(rng.randint(0, 3)),
                      'cdr3': residues(rng.randint(3, 6)),
                      'v_call': rng.choice(['IGHV1-2*01', 'IGHV1-2*02', 'IGHV3-23*01',
                                            '']),
                      'j_call': rng.choice(['IGHJ4*02', 'IGHJ6*01'])}
            if i is not None:
                values['id'] = 'r%d' % (i % 50)  # some ids span several rows
            return values

        return (frame([row() for _ in range(40)]),
                frame([row(i) for i in range(70)]))

    def _check(self, queries, reference, *args, **kwargs):
        expected = bruteForceAnnotation(queries, reference, *args, id_col='id',
                                        **kwargs)
        result = runSymscanAnnotation(queries, reference, *args, id_col='id',
                                      **kwargs)
        pandas.testing.assert_frame_equal(result, expected, check_dtype=False)

    def test_hamming(self):
        cases = [(['cdr3'], 0, {}), (['cdr3'], 1, {}),
                 (['cdr3'], 2, {'exact_cols': ['j_call']}),
                 (['cdr1', 'cdr3'], 2, {'gene_cols': ['v_call']}),
                 (['cdr3'], 0.34, {'threshold_is_fraction': True,
                                   'gene_cols': ['v_call', 'j_call'],
                                   'gene_match': 'allele'}),
                 (['cdr1', 'cdr3'], 0.25, {'threshold_is_fraction': True})]
        for seed in range(5):
            queries, reference = self._tables(seed)
            for cols, tolerance, kwargs in cases:
                with self.subTest(seed=seed, cols=cols, tolerance=tolerance, **kwargs):
                    self._check(queries, reference, cols, tolerance, **kwargs)

    def test_levenshtein(self):
        for seed in range(5):
            queries, reference = self._tables(seed)
            for max_edits, kwargs in ((1, {}), (2, {'gene_cols': ['v_call', 'j_call']})):
                with self.subTest(seed=seed, max_edits=max_edits):
                    self._check(queries, reference, ['cdr3'], max_edits,
                                distance_type='levenshtein', **kwargs)

    def test_chunked_queries_give_the_same_result(self):
        """A budget small enough to split every group's queries over many
        symscan calls, yet still fitting each single query."""
        queries, reference = self._tables(0)
        with mock.patch.object(Symscan, 'MAX_VARIANTS', 2000):
            self._check(queries, reference, ['cdr3'], 1)


if __name__ == '__main__':
    unittest.main()
