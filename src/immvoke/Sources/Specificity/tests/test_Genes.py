"""
Unit tests for V/J-gene call parsing and comparison (Genes.py)
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import unittest

import numpy
import pandas

# Immvoke imports
from immvoke.Sources.Specificity.Genes import (
    GeneGate,
    callsAgree,
    inferLocus,
    parseCall,
)


def agree(a, b):
    return callsAgree(parseCall(a), parseCall(b))


class TestParseCall(unittest.TestCase):
    """
    Tests for reading the gene names out of one gene call cell
    """

    def test_the_allele_is_dropped(self):
        self.assertEqual(parseCall('IGHV1-2*02'), {'IGHV1-2'})
        self.assertEqual(parseCall('IGHV1-2'), {'IGHV1-2'})

    def test_several_calls_in_one_cell(self):
        self.assertEqual(parseCall('IGHV5-9*01, IGHV5-9*03, IGHV5-9-5*01'),
                         {'IGHV5-9', 'IGHV5-9-5'})

    def test_a_gene_level_name_is_not_read_as_an_allele(self):
        """Not gene 'IGHV1' allele '2', as Change-O's allele pattern reads it."""
        self.assertEqual(parseCall('IGHV1-2'), {'IGHV1-2'})

    def test_names_change_o_parses(self):
        for call, gene in [('IGHV3/OR16-10*02', 'IGHV3/OR16-10'),
                           ('IGHV1-69D*01', 'IGHV1-69D'),
                           ('IGKV1D-39*01', 'IGKV1D-39'),
                           ('IGHV3-30-3*01', 'IGHV3-30-3'),
                           ('TRBV20-1*01', 'TRBV20-1')]:
            with self.subTest(call):
                self.assertEqual(parseCall(call), {gene})

    def test_blank_or_unrecognized_names_nothing(self):
        for call in ('', None, 'IGVL3-10', 'unknown'):
            with self.subTest(call):
                self.assertEqual(parseCall(call), frozenset())


class TestCallsAgree(unittest.TestCase):
    """
    Tests for whether two cells name the same gene, any allele
    """

    def test_same_gene_any_allele(self):
        self.assertTrue(agree('IGHV1-2*02', 'IGHV1-2*04'))
        self.assertTrue(agree('IGHV1-2*02', 'IGHV1-2'))
        self.assertFalse(agree('IGHV1-2*02', 'IGHV1-69*01'))
        self.assertFalse(agree('IGHV1-69*01', 'IGHV1-69D*01'))

    def test_any_shared_call_agrees(self):
        self.assertTrue(agree('IGHV3-23*01, IGHV1-2*02', 'IGHV1-2*04'))
        self.assertFalse(agree('IGHV3-23*01, IGHV1-2*02', 'IGHV1-69*01'))

    def test_blank_or_unrecognized_agrees_with_nothing(self):
        for other in ('IGHV1-2*02', '', 'IGVL3-10'):
            with self.subTest(other=other):
                self.assertFalse(agree('', other))
                self.assertFalse(agree('IGVL3-10', other))


class TestInferLocus(unittest.TestCase):
    """
    Tests for filling a blank locus from the chain's own gene calls
    """

    def test_the_v_call_names_the_locus(self):
        self.assertEqual(inferLocus('IGKV1-39*01', 'IGKJ1*01'), 'IGK')
        self.assertEqual(inferLocus('IGLV2-14*01, IGLV2-23*02', ''), 'IGL')

    def test_the_j_call_is_used_without_a_v_call(self):
        self.assertEqual(inferLocus('', 'IGLJ2*01'), 'IGL')

    def test_nothing_is_guessed(self):
        cases = [('', ''),                          # no gene call
                 ('IGVL3-10', ''),                  # unrecognized
                 ('IGKV1-39*01', 'IGLJ2*01'),       # V and J disagree
                 ('IGKV1-39*01, IGLV2-14*01', '')]  # V names two loci
        for v_call, j_call in cases:
            with self.subTest(v_call=v_call, j_call=j_call):
                self.assertEqual(inferLocus(v_call, j_call), '')


class TestGeneGate(unittest.TestCase):
    """
    Tests for the pairwise gate over several gene columns
    """

    def test_every_column_must_agree(self):
        queries = pandas.DataFrame([{'v_call': 'IGHV1-2*01', 'j_call': 'IGHJ4*02'}])
        reference = pandas.DataFrame([
            {'v_call': 'IGHV1-2*02', 'j_call': 'IGHJ4*02'},   # both agree
            {'v_call': 'IGHV1-2*02', 'j_call': 'IGHJ6*02'},   # J differs
            {'v_call': 'IGHV1-2', 'j_call': ''},              # J blank
        ])

        gate = GeneGate(queries, reference, ['v_call', 'j_call'])

        keep = gate.keep(numpy.array([0, 0, 0]), numpy.array([0, 1, 2]))
        self.assertEqual(keep.tolist(), [True, False, False])

    def test_unrecognized_values_are_logged(self):
        frame = pandas.DataFrame([{'v_call': 'IGVL3-10'}])
        with self.assertLogs('immvoke.Sources.Specificity.Genes', 'WARNING'):
            GeneGate(frame, frame, ['v_call'])


if __name__ == '__main__':
    unittest.main()
