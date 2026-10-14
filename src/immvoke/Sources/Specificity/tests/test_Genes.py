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
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity.Genes import (
    GeneGate,
    callsAgree,
    inferLocus,
    parseCall,
)


def agree(a, b, level):
    return callsAgree(parseCall(a, level), parseCall(b, level))


class TestParseCall(unittest.TestCase):
    """
    Tests for reading the names out of one gene call cell
    """

    def test_levels(self):
        call = 'IGHV1-2*02'
        self.assertEqual(parseCall(call, 'allele'), {'IGHV1-2': {'02'}})
        self.assertEqual(parseCall(call, 'gene'), {'IGHV1-2': {None}})
        self.assertEqual(parseCall(call, 'family'), {'IGHV1': {None}})

    def test_several_calls_in_one_cell(self):
        self.assertEqual(parseCall('IGHV5-9*01, IGHV5-9*03, IGHV5-9-5*01', 'allele'),
                         {'IGHV5-9': {'01', '03'}, 'IGHV5-9-5': {'01'}})

    def test_a_gene_level_name_has_no_allele(self):
        """Not gene 'IGHV1' allele '2', as Change-O's allele pattern reads it --
        and, at allele level, nothing to agree on."""
        self.assertEqual(parseCall('IGHV1-2', 'allele'), {'IGHV1-2': set()})

    def test_names_change_o_parses(self):
        for call, gene in [('IGHV3/OR16-10*02', 'IGHV3/OR16-10'),
                           ('IGHV1-69D*01', 'IGHV1-69D'),
                           ('IGKV1D-39*01', 'IGKV1D-39'),
                           ('IGHV3-30-3*01', 'IGHV3-30-3'),
                           ('TRBV20-1*01', 'TRBV20-1')]:
            with self.subTest(call):
                self.assertEqual(list(parseCall(call, 'gene')), [gene])

    def test_blank_or_unrecognized_is_unknown(self):
        for call in ('', None, 'IGVL3-10', 'unknown'):
            with self.subTest(call):
                self.assertEqual(parseCall(call, 'gene'), {})

    def test_unknown_level_raises(self):
        with self.assertRaises(AnnotationError):
            parseCall('IGHV1-2*02', 'segment')


class TestCallsAgree(unittest.TestCase):
    """
    Tests for whether two cells name the same gene at each level
    """

    def test_allele_level(self):
        self.assertTrue(agree('IGHV1-2*02', 'IGHV1-2*02', 'allele'))
        self.assertFalse(agree('IGHV1-2*02', 'IGHV1-2*04', 'allele'))
        # a gene level call is a gene level match, never an allele level one
        self.assertFalse(agree('IGHV1-2*02', 'IGHV1-2', 'allele'))
        self.assertFalse(agree('IGHV1-2', 'IGHV1-2', 'allele'))
        self.assertTrue(agree('IGHV1-2, IGHV1-2*02', 'IGHV1-2*02', 'allele'))

    def test_gene_level(self):
        self.assertTrue(agree('IGHV1-2*02', 'IGHV1-2*04', 'gene'))
        self.assertTrue(agree('IGHV1-2*02', 'IGHV1-2', 'gene'))
        self.assertFalse(agree('IGHV1-2*02', 'IGHV1-69*01', 'gene'))
        self.assertFalse(agree('IGHV1-69*01', 'IGHV1-69D*01', 'gene'))

    def test_family_level(self):
        self.assertTrue(agree('IGHV1-2*02', 'IGHV1-69*01', 'family'))
        self.assertFalse(agree('IGHV1-2*02', 'IGHV3-23*01', 'family'))

    def test_any_shared_call_agrees(self):
        self.assertTrue(agree('IGHV3-23*01, IGHV1-2*02', 'IGHV1-2*02', 'allele'))
        self.assertFalse(agree('IGHV3-23*01, IGHV1-2*02', 'IGHV1-69*01', 'gene'))

    def test_blank_or_unrecognized_agrees_with_nothing(self):
        for level in ('allele', 'gene', 'family'):
            for other in ('IGHV1-2*02', '', 'IGVL3-10'):
                with self.subTest(level=level, other=other):
                    self.assertFalse(agree('', other, level))
                    self.assertFalse(agree('IGVL3-10', other, level))


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

        gate = GeneGate(queries, reference, ['v_call', 'j_call'], 'gene')

        keep = gate.keep(numpy.array([0, 0, 0]), numpy.array([0, 1, 2]))
        self.assertEqual(keep.tolist(), [True, False, False])

    def test_unrecognized_values_are_logged(self):
        frame = pandas.DataFrame([{'v_call': 'IGVL3-10'}])
        with self.assertLogs('immvoke.Sources.Specificity.Genes', 'WARNING'):
            GeneGate(frame, frame, ['v_call'])


if __name__ == '__main__':
    unittest.main()
