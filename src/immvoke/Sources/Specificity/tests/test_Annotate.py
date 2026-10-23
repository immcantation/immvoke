"""
Unit tests for antigen specificity annotation (Annotate.py, AnnotateCli.py)
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import io
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy
import pandas
import yaml
from airr.schema import RearrangementSchema

# Immvoke imports
from immvoke.Cli import getArgParser
from immvoke.Convert import writeAirr
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity.Annotate import (
    ANNOTATION_FILE,
    ANNOTATION_METADATA,
    PAIRS_FILE,
    annotateSpecificity,
    runPairedExactAnnotation,
)
from immvoke.Sources.Specificity.AnnotateCli import handleAnnotate
from immvoke.Sources.Specificity.tests.BruteForce import bruteForceAnnotation

ANNOTATE = 'immvoke.Sources.Specificity.Annotate'


def frame(rows):
    return pandas.DataFrame(rows)


class TestRunPairedExactAnnotation(unittest.TestCase):
    """
    Tests for the heavy+light reconciliation (exact V-gene, J-gene and CDR3
    agreement, both chains against the same reference antibody)
    """

    def _antibodies(self):
        """
        5 query antibodies, one per classification, against 2 reference
        antibodies -- see the inline comments for which reference row each
        is built to agree or disagree with.
        """
        queries = frame([
            # q1: heavy and light both agree with refA -> both
            {'cell_id': 'q1', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'q1', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
             'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            # q2: heavy agrees with refA, light agrees with nothing -> heavy_only
            {'cell_id': 'q2', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'q2', 'locus': 'IGL', 'cdr3_aa': 'CDIFFER',
             'v_call': 'IGKV9', 'j_call': 'IGKJ9'},
            # q3: heavy agrees with nothing, light agrees with refA -> light_only
            {'cell_id': 'q3', 'locus': 'IGH', 'cdr3_aa': 'CDIFF2',
             'v_call': 'IGHV9', 'j_call': 'IGHJ9'},
            {'cell_id': 'q3', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
             'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            # q4: neither chain agrees with anything -> neither
            {'cell_id': 'q4', 'locus': 'IGH', 'cdr3_aa': 'CDIFF3',
             'v_call': 'IGHV8', 'j_call': 'IGHJ8'},
            {'cell_id': 'q4', 'locus': 'IGL', 'cdr3_aa': 'CDIFF4',
             'v_call': 'IGKV8', 'j_call': 'IGKJ8'},
            # q5: heavy agrees with refA, light agrees with refB -> conflicting
            {'cell_id': 'q5', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'q5', 'locus': 'IGL', 'cdr3_aa': 'CZZZZZZ',
             'v_call': 'IGKV7', 'j_call': 'IGKJ7'},
        ])
        reference = frame([
            {'cell_id': 'refA', 'locus': 'IGH', 'cdr3_aa': 'CARDYW',
             'v_call': 'IGHV1', 'j_call': 'IGHJ1'},
            {'cell_id': 'refA', 'locus': 'IGL', 'cdr3_aa': 'CQQYNS',
             'v_call': 'IGKV1', 'j_call': 'IGKJ1'},
            {'cell_id': 'refB', 'locus': 'IGH', 'cdr3_aa': 'CWWWWWW',
             'v_call': 'IGHV6', 'j_call': 'IGHJ6'},
            {'cell_id': 'refB', 'locus': 'IGL', 'cdr3_aa': 'CZZZZZZ',
             'v_call': 'IGKV7', 'j_call': 'IGKJ7'},
        ])
        return queries, reference

    def test_classifications(self):
        queries, reference = self._antibodies()

        result = runPairedExactAnnotation(queries, reference)

        by_cell_id = result.set_index('cell_id')['classification'].to_dict()
        self.assertEqual(by_cell_id, {
            'q1': 'both', 'q2': 'heavy_only', 'q3': 'light_only',
            'q4': 'neither', 'q5': 'conflicting',
        })

    def test_both_hit_ids_is_the_intersection_not_the_union(self):
        queries, reference = self._antibodies()

        result = runPairedExactAnnotation(queries, reference)

        row = result.set_index('cell_id').loc['q1']
        self.assertEqual(row['heavy_hit_ids'], 'refA')
        self.assertEqual(row['light_hit_ids'], 'refA')
        self.assertEqual(row['both_hit_ids'], 'refA')

        row5 = result.set_index('cell_id').loc['q5']
        self.assertEqual(row5['heavy_hit_ids'], 'refA')
        self.assertEqual(row5['light_hit_ids'], 'refB')
        self.assertEqual(row5['both_hit_ids'], '')

def airr(rows):
    """
    An AIRR rearrangement table: `rows`, plus every field the AIRR schema
    requires that they lack, blank -- what immvoke's AIRR writer produces.
    """
    table = pandas.DataFrame(rows)
    for field in RearrangementSchema.required:
        if field not in table.columns:
            table[field] = ''

    return table


#: sequence_region -> the columns it compares, as annotateSpecificity does.
REGION_COLUMNS = {'CDR3H': ['cdr3_aa'], 'VDJ': ['sequence_aa'],
                  'CDRs': ['cdr1_aa', 'cdr2_aa', 'cdr3_aa']}


def bruteForceSpecificity(queries, reference, sequence_region='CDR3H',
                          partition_by=None, identity_threshold=None,
                          distance_metric='hamming'):
    """
    annotateSpecificity's single-chain result table, by brute force: the
    same heavy chain filter, columns, gates and threshold, every pair
    compared directly (see BruteForce.py).
    """
    def heavyChains(rows):
        table = pandas.DataFrame(rows)
        if 'locus' in table.columns:
            table = table[table['locus'] == 'IGH']
        return table.reset_index(drop=True)

    queries, reference = heavyChains(queries), heavyChains(reference)
    gene_cols = {'V': ['v_call'], 'VJ': ['v_call', 'j_call'],
                 'VJL': ['v_call', 'j_call']}.get(partition_by)
    exact_cols = None
    if partition_by == 'VJL':
        queries = queries.assign(cdr3_length=queries['cdr3_aa'].str.len())
        reference = reference.assign(cdr3_length=reference['cdr3_aa'].str.len())
        exact_cols = ['cdr3_length']
    exact = identity_threshold is None or identity_threshold >= 1.0

    return bruteForceAnnotation(
        queries, reference, REGION_COLUMNS[sequence_region],
        0 if exact else 1 - identity_threshold,
        distance_type='hamming' if exact else distance_metric, id_col='sequence_id',
        exact_cols=exact_cols, threshold_is_fraction=not exact, gene_cols=gene_cols)


def heavy(sequence_id, cdr3_aa, v_call='IGHV1', j_call='IGHJ1', **extra):
    return dict(sequence_id=sequence_id, locus='IGH', cdr3_aa=cdr3_aa,
                v_call=v_call, j_call=j_call, **extra)


class TestAnnotateSpecificity(unittest.TestCase):
    """
    Tests for the high level entry point, step by step through its flow:
    paired?, heavy chain, sequence_region, partition_by,
    identity_threshold/distance_metric, symscan, the outdir bundle
    """

    def _hits(self, queries, reference, **kwargs):
        result = annotateSpecificity(frame(queries), airr(reference), **kwargs)
        return result['hit_ids'].tolist()

    def test_default_is_an_exact_cdr3_match_reporting_every_query_row(self):
        result = annotateSpecificity(
            frame([heavy('q1', 'CARDYW'), heavy('q2', 'CARDYY')]),
            airr([heavy('ref1', 'CARDYW', v_call='IGHV3')]))

        self.assertEqual(result['hit_ids'].tolist(), ['ref1', ''])
        self.assertEqual(result['n_hits_total'].tolist(), [1, 0])
        self.assertEqual(result.columns[:5].tolist(),
                         ['query_row', 'n_hits_total', 'min_mismatches',
                          'is_unique_hit', 'hit_ids'])

    def test_only_heavy_chain_rows_are_compared(self):
        light = dict(sequence_id='q2', locus='IGK', cdr3_aa='CQQYNS',
                     v_call='IGKV1', j_call='IGKJ1')
        reference = [heavy('ref1', 'CARDYW'), dict(light, sequence_id='ref2')]

        result = annotateSpecificity(frame([heavy('q1', 'CARDYW'), light]),
                                     airr(reference))

        self.assertEqual(result['hit_ids'].tolist(), ['ref1'])

    def test_sequence_region_picks_the_compared_columns(self):
        queries = [heavy('q1', 'CARDYW', sequence_aa='EVQLCARDYW',
                         cdr1_aa='GFT', cdr2_aa='ISY')]
        reference = [heavy('ref1', 'CARDYW', sequence_aa='QVQLCARDYW',
                           cdr1_aa='GFT', cdr2_aa='ISS')]

        self.assertEqual(self._hits(queries, reference, sequence_region='CDR3H'),
                         ['ref1'])
        self.assertEqual(self._hits(queries, reference, sequence_region='VDJ'), [''])
        self.assertEqual(self._hits(queries, reference, sequence_region='CDRs'), [''])

    def test_cdrs_pools_the_threshold_across_cdr1_to_cdr3(self):
        """12 pooled residues at 0.85 identity -> floor(1.8) == 1 substitution,
        here in CDR2 -- no single column is compared on its own."""
        queries = [heavy('q1', 'CARDYW', cdr1_aa='GFT', cdr2_aa='ISY')]
        reference = [heavy('ref1', 'CARDYW', cdr1_aa='GFT', cdr2_aa='ISS'),
                     heavy('ref2', 'CARDYQ', cdr1_aa='GFA', cdr2_aa='ISY')]

        self.assertEqual(self._hits(queries, reference, sequence_region='CDRs',
                                    identity_threshold=0.85), ['ref1'])

    def test_partition_by_gates_on_genes_and_cdr3_length(self):
        queries = [heavy('q1', 'CARDYW', sequence_aa='EVQLCARDYW')]
        reference = [
            heavy('sameVJ', 'CARDDYW', sequence_aa='EVQLCARDDYW'),  # 1 insertion
            heavy('otherJ', 'CARDYW', j_call='IGHJ4', sequence_aa='EVQLCARDYW'),
            heavy('otherV', 'CARDYW', v_call='IGHV3', sequence_aa='EVQLCARDYW'),
        ]
        loose = dict(sequence_region='VDJ', identity_threshold=0.9,
                     distance_metric='levenshtein')

        def hits(partition_by):
            ids = self._hits(queries, reference, partition_by=partition_by, **loose)
            return set(ids[0].split(';')) - {''}

        self.assertEqual(hits(None), {'sameVJ', 'otherJ', 'otherV'})
        self.assertEqual(hits('V'), {'sameVJ', 'otherJ'})
        self.assertEqual(hits('VJ'), {'sameVJ'})
        self.assertEqual(hits('VJL'), set())

    def test_identity_threshold_scales_with_query_length(self):
        """length 10 at 0.85 identity -> floor(1.5) == 1 substitution."""
        queries = [heavy('q1', 'AAAAAAAAAA')]
        reference = [heavy('close', 'AAAAAAAAAB'), heavy('far', 'AABBAAAAAB')]

        self.assertEqual(self._hits(queries, reference, identity_threshold=0.85),
                         ['close'])

    def test_levenshtein_tolerates_an_indel_hamming_does_not(self):
        queries = [heavy('q1', 'CARDYWAAAA')]
        reference = [heavy('ref1', 'CARDXYWAAAA')]  # 1 insertion

        self.assertEqual(self._hits(queries, reference, identity_threshold=0.9,
                                    distance_metric='levenshtein'), ['ref1'])
        self.assertEqual(self._hits(queries, reference, identity_threshold=0.9), [''])

    def test_matches_brute_force(self):
        queries = [heavy('q%d' % i, cdr3, cdr1_aa='GFT', cdr2_aa='ISY')
                   for i, cdr3 in enumerate(['CARDYW', 'CARXYW', 'CARDDYW', 'CQ'])]
        reference = [heavy('r%d' % i, cdr3, cdr1_aa='GFT', cdr2_aa='ISS')
                     for i, cdr3 in enumerate(['CARDYW', 'CARQYW', 'CXRDYW', 'CQ'])]
        modes = [{}, {'identity_threshold': 0.8},
                 {'identity_threshold': 0.8, 'distance_metric': 'levenshtein'},
                 {'sequence_region': 'CDRs', 'identity_threshold': 0.9},
                 {'partition_by': 'VJL', 'identity_threshold': 0.8}]
        for mode in modes:
            with self.subTest(**mode):
                result = annotateSpecificity(frame(queries), airr(reference), **mode)
                pandas.testing.assert_frame_equal(
                    result, bruteForceSpecificity(queries, reference, **mode),
                    check_dtype=False)

    def test_x_matches_any_residue(self):
        self.assertEqual(self._hits([heavy('q1', 'CARXYW'), heavy('q2', 'CAKDYW')],
                                    [heavy('ref1', 'CARDYW')]), ['ref1', ''])

    def test_a_threshold_beyond_symscan_is_skipped_with_a_warning(self):
        """20 residues at 0.25 identity -> 15 differences, above symscan's 14:
        symscan is never called with it, and the row is reported."""
        with mock.patch('symscan.get_neighbors_across',
                        side_effect=AssertionError), \
             self.assertLogs('immvoke.Sources.Specificity.Symscan', 'WARNING'):
            hits = self._hits([heavy('q1', 'A' * 20)],
                              [heavy('ref1', 'B' * 15 + 'A' * 5)],
                              identity_threshold=0.25)

        self.assertEqual(hits, [''])

    def test_paired_classifies_each_antibody(self):
        rows = [dict(cell_id='ab1', locus='IGH', cdr3_aa='CARDYW', v_call='IGHV1',
                     j_call='IGHJ1'),
                dict(cell_id='ab1', locus='IGK', cdr3_aa='CQQYNS', v_call='IGKV1',
                     j_call='IGKJ1')]
        reference = [dict(row, cell_id='refA') for row in rows]

        result = annotateSpecificity(frame(rows), airr(reference), paired=True)
        self.assertEqual(result['classification'].tolist(), ['both'])
        self.assertEqual(result['both_hit_ids'].tolist(), ['refA'])

        with self.assertLogs(ANNOTATE, 'WARNING'):
            annotateSpecificity(frame(rows), airr(reference), paired=True,
                                identity_threshold=0.9)

    def test_reads_paths_and_writes_the_bundle(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            frame([heavy('q1', 'CARDYW'), heavy('q2', 'CQQQQQ')]).to_csv(
                tmp / 'q.tsv', sep='\t', index=False)
            airr([heavy('ref1', 'CARDYW')]).to_csv(tmp / 'r.tsv', sep='\t',
                                                    index=False)

            result = annotateSpecificity(str(tmp / 'q.tsv'), tmp / 'r.tsv',
                                         partition_by='VJ', outdir=tmp / 'sub' / 'run')

            run = tmp / 'sub' / 'run'
            self.assertEqual(sorted(p.name for p in run.iterdir()),
                             [ANNOTATION_FILE, ANNOTATION_METADATA, PAIRS_FILE])
            written = pandas.read_csv(run / ANNOTATION_FILE, sep='\t',
                                      keep_default_na=False)
            self.assertEqual(written['hit_ids'].tolist(), result['hit_ids'].tolist())
            metadata = yaml.safe_load((run / ANNOTATION_METADATA).read_text())
        self.assertEqual(metadata['query'], {'path': str(tmp / 'q.tsv'), 'rows': 2})
        self.assertEqual(metadata['reference'], {'path': str(tmp / 'r.tsv'), 'rows': 1})
        self.assertEqual(metadata['options']['partition_by'], 'VJ')
        self.assertEqual(metadata['results'], {
            'query_rows_compared': 2, 'query_rows_with_hits': 1, 'pairs': 1,
            'skipped_query_rows': 0})
        self.assertNotIn('command', metadata)  # called from Python, not the commandline

    def test_metadata_names_the_reference_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            database = Path(tmp) / 'databases' / 'iedb'
            database.mkdir(parents=True)
            airr([heavy('ref1', 'CARDYW')]).to_csv(database / 'bcr.tsv', sep='\t',
                                                    index=False)
            (database / 'download_metadata.yml').write_text(yaml.safe_dump(
                {'source': 'iedb', 'runs': [{'finished': '2026-10-14T17:52:12Z'}]}))

            annotateSpecificity(frame([heavy('q1', 'CARDYW')]), database / 'bcr.tsv',
                                outdir=Path(tmp) / 'run')

            metadata = yaml.safe_load((Path(tmp) / 'run' / ANNOTATION_METADATA).read_text())
        self.assertEqual(metadata['query']['path'], None)
        self.assertEqual(metadata['reference']['source'], 'iedb')
        self.assertEqual(metadata['reference']['downloaded'], '2026-10-14T17:52:12Z')

    def test_metadata_counts_the_skipped_rows(self):
        with tempfile.TemporaryDirectory() as tmp, \
             self.assertLogs('immvoke.Sources.Specificity.Symscan', 'WARNING'):
            annotateSpecificity(frame([heavy('q1', 'A' * 20), heavy('q2', 'CARDYW')]),
                                airr([heavy('ref1', 'CARDYW'),
                                      heavy('ref2', 'B' * 15 + 'A' * 5)]),
                                identity_threshold=0.25,
                                outdir=tmp)

            metadata = yaml.safe_load((Path(tmp) / ANNOTATION_METADATA).read_text())
        self.assertEqual(metadata['results']['skipped_query_rows'], 1)

    def test_metadata_takes_numpy_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            annotateSpecificity(frame([heavy('q1', 'CARDYW')]),
                                airr([heavy('ref1', 'CARDYW')]),
                                identity_threshold=numpy.float64(0.9), outdir=tmp)

            metadata = yaml.safe_load((Path(tmp) / ANNOTATION_METADATA).read_text())
        self.assertEqual(metadata['options']['identity_threshold'], 0.9)

    def test_a_run_without_hits_fills_text_columns_cleanly(self):
        """No hit at all leaves hit_ids and mismatch_positions all missing,
        which pandas reads as float: filling those with '' must not warn
        (a future pandas raises instead)."""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            result = annotateSpecificity(frame([heavy('q1', 'CQQQQQ')]),
                                         airr([heavy('ref1', 'CARDYW')]))

        self.assertEqual(result[['hit_ids', 'mismatch_positions']].values.tolist(),
                         [['', '']])

    def test_pairs_are_written_one_row_per_hit(self):
        queries = [heavy('q1', 'CARDYW'), heavy('q2', 'CQQQQQ'), heavy('q3', 'CARDYY')]
        reference = [heavy('ref1', 'CARDYW'), heavy('ref2', 'CARDYY')]
        with tempfile.TemporaryDirectory() as tmp:
            result = annotateSpecificity(frame(queries), airr(reference),
                                         identity_threshold=0.8, outdir=tmp)

            pairs = pandas.read_csv(Path(tmp) / PAIRS_FILE, sep='\t',
                                    keep_default_na=False)
        self.assertEqual(result['hit_ids'].tolist(), ['ref1;ref2', '', 'ref2;ref1'])
        self.assertEqual(pairs[['query_row', 'hit_id', 'mismatches']].values.tolist(),
                         [[0, 'ref1', 0], [0, 'ref2', 1], [2, 'ref2', 0], [2, 'ref1', 1]])
        self.assertEqual(pairs['identity'].round(4).tolist(),
                         [1.0, 0.8333, 1.0, 0.8333])

    def test_paired_pairs_name_the_matched_chains(self):
        def antibody(cell_id, heavy_cdr3, light_cdr3):
            return [dict(cell_id=cell_id, locus='IGH', cdr3_aa=heavy_cdr3,
                         v_call='IGHV1', j_call='IGHJ1'),
                    dict(cell_id=cell_id, locus='IGK', cdr3_aa=light_cdr3,
                         v_call='IGKV1', j_call='IGKJ1')]
        queries = antibody('ab1', 'CARDYW', 'CQQYNS') + antibody('ab2', 'CAKKKK', 'CQQAAA')
        reference = antibody('refA', 'CARDYW', 'CQQYNS') + antibody('refB', 'CARDYW',
                                                                    'CQQTTT')
        with tempfile.TemporaryDirectory() as tmp:
            annotateSpecificity(frame(queries), airr(reference), paired=True,
                                outdir=tmp)

            pairs = pandas.read_csv(Path(tmp) / PAIRS_FILE, sep='\t')
            metadata = yaml.safe_load((Path(tmp) / ANNOTATION_METADATA).read_text())
        self.assertEqual(pairs.values.tolist(), [['ab1', 'refA', 'both'],
                                                 ['ab1', 'refB', 'heavy']])
        self.assertEqual(metadata['results'], {
            'query_antibodies': 2, 'classification': {'both': 1, 'neither': 1},
            'pairs': 2, 'skipped_query_rows': 0})

    def test_genes_match_by_name_any_allele(self):
        queries = [heavy('q1', 'CARDYW', v_call='IGHV1-2*02', j_call='IGHJ4*02')]
        reference = [
            heavy('sameAllele', 'CARDYW', v_call='IGHV1-2*02', j_call='IGHJ4*02'),
            heavy('otherAllele', 'CARDYW', v_call='IGHV1-2*04', j_call='IGHJ4*02'),
            heavy('geneOnly', 'CARDYW', v_call='IGHV1-2', j_call='IGHJ4'),
            heavy('several', 'CARDYW', v_call='IGHV3-23*01, IGHV1-2*04',
                  j_call='IGHJ4*02'),
            heavy('blankJ', 'CARDYW', v_call='IGHV1-2*04', j_call=''),
            heavy('sameFamily', 'CARDYW', v_call='IGHV1-69*01', j_call='IGHJ4*02'),
            heavy('otherFamily', 'CARDYW', v_call='IGHV3-23*01', j_call='IGHJ4*02'),
        ]

        hits = self._hits(queries, reference, partition_by='VJ')

        # blankJ: a blank call agrees with nothing; sameFamily: another gene
        self.assertEqual(set(hits[0].split(';')),
                         {'sameAllele', 'otherAllele', 'geneOnly', 'several'})

    def test_messy_gene_calls_match_brute_force(self):
        calls = ['IGHV1-2*02', 'IGHV1-2*04', 'IGHV1-2', '', 'IGVL3-10',
                 'IGHV3-23*01, IGHV1-2*02', 'IGHV1-69*01']
        queries = [heavy('q%d' % i, cdr3, v_call=calls[i % len(calls)],
                         j_call=('IGHJ4*02', '', 'IGHJ6')[i % 3])
                   for i, cdr3 in enumerate(['CARDYW', 'CARQYW', 'CARDDYW', 'CAKDYW'] * 3)]
        reference = [heavy('r%d' % i, cdr3, v_call=calls[(i + 2) % len(calls)],
                           j_call=('IGHJ4', 'IGHJ4*02', '')[i % 3])
                     for i, cdr3 in enumerate(['CARDYW', 'CARQYW', 'CARDYY', 'CQ'] * 3)]
        for mode in ({}, {'identity_threshold': 0.8},
                     {'identity_threshold': 0.8, 'distance_metric': 'levenshtein'}):
            with self.subTest(**mode):
                kwargs = dict(partition_by='VJ', **mode)
                result = annotateSpecificity(frame(queries), airr(reference), **kwargs)
                pandas.testing.assert_frame_equal(
                    result, bruteForceSpecificity(queries, reference, **kwargs),
                    check_dtype=False)

    def _pairedCell(self, heavies, lights=None, **kwargs):
        """
        Classify one query cell `ab1` against reference antibody `refA`
        (heavy CDR3 CARDYW, light CDR3 CQQYNS). `heavies`/`lights` are
        (cdr3_aa, extra fields) per query chain.
        """
        def chain(locus, cdr3, v, j, extra):
            return dict(cell_id='ab1', locus=locus, cdr3_aa=cdr3, v_call=v,
                        j_call=j, **extra)
        lights = lights if lights is not None else [('CQQYNS', {})]
        queries = ([chain('IGH', c, 'IGHV1-2*02', 'IGHJ4*02', e) for c, e in heavies]
                   + [chain('IGK', c, 'IGKV1-39*01', 'IGKJ1*01', e) for c, e in lights])
        reference = [chain('IGH', 'CARDYW', 'IGHV1-2*02', 'IGHJ4*02', {}),
                     chain('IGK', 'CQQYNS', 'IGKV1-39*01', 'IGKJ1*01', {})]
        reference = [dict(row, cell_id='refA') for row in reference]

        result = annotateSpecificity(frame(queries), airr(reference), paired=True,
                                     **kwargs)
        return result['classification'].tolist()

    def test_several_heavy_chains_keep_the_most_abundant(self):
        self.assertEqual(self._pairedCell([('CAKQQQ', {'duplicate_count': '5'}),
                                           ('CARDYW', {'duplicate_count': '50'})]),
                         ['both'])
        self.assertEqual(self._pairedCell([('CAKQQQ', {'duplicate_count': '50'}),
                                           ('CARDYW', {'duplicate_count': '5'})]),
                         ['light_only'])

    def test_several_light_chains_keep_the_most_abundant(self):
        heavy = [('CARDYW', {'duplicate_count': '9'})]
        self.assertEqual(self._pairedCell(heavy, [('CQQAAA', {'duplicate_count': '3'}),
                                                  ('CQQYNS', {'duplicate_count': '7'})]),
                         ['both'])
        self.assertEqual(self._pairedCell(heavy, [('CQQAAA', {'duplicate_count': '7'}),
                                                  ('CQQYNS', {'duplicate_count': '3'})]),
                         ['heavy_only'])

    def test_ties_and_missing_counts_keep_the_first_row(self):
        tie = [('CAKQQQ', {'duplicate_count': '5'}), ('CARDYW', {'duplicate_count': '5'})]
        no_column = [('CAKQQQ', {}), ('CARDYW', {})]
        blank = [('CARDYW', {'duplicate_count': ''}), ('CAKQQQ', {'duplicate_count': '1'})]
        self.assertEqual(self._pairedCell(tie), ['light_only'])
        self.assertEqual(self._pairedCell(no_column), ['light_only'])
        # A blank count ranks below any number.
        self.assertEqual(self._pairedCell(blank), ['light_only'])

    def test_the_count_field_can_be_named(self):
        heavies = [('CAKQQQ', {'duplicate_count': '50', 'umi_count': '1'}),
                   ('CARDYW', {'duplicate_count': '5', 'umi_count': '9'})]
        self.assertEqual(self._pairedCell(heavies), ['light_only'])
        self.assertEqual(self._pairedCell(heavies, count_field='umi_count'), ['both'])

    def test_paired_leaves_out_rows_without_a_cell_id(self):
        """Blank cell_ids are not one antibody: each such chain has no
        partner, so it is left out rather than paired with the others."""
        def chain(cell_id, locus, cdr3):
            v, j = ('IGHV1', 'IGHJ1') if locus == 'IGH' else ('IGKV1', 'IGKJ1')
            return dict(cell_id=cell_id, locus=locus, cdr3_aa=cdr3, v_call=v, j_call=j)
        queries = [chain('', 'IGH', 'CARDYW'), chain('', 'IGH', 'CAKKKK'),
                   chain('', 'IGK', 'CQQYNS'), chain('ab1', 'IGH', 'CAKKKK'),
                   chain('ab1', 'IGK', 'CQQYNS')]
        reference = [chain('refA', 'IGH', 'CAKKKK'), chain('refA', 'IGK', 'CQQYNS')]

        with self.assertLogs(ANNOTATE, 'WARNING') as logs:
            result = annotateSpecificity(frame(queries), airr(reference), paired=True)

        self.assertIn('query: 3 rows have no cell_id', logs.output[0])
        self.assertEqual(result[['cell_id', 'classification']].values.tolist(),
                         [['ab1', 'both']])

    def test_paired_compares_genes_by_name(self):
        rows = [dict(cell_id='ab1', locus='IGH', cdr3_aa='CARDYW', v_call='IGHV1-2*02',
                     j_call='IGHJ4*02'),
                dict(cell_id='ab1', locus='IGK', cdr3_aa='CQQYNS', v_call='IGKV1-39*01',
                     j_call='IGKJ1*01')]
        same_gene = [dict(rows[0], cell_id='refA', v_call='IGHV1-2*04'),
                     dict(rows[1], cell_id='refA')]
        other_gene = [dict(rows[0], cell_id='refA', v_call='IGHV1-69*01'),
                      dict(rows[1], cell_id='refA')]

        for reference, classification in ((same_gene, 'both'),
                                          (other_gene, 'light_only')):
            with self.subTest(classification):
                result = annotateSpecificity(frame(rows), airr(reference), paired=True)
                self.assertEqual(result['classification'].tolist(), [classification])
        reference = same_gene

        # A blank J-gene call agrees with nothing, so the light chain cannot match.
        blank_j = [rows[0], dict(rows[1], j_call='')]
        result = annotateSpecificity(frame(blank_j), airr(reference), paired=True)
        self.assertEqual(result['classification'].tolist(), ['heavy_only'])

    def test_invalid_arguments_raise(self):
        queries = [heavy('q1', 'CARDYW', cdr1_aa='GFT', cdr2_aa='ISY')]
        cases = [{'sequence_region': 'CDR4'}, {'partition_by': 'J'},
                 {'distance_metric': 'blosum'},
                 {'identity_threshold': 0}, {'identity_threshold': 1.5},
                 {'sequence_region': 'CDRs', 'identity_threshold': 0.9,
                  'distance_metric': 'levenshtein'},
                 {'sequence_region': 'VDJ'}]  # no sequence_aa column
        for kwargs in cases:
            with self.subTest(**kwargs), self.assertRaises(AnnotationError):
                self._hits(queries, queries, **kwargs)


class TestReferenceContract(unittest.TestCase):
    """
    The reference is always an AIRR rearrangement table, from whichever
    database -- annotate must accept any such table and refuse anything else
    """

    def test_a_non_airr_reference_is_refused(self):
        query = frame([heavy('q1', 'CARDYW', cell_id='ab1')])
        cov_abdab = frame([{'Name': 'ab1', 'VHorVHH': 'EVQL', 'CDRH3': 'CARDYW',
                            'Heavy V Gene': 'IGHV1-2 (Human)'}])
        for paired in (False, True):
            with self.subTest(paired=paired), \
                 self.assertRaisesRegex(AnnotationError,
                                        'reference is not an AIRR rearrangement'):
                annotateSpecificity(query, cov_abdab, paired=paired)

    def test_a_missing_required_field_is_named(self):
        reference = airr([heavy('ref1', 'CARDYW')]).drop(columns='productive')
        with self.assertRaisesRegex(AnnotationError, 'missing required fields '
                                                     r'\(productive\)'):
            annotateSpecificity(frame([heavy('q1', 'CARDYW')]), reference)

    def test_the_query_needs_only_the_columns_it_uses(self):
        query = frame([{'locus': 'IGH', 'cdr3_aa': 'CARDYW'}])
        result = annotateSpecificity(query, airr([heavy('ref1', 'CARDYW')]))
        self.assertEqual(result['hit_ids'].tolist(), ['ref1'])

    def test_a_reference_written_by_immvoke_is_accepted(self):
        """Every immvoke source that writes AIRR goes through writeAirr."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'reference.tsv'
            writeAirr([frame([heavy('ref1', 'CARDYW', sequence_aa='EVQLCARDYW')])],
                      path)

            result = annotateSpecificity(frame([heavy('q1', 'CARDYW')]), path)

        self.assertEqual(result['hit_ids'].tolist(), ['ref1'])

    def test_any_airr_reference_finds_its_own_rows(self):
        """
        References shaped like different databases: IEDB (allele calls, blank
        J-genes, provenance columns), a user's own IgBLAST/Change-O table
        (several calls in one cell, extra AIRR fields), and the bare minimum.
        Every query row is in its reference, so must find itself, and the
        whole table must match brute force.
        """
        cdr3s = ['CARDYW', 'CARQYW', 'CAKDGYSSGWYFDYW', 'CTTDLW']
        references = {
            'iedb': [heavy('%d_1' % i, c, v_call='IGHV1-2*0%d' % (i % 2 + 1),
                           j_call='' if i % 2 else 'IGHJ4*02',
                           epitope_name='spike', immvoke_source='iedb')
                     for i, c in enumerate(cdr3s)],
            'own': [heavy('seq%d' % i, c, v_call='IGHV3-23*01, IGHV3-23D*01',
                          j_call='IGHJ4*02', productive='T', clone_id=str(i))
                    for i, c in enumerate(cdr3s)],
            'minimal': [{'sequence_id': 'r%d' % i, 'cdr3_aa': c}
                        for i, c in enumerate(cdr3s)],
        }
        modes = [{}, {'identity_threshold': 0.8},
                 {'identity_threshold': 0.8, 'distance_metric': 'levenshtein'}]
        for name, rows in references.items():
            reference = airr(rows)
            queries = frame(rows)
            gated = [{'partition_by': 'VJ'}] if 'v_call' in rows[0] else []
            for mode in modes + gated:
                with self.subTest(reference=name, **mode):
                    result = annotateSpecificity(queries, reference, **mode)
                    pandas.testing.assert_frame_equal(
                        result, bruteForceSpecificity(rows, rows, **mode),
                        check_dtype=False)
                    for row, hits in zip(rows, result['hit_ids']):
                        # Under a partition, a blank J-gene call agrees with
                        # nothing -- not even the row's own.
                        if mode.get('partition_by') and not row['j_call']:
                            self.assertNotIn(row['sequence_id'], hits.split(';'))
                        else:
                            self.assertIn(row['sequence_id'], hits.split(';'))

    def test_no_heavy_chain_rows_left_raises(self):
        light = dict(sequence_id='k1', locus='IGK', cdr3_aa='CQQYNS',
                     v_call='IGKV1-39*01', j_call='IGKJ1*01')
        tcr = dict(light, sequence_id='b1', locus='TRB', v_call='TRBV20-1*01',
                   j_call='TRBJ2-7*01')
        for query, reference in (([light], [heavy('ref1', 'CARDYW')]),
                                 ([heavy('q1', 'CARDYW')], [tcr])):
            with self.subTest(query=query[0]['locus'],
                              reference=reference[0]['locus']), \
                 self.assertRaisesRegex(AnnotationError, 'no heavy chain rows'):
                annotateSpecificity(frame(query), airr(reference))


class TestAnnotateCli(unittest.TestCase):
    """
    Tests for `immvoke annotate`, a thin wrapper over annotateSpecificity
    """

    BASE = ['annotate', '--query', 'q.tsv', '--reference', 'r.tsv', '--outdir', 'run']

    def test_annotate_is_a_top_level_subcommand(self):
        with mock.patch('sys.stdout', new_callable=io.StringIO) as out:
            with self.assertRaises(SystemExit):
                getArgParser().parse_args(['annotate', '--help'])

        text = out.getvalue()
        for flag in ('--query', '--reference', '--sequence-region', '--partition-by',
                     '--paired', '--identity-threshold', '--distance-metric',
                     '--count-field', '--threads', '--outdir'):
            self.assertIn(flag, text)
        for gone in ('--backend', '--gene-match', 'compairr', 'python', '-o,', '--out '):
            self.assertNotIn(gone, text.lower())

    def test_defaults(self):
        args = getArgParser().parse_args(self.BASE)

        self.assertEqual((args.sequence_region, args.partition_by, args.paired,
                          args.identity_threshold, args.distance_metric),
                         ('CDR3H', None, False, None, 'hamming'))

    def test_invalid_choices_are_rejected(self):
        for extra in (['--backend', 'symscan'], ['--partition-by', 'J'],
                      ['--sequence-region', 'CDR4'], ['--distance-metric', 'x'],
                      ['--gene-match', 'gene']):
            with self.subTest(extra=extra), \
                 mock.patch('sys.stderr', new_callable=io.StringIO), \
                 self.assertRaises(SystemExit):
                getArgParser().parse_args(self.BASE + extra)

    def test_required_flags_are_enforced(self):
        with mock.patch('sys.stderr', new_callable=io.StringIO):
            with self.assertRaises(SystemExit):
                getArgParser().parse_args(['annotate'])

    def test_handle_annotate_forwards_every_flag(self):
        args = getArgParser().parse_args(self.BASE + [
            '--sequence-region', 'CDRs', '--partition-by', 'VJL',
            '--identity-threshold', '0.9', '--threads', '4',
            '--count-field', 'umi_count'])

        with mock.patch('immvoke.Sources.Specificity.AnnotateCli.annotateSpecificity') \
                as annotate:
            self.assertEqual(handleAnnotate(args), 0)

        annotate.assert_called_once_with(
            query=Path('q.tsv'), reference=Path('r.tsv'), sequence_region='CDRs',
            partition_by='VJL', paired=False, identity_threshold=0.9,
            distance_metric='hamming', threads=4,
            outdir=Path('run'), count_field='umi_count')

    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            frame([heavy('q1', 'CARDYW'), heavy('q2', 'CQQQQQ')]).to_csv(
                tmp / 'q.tsv', sep='\t', index=False)
            airr([heavy('ref1', 'CARQYW')]).to_csv(tmp / 'r.tsv', sep='\t',
                                                    index=False)
            argv = ['annotate', '--query', str(tmp / 'q.tsv'), '--reference',
                    str(tmp / 'r.tsv'), '--outdir', str(tmp / 'run'),
                    '--partition-by', 'VJ', '--identity-threshold', '0.8']
            args = getArgParser().parse_args(argv)

            with mock.patch('sys.argv', ['/venv/bin/immvoke'] + argv):
                self.assertEqual(handleAnnotate(args), 0)

            written = pandas.read_csv(tmp / 'run' / ANNOTATION_FILE, sep='\t',
                                      keep_default_na=False)
            metadata = yaml.safe_load((tmp / 'run' / ANNOTATION_METADATA).read_text())
        self.assertEqual(written['hit_ids'].tolist(), ['ref1', ''])
        self.assertEqual(metadata['command'], ' '.join(['immvoke'] + argv))
