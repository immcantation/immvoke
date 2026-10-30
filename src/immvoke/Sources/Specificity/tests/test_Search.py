"""
Unit tests for `immvoke specificity <db> search` and the filter flags' help
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import json
import os
import tempfile
import time
import unittest
from argparse import ArgumentParser
from pathlib import Path
from unittest import mock

# Immvoke imports
from immvoke.Schema import loadSchema
from immvoke.Sources.Base import DataUnit
from immvoke.Sources.Specificity.Iedb import BULK_URL, CACHE_MAX_AGE, IedbSource
from immvoke.Sources.Specificity.Search import (
    SpecificityHelpFormatter,
    _flow,
    addFilterGroup,
    formatSummary,
    searchTable,
    searchTables,
)
from immvoke.Sources.Specificity.tests.test_Specificity import (
    iedbAirrRows,
    isolateCache,
    makeClient,
    makeZip,
)
from tests.FakeHttp import FakeResponse, rangeHandler


def iedbZip():
    """A bulk ZIP holding a small bcr table: a human and a mouse antibody."""
    bcr = iedbAirrRows(
        {'cell_id': 'h', 'sequence_id': 'h_1', 'locus': 'IGH',
         'locus_species': 'NCBITAXON:9606', 'iedb_receptor_group': 'IEDB_RECEPTOR:1',
         'reactivity_ref': 'IEDB_EPITOPE:7', 'iedb_assay_outcome': 'positive'},
        {'cell_id': 'h', 'sequence_id': 'h_2', 'locus': 'IGK',
         'locus_species': 'NCBITAXON:9606', 'iedb_receptor_group': 'IEDB_RECEPTOR:1',
         'reactivity_ref': 'IEDB_EPITOPE:7', 'iedb_assay_outcome': 'positive'},
        {'cell_id': 'm', 'sequence_id': 'm_1', 'locus': 'IGH',
         'locus_species': 'NCBITAXON:10000000', 'iedb_receptor_group': 'IEDB_RECEPTOR:2',
         'reactivity_ref': 'IEDB_EPITOPE:8', 'iedb_assay_outcome': 'negative'})

    return makeZip({'bcr_full_airr.tsv': bcr.to_csv(sep='\t', index=False),
                    'tcr_full_airr.tsv': 'cell_id\n'})


#: epitope_search's answer for iedbZip's two epitopes.
EPITOPES = [
    {'structure_id': 7, 'source_organism_iri_search': ['NCBITaxon:11320'],
     'structure_type': 'Linear peptide',
     'curated_source_antigens': [{'name': 'hemagglutinin',
                                  'source_organism_name': 'Influenza A virus'}]},
    {'structure_id': 8, 'source_organism_iri_search': ['NCBITaxon:9606'],
     'structure_type': 'Non-peptidic',
     'curated_source_antigens': [{'name': 'Angiogenin',
                                  'source_organism_name': 'Homo sapiens (human)'}]},
]


def iedbHandler(zip_bytes):
    """Serve the bulk ZIP, and EPITOPES to any epitope_search request."""
    serve_zip = rangeHandler(zip_bytes)

    def handler(method, url, headers, index):
        if 'epitope_search' in url:
            return FakeResponse(200, json.dumps(EPITOPES).encode('utf-8'))
        return serve_zip(method, url, headers, index)

    return handler


class CountingHandler:
    """Serve a body like rangeHandler, counting the requests made."""

    def __init__(self, body):
        self.serve = rangeHandler(body)
        self.calls = 0

    def __call__(self, *args):
        self.calls += 1
        return self.serve(*args)


class TestZipCache(unittest.TestCase):
    """
    Tests for the bulk ZIP cache download and search share
    """

    def setUp(self):
        self.cache = isolateCache(self)

    def test_the_zip_is_reused_then_refetched_once_stale(self):
        handler = CountingHandler(iedbZip())
        source = IedbSource(makeClient(handler))
        unit = DataUnit(unit_id='bcr_full_airr.tsv', collection='bcr', url=BULK_URL)

        with tempfile.TemporaryDirectory() as outdir:
            source.fetchUnit(unit, Path(outdir), progress=False)
            fetched = handler.calls
            source.fetchUnit(unit, Path(outdir), progress=False)
            self.assertEqual(handler.calls, fetched)

            stale = time.time() - CACHE_MAX_AGE - 60
            os.utime(self.cache / 'airr_full.zip', (stale, stale))
            source.fetchUnit(unit, Path(outdir), progress=False)
            self.assertGreater(handler.calls, fetched)


class TestFilterHelp(unittest.TestCase):
    """
    Tests for the filter flags' help
    """

    def test_flow_never_splits_a_value(self):
        values = ['Mus musculus', 'Mus musculus BALB/c', 'Homo sapiens', 'Rattus']
        lines = _flow(values, 30)

        self.assertTrue(all(len(line) <= 30 for line in lines))
        self.assertEqual(' '.join(x.strip() for x in lines),
                         'Mus musculus, Mus musculus BALB/c, Homo sapiens, Rattus')

    def test_filters_are_grouped_and_described(self):
        parser = ArgumentParser(formatter_class=SpecificityHelpFormatter)
        addFilterGroup(parser, loadSchema('iedb'), IedbSource, 'bcr')
        with mock.patch.dict(os.environ, {'COLUMNS': '100'}):
            text = parser.format_help()

        self.assertIn('filters:', text)
        self.assertIn('--air_species SPECIES', text)
        self.assertIn('--antigen_species ORGANISM', text)
        self.assertIn('--antigen_structure_type TYPE', text)
        self.assertIn('--assay_type ASSAY', text)
        self.assertIn('--disease DISEASE', text)
        self.assertIn('--paired_only ', text)
        self.assertIn('one of: Positive, Negative', text)
        self.assertIn('Mus musculus BALB/c', text)
        self.assertNotIn('default: None', text)
        args = parser.parse_args(['--air_species', 'Homo sapiens', '--paired_only'])
        self.assertEqual(args.filter_air_species, 'Homo sapiens')
        self.assertIs(args.filter_paired_only, True)
        self.assertIsNone(parser.parse_args([]).filter_paired_only)


class TestSearch(unittest.TestCase):
    """
    Tests for counting what a table holds
    """

    def setUp(self):
        isolateCache(self)
        self.source = IedbSource(makeClient(iedbHandler(iedbZip())))

    def test_counts_rows_distinct_values_and_top_values(self):
        with tempfile.TemporaryDirectory() as workdir:
            summary = searchTable(self.source, 'bcr', {}, Path(workdir))

        self.assertEqual((summary['rows_in'], summary['rows']), (3, 3))
        self.assertEqual(summary['distinct'], {'receptor groups': 2, 'epitopes': 2})
        self.assertEqual(dict(summary['top']['locus']), {'IGH': 2, 'IGK': 1})
        # taxon codes are named, strains by their own name
        self.assertEqual(dict(summary['top']['receptor species']),
                         {'Homo sapiens': 2, 'Mus musculus BALB/c': 1})
        self.assertEqual(dict(summary['top']['antigen species']),
                         {'Influenza A virus': 2, 'Homo sapiens (human)': 1})
        self.assertEqual(dict(summary['top']['antigen structure type']),
                         {'Linear peptide': 2, 'Non-peptidic': 1})

    def test_filters_narrow_the_counts(self):
        with tempfile.TemporaryDirectory() as workdir:
            summary = searchTable(self.source, 'bcr', {'air_species': 'Mus musculus'},
                                  Path(workdir))

        self.assertEqual(summary['filters'], {'air_species': 'Mus musculus'})
        self.assertEqual((summary['rows_in'], summary['rows']), (3, 1))
        self.assertEqual(dict(summary['top']['outcome']), {'negative': 1})

    def test_summary_reads_well(self):
        with tempfile.TemporaryDirectory() as workdir:
            summary = searchTable(self.source, 'bcr', {'outcome': 'Positive'},
                                  Path(workdir))
        text = formatSummary(self.source, summary, top=1)

        self.assertIn('filters: outcome = Positive', text)
        self.assertIn('rows             2 of 3', text)
        self.assertIn('locus (top 1 of 2)', text)

    def test_search_writes_nothing_but_the_counts_asked_for(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'counts.tsv'
            with mock.patch('builtins.print'):
                searchTables(self.source, ['bcr'], {}, out=out)

            self.assertEqual(sorted(p.name for p in Path(tmp).iterdir()), ['counts.tsv'])
            lines = out.read_text().splitlines()

        self.assertEqual(lines[0], 'table\tfield\tvalue\trows')
        self.assertIn('bcr\trows\t\t3', lines)
        self.assertIn('bcr\treceptor species\tHomo sapiens\t2', lines)


if __name__ == '__main__':
    unittest.main()
