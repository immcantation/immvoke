"""
Unit tests for specificity sources (IEDB) and the Range-header paginator
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import json
import tempfile
import unittest
import zipfile
from io import BytesIO
from pathlib import Path

import pandas

# Immvoke imports
from immvoke.Exceptions import HttpError, IedbParseError
from immvoke.Http import HttpClient
from immvoke.Sources.Base import DataUnit, Query
from immvoke.Sources.Specificity import REGISTRY, getSpecificitySource
from immvoke.Sources.Specificity.Iedb import (
    API_BASE,
    BULK_URL,
    IedbSource,
)
from immvoke.Sources.Specificity.Paginate import pageByRange
from tests.FakeHttp import FakeResponse, FakeSession, rangeHandler

#: The columns of IEDB's AIRR export, bcr_full_airr.tsv, in file order.
IEDB_AIRR_COLUMNS = (
    'cell_id', 'sequence_id', 'iedb_receptor_group', 'reactivity_ref',
    'iedb_epitope_name', 'iedb_assay_outcome', 'iedb_assay_ids',
    'iedb_mhc_alleles', 'productive', 'locus', 'locus_species', 'v_call',
    'd_call', 'j_call', 'sequence', 'sequence_aa', 'cdr1_aa', 'cdr2_aa',
    'cdr3_aa', 'junction_aa', 'rev_comp', 'sequence_alignment',
    'germline_alignment', 'junction', 'cdr3', 'v_cigar', 'd_cigar', 'j_cigar',
)


def iedbAirrRows(*rows):
    """
    Rows shaped like IEDB's AIRR export: every column of IEDB_AIRR_COLUMNS,
    blank unless a row gives it.
    """
    return pandas.DataFrame([{c: row.get(c, '') for c in IEDB_AIRR_COLUMNS}
                             for row in rows])


def makeClient(handler):
    """Build an HttpClient with no politeness delay over a scripted session."""
    return HttpClient(delay=0, backoff=0, session=FakeSession(handler))


def itemRangeHandler(pages):
    """
    Build a handler serving JSON pages via item Range headers, PostgREST style.

    Arguments:
      pages (list): one list of records per page, in the order they are served.
        A request past the last page gets a 416, the way IEDB's API ends
        pagination.

    Returns:
      callable: a handler for FakeSession.
    """
    total = sum(len(page) for page in pages)

    def handler(method, url, headers, index):
        if index >= len(pages):
            return FakeResponse(416, b'')

        page = pages[index]
        start = sum(len(p) for p in pages[:index])
        end = start + len(page) - 1
        body = json.dumps(page).encode('utf-8')

        return FakeResponse(206, body, {
            'Content-Range': 'items %d-%d/%d' % (start, end, total)})

    return handler


class TestPageByRange(unittest.TestCase):
    """
    Tests for the PostgREST Range-header paginator
    """

    def test_pages_until_content_range_total_is_reached(self):
        client = makeClient(itemRangeHandler([[{'id': 1}, {'id': 2}], [{'id': 3}]]))
        records = [r for batch in pageByRange(client, 'https://x/y', page_size=2)
                  for r in batch]

        self.assertEqual(records, [{'id': 1}, {'id': 2}, {'id': 3}])

    def test_a_short_final_page_also_ends_pagination(self):
        """A page shorter than page_size ends iteration even with no total."""
        def handler(method, url, headers, index):
            if index == 0:
                return FakeResponse(206, json.dumps([{'id': 1}]).encode())
            raise AssertionError('paginator should have stopped after one page')

        client = makeClient(handler)
        records = [r for batch in pageByRange(client, 'https://x/y', page_size=5)
                  for r in batch]

        self.assertEqual(records, [{'id': 1}])

    def test_416_past_the_last_row_ends_pagination(self):
        client = makeClient(itemRangeHandler([[{'id': 1}]]))
        records = [r for batch in pageByRange(client, 'https://x/y', page_size=1)
                  for r in batch]

        self.assertEqual(records, [{'id': 1}])

    def test_unexpected_status_raises(self):
        client = makeClient(lambda method, url, headers, index: FakeResponse(404))

        with self.assertRaises(HttpError):
            list(pageByRange(client, 'https://x/y'))


def makeZip(members):
    """
    Build ZIP bytes containing the given filename to text-content mapping.
    """
    buffer = BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        for name, content in members.items():
            archive.writestr(name, content)

    return buffer.getvalue()


class TestIedbSource(unittest.TestCase):
    """
    Tests for the IEDB specificity source
    """

    def test_registered(self):
        self.assertIs(REGISTRY['iedb'], IedbSource)
        client = makeClient(lambda *a: FakeResponse(200))
        self.assertIsInstance(getSpecificitySource('iedb', client), IedbSource)

    def test_unknown_source_raises(self):
        client = makeClient(lambda *a: FakeResponse(200))
        with self.assertRaises(KeyError):
            getSpecificitySource('not-a-real-db', client)

    def test_search_units_returns_one_unit_per_table(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))

        for table in source.collections:
            units = source.searchUnits(Query(collection=table, filters={}))
            self.assertEqual(len(units), 1)
            self.assertEqual(units[0].collection, table)

        bcr = source.searchUnits(Query(collection='bcr', filters={}))[0]
        self.assertEqual(bcr.url, BULK_URL)
        bcell = source.searchUnits(
            Query(collection='bcell', filters={'outcome': 'Positive'}))[0]
        self.assertEqual(bcell.url, API_BASE + '/bcell_search')
        self.assertEqual(bcell.metadata, {'outcome': 'Positive'})

    def test_fetch_unit_extracts_the_right_member_from_the_shared_zip(self):
        """bcr and tcr are both members of one AIRR ZIP, extracted by name."""
        zip_bytes = makeZip({
            'tcr_full_airr.tsv': 'sequence_id\tlocus\nt1\tTRB\n',
            'bcr_full_airr.tsv': 'sequence_id\tlocus\nb1\tIGH\n',
        })
        client = makeClient(rangeHandler(zip_bytes))
        source = IedbSource(client)
        unit = DataUnit(unit_id='bcr_full_airr.tsv', collection='bcr', url=BULK_URL)

        with tempfile.TemporaryDirectory() as outdir:
            result = source.fetchUnit(unit, Path(outdir))
            self.assertEqual(result.path.read_text(), 'sequence_id\tlocus\nb1\tIGH\n')
            # The ZIP itself is cached alongside the extracted table, so a
            # second table sharing it does not re-download it.
            self.assertTrue((Path(outdir) / 'airr_full.zip').exists())

    def test_fetch_unit_raises_when_the_zip_member_is_missing(self):
        zip_bytes = makeZip({'tcr_full_airr.tsv': 'sequence_id\nt1\n'})
        client = makeClient(rangeHandler(zip_bytes))
        source = IedbSource(client)
        unit = DataUnit(unit_id='bcr_full_airr.tsv', collection='bcr', url=BULK_URL)

        with tempfile.TemporaryDirectory() as outdir:
            with self.assertRaises(IedbParseError):
                source.fetchUnit(unit, Path(outdir))

    def test_fetch_unit_pages_the_api_and_writes_json(self):
        client = makeClient(itemRangeHandler([[{'receptor_group_id': 1,
                                               'bcell_id': 2}]]))
        source = IedbSource(client)
        unit = DataUnit(unit_id='bcr_to_bcell.json', collection='bcr_to_bcell',
                       url=API_BASE + '/bcr_to_bcell')

        with tempfile.TemporaryDirectory() as outdir:
            result = source.fetchUnit(unit, Path(outdir))
            written = json.loads(result.path.read_text())
            self.assertEqual(written, [{'receptor_group_id': 1, 'bcell_id': 2}])

    def test_normalize_chunk_adds_provenance_columns(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_to_bcell.json', collection='bcr_to_bcell',
                       url=API_BASE + '/bcr_to_bcell', metadata={})
        chunk = pandas.DataFrame([{'receptor_group_id': 1, 'bcell_id': 2}])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['immvoke_source'].tolist(), ['iedb'])
        self.assertEqual(frame['immvoke_collection'].tolist(), ['bcr_to_bcell'])
        self.assertEqual(frame['immvoke_unit_id'].tolist(), ['bcr_to_bcell.json'])
        self.assertEqual(len(frame['immvoke_row_hash'].iloc[0]), 12)
        self.assertEqual(report['rows_in'], 1)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_applies_the_outcome_filter(self):
        """
        bcell_search has no server side filter for qualitative_measure, so
        --outcome is applied client side during normalization, matching the
        grade exactly: Positive does not take in Positive-High.
        """
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcell_search.json', collection='bcell',
                       url=API_BASE + '/bcell_search',
                       metadata={'outcome': 'Positive'})
        chunk = pandas.DataFrame([
            {'qualitative_measure': 'Positive'},
            {'qualitative_measure': 'Negative'},
            {'qualitative_measure': 'Positive-High'},
        ])
        report = {'rows_in': 0, 'rows_out': 0}

        frame = source.normalizeChunk({}, chunk, unit, 0, report)

        self.assertEqual(frame['qualitative_measure'].tolist(), ['Positive'])
        self.assertEqual(report['rows_in'], 3)
        self.assertEqual(report['rows_out'], 1)

    def test_normalize_chunk_keeps_everything_for_the_wildcard(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcell_search.json', collection='bcell',
                       url=API_BASE + '/bcell_search',
                       metadata={'outcome': '*'})
        chunk = pandas.DataFrame([
            {'qualitative_measure': 'Positive'},
            {'qualitative_measure': 'Negative'},
        ])

        frame = source.normalizeChunk({}, chunk, unit, 0,
                                      {'rows_in': 0, 'rows_out': 0})

        self.assertEqual(len(frame), 2)

    def test_read_unit_bulk_reads_the_airr_tsv_as_text(self):
        """Blank cells stay '' rather than NaN, and nothing turns numeric."""
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_full_airr.tsv', collection='bcr', url=BULK_URL)
        frame_in = iedbAirrRows({'cell_id': '1', 'sequence_id': '1_1', 'locus': 'IGH',
                                 'locus_species': 'NCBITAXON:9606',
                                 'cdr3_aa': 'ARDYW', 'iedb_assay_ids': '32'})

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr_full_airr.tsv'
            frame_in.to_csv(path, sep='\t', index=False)
            metadata, chunks = source.readUnit(path, unit)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata, {'kept_cells': None, 'kept_epitopes': None})
        self.assertEqual(list(frame.columns), list(IEDB_AIRR_COLUMNS))
        self.assertEqual(frame['cdr3_aa'].tolist(), ['ARDYW'])
        self.assertEqual(frame['iedb_assay_ids'].tolist(), ['32'])
        self.assertEqual(frame['j_call'].tolist(), [''])

    def test_read_unit_api_json(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = DataUnit(unit_id='bcr_to_bcell.json', collection='bcr_to_bcell',
                       url=API_BASE + '/bcr_to_bcell')

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr_to_bcell.json'
            path.write_text(json.dumps([{'receptor_group_id': 1, 'bcell_id': 2}]))
            metadata, chunks = source.readUnit(path, unit)
            frame = pandas.concat(list(chunks), ignore_index=True)

        self.assertEqual(metadata, {})
        self.assertEqual(frame.to_dict('records'),
                         [{'receptor_group_id': 1, 'bcell_id': 2}])


class TestIedbAirrExport(unittest.TestCase):
    """
    Tests for the receptor tables, read from IEDB's own AIRR export
    """

    def test_bcr_and_tcr_are_airr_collections(self):
        self.assertEqual(IedbSource.airr_collections, frozenset(['bcr', 'tcr']))
        self.assertEqual(IedbSource.airr_collections & set(IedbSource.collections),
                         frozenset(['bcr', 'tcr']))

    def test_bcr_comes_from_the_airr_export(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        unit = source.searchUnits(Query(collection='bcr', filters={}))[0]

        self.assertEqual(unit.unit_id, 'bcr_full_airr.tsv')
        self.assertTrue(unit.url.endswith('airr_full.zip'))

    def _normalize(self, frame, air_species=None, outcome=None,
                   antigen_species=None, handler=None):
        """Read and normalize `frame` the way `immvoke specificity iedb
        download bcr` does, from a downloaded TSV, with --air_species,
        --outcome and --antigen_species if given; `handler` serves IEDB's
        API."""
        source = IedbSource(makeClient(handler or (lambda *a: FakeResponse(200))))
        filters = {k: v for k, v in (('air_species', air_species),
                                     ('outcome', outcome),
                                     ('antigen_species', antigen_species))
                   if v is not None}
        unit = DataUnit(unit_id='bcr_full_airr.tsv', collection='bcr', url=BULK_URL,
                        metadata=filters)
        report = {'rows_in': 0, 'rows_out': 0}

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr_full_airr.tsv'
            frame.to_csv(path, sep='\t', index=False)
            metadata, chunks = source.readUnit(path, unit, chunksize=1)
            normalized = [source.normalizeChunk(metadata, chunk, unit, 0, report)
                          for chunk in chunks]

        return pandas.concat(normalized, ignore_index=True), report

    def _antibody(self):
        """One antibody, both chains, as IEDB's AIRR export writes it."""
        return iedbAirrRows(
            {'cell_id': '1284_2134780p', 'sequence_id': '1284_1', 'locus': 'IGH',
             'locus_species': 'NCBITAXON:9606',
             'productive': 'T', 'v_call': 'IGHV5-9*01, IGHV5-9*03',
             'j_call': 'IGHJ4*02', 'cdr3_aa': 'TRLGDYGYAYTMDY', 'iedb_epitope_name': 'G58, L59',
             'iedb_assay_outcome': 'positive'},
            {'cell_id': '1284_2134780p', 'sequence_id': '1284_2', 'locus': 'IGK',
             'productive': 'T', 'v_call': 'IGKV3-2*01', 'j_call': 'IGKJ2*01',
             'cdr3_aa': 'QQSKEVPLT',
             'iedb_epitope_name': 'G58, L59', 'iedb_assay_outcome': 'positive'})

    def _mixedSpecies(self):
        """Receptors of several species, one of a BALB/c mouse."""
        human, mouse, construct = 'NCBITAXON:9606', 'NCBITAXON:10090', 'NCBITAXON:32630'
        return iedbAirrRows(
            {'cell_id': 'human', 'sequence_id': '1_1', 'locus_species': human},
            {'cell_id': 'human', 'sequence_id': '1_2', 'locus_species': human},
            {'cell_id': 'humanBlank', 'sequence_id': '2_1', 'locus_species': human},
            {'cell_id': 'humanBlank', 'sequence_id': '2_2', 'locus_species': ''},
            {'cell_id': 'mouse', 'sequence_id': '3_1', 'locus_species': mouse},
            {'cell_id': 'mixed', 'sequence_id': '4_1', 'locus_species': human},
            {'cell_id': 'mixed', 'sequence_id': '4_2', 'locus_species': construct},
            {'cell_id': 'blank', 'sequence_id': '5_1', 'locus_species': ''},
            {'cell_id': 'balbc', 'sequence_id': '6_1',
             'locus_species': 'NCBITAXON:10000000'})

    def test_every_species_is_kept_without_air_species(self):
        kept, report = self._normalize(self._mixedSpecies())

        self.assertEqual(len(kept), 9)
        self.assertEqual((report['rows_in'], report['rows_out']), (9, 9))

    def test_air_species_keeps_only_that_species(self):
        """Per receptor: a chain with no species recorded rides along with
        one of the asked species; any chain of another species drops the
        whole receptor."""
        kept, report = self._normalize(self._mixedSpecies(), 'Homo sapiens')

        self.assertEqual(kept['sequence_id'].tolist(), ['1_1', '1_2', '2_1', '2_2'])
        self.assertEqual((report['rows_in'], report['rows_out']), (9, 4))

    def test_air_species_includes_its_strains(self):
        kept, _ = self._normalize(self._mixedSpecies(), 'Mus musculus')
        self.assertEqual(kept['sequence_id'].tolist(), ['3_1', '6_1'])

        kept, _ = self._normalize(self._mixedSpecies(), 'Mus musculus BALB/c')
        self.assertEqual(kept['sequence_id'].tolist(), ['6_1'])

    def test_outcome_keeps_only_that_outcome(self):
        """IEDB's AIRR export writes the outcome in lower case; --outcome
        matches it either way, and combines with --air_species."""
        frame = iedbAirrRows(
            {'cell_id': 'a', 'sequence_id': '1_1', 'iedb_assay_outcome': 'positive',
             'locus_species': 'NCBITAXON:9606'},
            {'cell_id': 'b', 'sequence_id': '2_1', 'iedb_assay_outcome': 'negative',
             'locus_species': 'NCBITAXON:9606'},
            {'cell_id': 'c', 'sequence_id': '3_1', 'iedb_assay_outcome': 'positive',
             'locus_species': 'NCBITAXON:10090'})

        kept, _ = self._normalize(frame, outcome='Negative')
        self.assertEqual(kept['sequence_id'].tolist(), ['2_1'])

        kept, report = self._normalize(frame, air_species='Homo sapiens',
                                       outcome='Positive')
        self.assertEqual(kept['sequence_id'].tolist(), ['1_1'])
        self.assertEqual((report['rows_in'], report['rows_out']), (3, 1))

    def test_air_species_is_a_bcr_and_tcr_filter(self):
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))

        for table in ('bcr', 'tcr'):
            query = source.validateQuery(table, {'air_species': 'Homo sapiens'})
            self.assertEqual(query.filters, {'air_species': 'Homo sapiens',
                                             'antigen_species': '*',
                                             'outcome': '*'})
        with self.assertRaises(Exception):
            source.validateQuery('bcr', {'air_species': 'human'})

    #: epitope_search records: a flu strain, the flu species itself, a human
    #: self antigen, and one epitope with two organisms whose name and id
    #: lists are not in the same order.
    EPITOPES = [
        {'structure_id': 1,
         'source_organism_names': ['Influenza A virus (A/Puerto Rico/8/1934(H1N1))'],
         'source_organism_iris': ['NCBITaxon:211044'],
         'source_organism_iri_search': ['NCBITaxon:1', 'NCBITaxon:11320',
                                        'NCBITaxon:211044']},
        {'structure_id': 2,
         'source_organism_names': ['Influenza A virus'],
         'source_organism_iris': ['NCBITaxon:11320'],
         'source_organism_iri_search': ['NCBITaxon:1', 'NCBITaxon:11320']},
        {'structure_id': 3,
         'source_organism_names': ['Homo sapiens (human)'],
         'source_organism_iris': ['NCBITaxon:9606'],
         'source_organism_iri_search': ['NCBITaxon:1', 'NCBITaxon:9606']},
        {'structure_id': 4,
         'source_organism_names': ['Borreliella burgdorferi B31',
                                   'Borreliella burgdorferi'],
         'source_organism_iris': ['NCBITaxon:139', 'NCBITaxon:224326'],
         'source_organism_iri_search': ['NCBITaxon:1', 'NCBITaxon:139',
                                        'NCBITaxon:224326']},
    ]

    def test_antigen_name_resolves_to_iedb_organism(self):
        from immvoke.Exceptions import SchemaError
        from immvoke.Sources.Specificity.Iedb import _resolveAntigen

        # The full name wins over a strain sharing it before the brackets.
        self.assertEqual(_resolveAntigen('influenza A VIRUS', self.EPITOPES),
                         ('Influenza A virus', 'NCBITaxon:11320'))
        # Without the bracketed common name, when that leaves one organism.
        self.assertEqual(_resolveAntigen('Homo sapiens', self.EPITOPES),
                         ('Homo sapiens (human)', 'NCBITaxon:9606'))
        # Names only from an epitope with two organisms cannot be paired.
        with self.assertRaises(SchemaError):
            _resolveAntigen('Borreliella burgdorferi', self.EPITOPES)
        with self.assertRaisesRegex(SchemaError, 'closest IEDB names: Influenza A virus'):
            _resolveAntigen('Influenza A', self.EPITOPES)

    def test_antigen_species_keeps_that_organism_and_those_under_it(self):
        """'Influenza A virus' keeps the species' own epitope and its
        strain's, through IEDB's source_organism_iri_search."""
        frame = iedbAirrRows(*[
            {'cell_id': str(n), 'sequence_id': '%d_1' % n,
             'reactivity_ref': 'IEDB_EPITOPE:%d' % n} for n in (1, 2, 3, 4)])

        kept, report = self._normalize(frame, antigen_species='Influenza A virus',
                                       handler=itemRangeHandler([self.EPITOPES]))

        self.assertEqual(kept['sequence_id'].tolist(), ['1_1', '2_1'])
        self.assertEqual((report['rows_in'], report['rows_out']), (4, 2))

    def test_antigen_species_accepts_any_name(self):
        """The schema lets any name through; IEDB is asked when fetching."""
        source = IedbSource(makeClient(lambda *a: FakeResponse(200)))
        query = source.validateQuery('tcr', {'antigen_species': 'Yellow fever virus'})
        self.assertEqual(query.filters['antigen_species'], 'Yellow fever virus')

    def test_every_species_name_has_codes(self):
        from immvoke.Sources.Specificity.Iedb import (
            RECEPTOR_SPECIES,
            _speciesCodes,
        )

        for name, parent in RECEPTOR_SPECIES.values():
            self.assertTrue(_speciesCodes(name), name)
            self.assertTrue(parent is None or parent in RECEPTOR_SPECIES, name)

    def test_records_pass_through_with_provenance(self):
        frame, report = self._normalize(self._antibody())

        self.assertEqual(list(frame.columns),
                         list(IEDB_AIRR_COLUMNS) + ['immvoke_source',
                                                    'immvoke_collection',
                                                    'immvoke_unit_id',
                                                    'immvoke_row_hash'])
        self.assertEqual(frame['sequence_id'].tolist(), ['1284_1', '1284_2'])
        self.assertEqual(frame['cell_id'].tolist(), ['1284_2134780p'] * 2)
        self.assertEqual(frame['v_call'].tolist(),
                         ['IGHV5-9*01, IGHV5-9*03', 'IGKV3-2*01'])
        self.assertEqual(frame['immvoke_source'].tolist(), ['iedb'] * 2)
        self.assertEqual((report['rows_in'], report['rows_out']), (2, 2))

    def test_writes_as_a_header_and_row_valid_rearrangement(self):
        """What `immvoke specificity iedb download bcr` writes is valid AIRR, as is."""
        from immvoke import Convert

        frame, _ = self._normalize(self._antibody())

        with tempfile.TemporaryDirectory() as outdir:
            validation = Convert.writeAirr(iter([frame]), Path(outdir) / 'bcr.tsv')

        self.assertTrue(validation['header_valid'])
        self.assertEqual(validation['rows_invalid'], 0)

    def test_written_file_keeps_iedb_columns_and_every_required_field(self):
        """
        The way `immvoke specificity iedb download bcr` writes it (strict_airr
        off): IEDB's own columns survive, and a required field the export
        lacks is still written, blank -- so the file is always a full AIRR
        table.
        """
        from airr.schema import RearrangementSchema

        from immvoke import Convert

        frame, _ = self._normalize(self._antibody().drop(columns=['d_cigar',
                                                                  'rev_comp']))

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr.tsv'
            validation = Convert.writeAirr(iter([frame]), path, strict=False)
            written = pandas.read_csv(path, sep='\t', dtype=str, na_filter=False)

        self.assertTrue(validation['header_valid'])
        self.assertEqual([f for f in RearrangementSchema.required
                          if f not in written.columns], [])
        self.assertEqual(written['d_cigar'].tolist(), ['', ''])
        for column in ('iedb_epitope_name', 'iedb_assay_outcome', 'locus_species',
                       'immvoke_source'):
            self.assertIn(column, written.columns)
        self.assertEqual(written['iedb_epitope_name'].tolist(), ['G58, L59'] * 2)

    def test_annotate_accepts_it_as_the_reference(self):
        from immvoke import Convert
        from immvoke.Sources.Specificity.Annotate import annotateSpecificity

        frame, _ = self._normalize(self._antibody())

        with tempfile.TemporaryDirectory() as outdir:
            path = Path(outdir) / 'bcr.tsv'
            Convert.writeAirr(iter([frame]), path)
            single = annotateSpecificity(frame, path)
            paired = annotateSpecificity(frame, path, paired=True)

        self.assertEqual(single['hit_ids'].tolist(), ['1284_1'])
        self.assertEqual(paired['classification'].tolist(), ['both'])

    def test_download_removes_the_raw_mirror(self):
        """
        `immvoke specificity iedb download` keeps only the normalized tables:
        the bulk ZIP and the tables extracted from it are deleted afterwards.
        """
        from argparse import Namespace
        from unittest import mock

        from immvoke import Cli

        zip_bytes = makeZip({
            'bcr_full_airr.tsv': self._antibody().to_csv(sep='\t', index=False),
            'tcr_full_airr.tsv': 'sequence_id\tlocus\nt1\tTRB\n',
        })
        args = Namespace(db='iedb', table='bcr', limit=None, dry_run=False,
                         no_resume=False, strict_airr=False)

        with tempfile.TemporaryDirectory() as outdir:
            args.outdir = Path(outdir)
            with mock.patch.object(Cli, 'makeClient',
                                   return_value=makeClient(rangeHandler(zip_bytes))):
                Cli.handleSpecificityDownload(args)

            self.assertFalse((Path(outdir) / 'iedb' / 'raw').exists())
            self.assertTrue((Path(outdir) / 'iedb' / 'bcr.tsv').exists())
            self.assertTrue((Path(outdir) / 'iedb' / 'download_metadata.yml').exists())
