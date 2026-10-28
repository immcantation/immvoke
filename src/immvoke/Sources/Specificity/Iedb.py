"""
IEDB - B-cell/antibody receptor sequences and assay annotations

IEDB (https://www.iedb.org/) publishes receptor sequences as a single bulk
ZIP export and assay data through a PostgREST API
(https://query-api.iedb.org/api/v1). A table is not an index into many
separately downloadable files: the table itself *is* the data unit, so each
table produces exactly one DataUnit.

IEDB has no HTML search form to scrape, so harvestSchema returns a hand
curated snapshot rather than one built from a live page.

The receptor tables (bcr, tcr) come from IEDB's own AIRR export,
airr_full.zip: one AIRR rearrangement TSV per table, a record per chain
linked by cell_id, built from IEDB's calculated (IMGT standardized) fields
and carrying the epitope and assay context alongside. They are already AIRR,
so normalizeChunk passes them through, adding only provenance columns. Every
receptor species is kept unless --air_species names one (see
RECEPTOR_SPECIES and _speciesCells).
bcell and bcr_to_bcell carry no sequences -- one is assay results, the other
a join table -- so they keep their own column names.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import hashlib
import json
import logging
import zipfile
from datetime import UTC
from pathlib import Path

import pandas

# Immvoke imports
from immvoke.Convert import coerceAirrTypes
from immvoke.Exceptions import IedbParseError
from immvoke.Http import hashFile
from immvoke.Sources.Base import DataUnit, DownloadResult, SourceBase
from immvoke.Sources.Specificity.Paginate import pageByRange

log = logging.getLogger(__name__)

#: Endpoints.
API_BASE = 'https://query-api.iedb.org/api/v1'
BULK_URL = 'https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip'
#: Where the bulk ZIP is cached under the output directory.
BULK_ZIP = 'airr_full.zip'

#: Rows requested per page from the PostgREST API.
API_PAGE_SIZE = 2000

#: table -> (unit_id, url). bcr and tcr are members of one AIRR bulk ZIP;
#: bcell and bcr_to_bcell are each their own API endpoint.
_TABLE_SOURCE = {
    'bcr': ('bcr_full_airr.tsv', BULK_URL),
    'tcr': ('tcr_full_airr.tsv', BULK_URL),
    'bcell': ('bcell_search.json', API_BASE + '/bcell_search'),
    'bcr_to_bcell': ('bcr_to_bcell.json', API_BASE + '/bcr_to_bcell'),
}

#: Tables served from the shared AIRR bulk ZIP rather than the API.
_BULK_TABLES = frozenset(['bcr', 'tcr'])

#: The receptor species IEDB's AIRR export records in `locus_species`, as
#: NCBI taxon code -> (name, parent code). Names are IEDB's own
#: (receptor_species_names on its bcr_search/tcr_search API) without the
#: common name in brackets, so 'Homo sapiens (human)' is 'Homo sapiens'.
#: Codes from 10000000 up are IEDB's own strain and line codes (its ONTIE
#: terms), not real NCBI taxa; each names its parent species, so asking for
#: 'Mus musculus' also keeps BALB/c, C57BL/6, Xenomouse and the rest.
#: Built by joining the export's codes to those names on
#: iedb_receptor_group (export of 2026-10-20); a code IEDB adds later is
#: kept when no species is asked for, but cannot be asked for by name until
#: it is added here.
RECEPTOR_SPECIES = {
    'NCBITAXON:7797': ('Squalus acanthias', None),
    'NCBITAXON:7801': ('Ginglymostoma cirratum', None),
    'NCBITAXON:9031': ('Gallus gallus', None),
    'NCBITAXON:9534': ('Chlorocebus aethiops', None),
    'NCBITAXON:9541': ('Macaca fascicularis', None),
    'NCBITAXON:9544': ('Macaca mulatta', None),
    'NCBITAXON:9598': ('Pan troglodytes', None),
    'NCBITAXON:9606': ('Homo sapiens', None),
    'NCBITAXON:9823': ('Sus scrofa', None),
    'NCBITAXON:9835': ('Camelid', None),
    'NCBITAXON:9837': ('Camelus bactrianus', None),
    'NCBITAXON:9838': ('Camelus dromedarius', None),
    'NCBITAXON:9844': ('Lama glama', None),
    'NCBITAXON:9913': ('Bos taurus', None),
    'NCBITAXON:9986': ('Oryctolagus cuniculus', None),
    'NCBITAXON:10090': ('Mus musculus', None),
    'NCBITAXON:10114': ('Rattus', None),
    'NCBITAXON:10116': ('Rattus norvegicus', None),
    'NCBITAXON:10117': ('Rattus rattus', None),
    'NCBITAXON:30538': ('Vicugna pacos', None),
    'NCBITAXON:32630': ('synthetic construct', None),
    'NCBITAXON:36176': ('Chiloscyllium plagiosum', None),
    'NCBITAXON:168098': ('Orectolobus maculatus', None),
    'NCBITAXON:419612': ('Camelus ferus', None),
    'NCBITAXON:447135': ('Myodes glareolus', None),
    'NCBITAXON:3122392': ('Nothocricetulus migratorius', None),
    'NCBITAXON:10000000': ('Mus musculus BALB/c', 'NCBITAXON:10090'),
    'NCBITAXON:10000038': ('Mus musculus C3H/He', 'NCBITAXON:10090'),
    'NCBITAXON:10000061': ('Mus musculus FVB', 'NCBITAXON:10090'),
    'NCBITAXON:10000067': ('Mus musculus C57BL/6', 'NCBITAXON:10090'),
    'NCBITAXON:10000206': ('Oryctolagus cuniculus New Zealand White',
                           'NCBITAXON:9986'),
    'NCBITAXON:10000228': ('Mus musculus CBA', 'NCBITAXON:10090'),
    'NCBITAXON:10000236': ('Mus musculus CBA/J', 'NCBITAXON:10090'),
    'NCBITAXON:10000248': ('Mus musculus 129 GIX+', 'NCBITAXON:10090'),
    'NCBITAXON:10000275': ('Mus musculus B10.RIII', 'NCBITAXON:10090'),
    'NCBITAXON:10000617': ('Mus musculus HuMAb (Medarex)', 'NCBITAXON:10090'),
    'NCBITAXON:10000645': ('Mus musculus PrP null', 'NCBITAXON:10090'),
    'NCBITAXON:10000655': ('Mus musculus Swiss Webster', 'NCBITAXON:10090'),
    'NCBITAXON:10000657': ('Mus musculus Xenomouse', 'NCBITAXON:10090'),
    'NCBITAXON:10000667': ('Rattus norvegicus Sprague-Dawley',
                           'NCBITAXON:10116'),
    'NCBITAXON:10000915': ('Mus musculus Kunming', 'NCBITAXON:10090'),
    'NCBITAXON:10001602': ('Mus musculus RBF', 'NCBITAXON:10090'),
    'NCBITAXON:10002257': ('Mus musculus Kymouse', 'NCBITAXON:10090'),
}

#: --air_species's vocabulary: every name in RECEPTOR_SPECIES, sorted.
_AIR_SPECIES = tuple(sorted(name for name, _ in RECEPTOR_SPECIES.values()))

#: qualitative_measure has no server-side filter on bcell_search, so it is
#: applied client-side after fetching; this is its known controlled
#: vocabulary.
_QUALITATIVE_MEASURES = ('Positive', 'Positive-High', 'Positive-Intermediate',
                         'Positive-Low', 'Negative')

def _speciesCodes(name):
    """
    The locus_species codes an --air_species name stands for: the species'
    own code and those of its strains.

    Arguments:
      name (str): a name from RECEPTOR_SPECIES.

    Returns:
      frozenset: the matching codes.
    """
    return frozenset(code for code, (own, parent) in RECEPTOR_SPECIES.items()
                     if own == name
                     or (parent is not None and RECEPTOR_SPECIES[parent][0] == name))


def _speciesCells(path, chunksize, codes):
    """
    The receptors (cell_id) of a receptor AIRR table to keep: those with a
    chain of the wanted species and no chain of another. A chain with no
    species recorded does not count against its receptor, so a heavy chain
    of the wanted species keeps its light chain even when only the heavy one
    names a species. Decided per receptor, not per chain, so a pair is never
    split.

    Arguments:
      path (Path): the downloaded AIRR TSV.
      chunksize (int): rows per chunk while scanning it.
      codes (frozenset): the wanted locus_species codes (see _speciesCodes).

    Returns:
      frozenset: the cell_id values to keep.
    """
    wanted, other = set(), set()
    for chunk in pandas.read_csv(path, sep='\t', dtype=str, na_filter=False,
                                 usecols=['cell_id', 'locus_species'],
                                 chunksize=chunksize):
        species = chunk['locus_species']
        matched = species.isin(codes)
        wanted.update(chunk.loc[matched, 'cell_id'])
        other.update(chunk.loc[~matched & (species != ''), 'cell_id'])
    log.info('%s: %d receptors of the asked species kept, %d with a chain of '
             'another species dropped', Path(path).name, len(wanted - other),
             len(wanted & other))

    return frozenset(wanted - other)


def _rowHash(row):
    """
    Build a short content hash for an annotation row.

    Arguments:
      row (pandas.Series): a row of upstream columns, before provenance
        columns are added.

    Returns:
      str: the first 12 hex characters of a SHA-256 digest.
    """
    key = '|'.join(str(x) for x in row.tolist())

    return hashlib.sha256(key.encode('utf-8')).hexdigest()[:12]


def _findZipMember(archive, filename):
    """
    Locate a member of the AIRR bulk ZIP by its filename.

    Arguments:
      archive (zipfile.ZipFile): the opened archive.
      filename (str): the filename to find, ignoring any directory prefix.

    Returns:
      str: the matching member name.

    Raises:
      IedbParseError: if no member has that filename.
    """
    for member in archive.namelist():
        if Path(member).name == filename:
            return member

    raise IedbParseError(
        "'%s' not found in %s; the AIRR bulk export layout has changed"
        % (filename, BULK_URL))


class IedbSource(SourceBase):
    """
    The IEDB specificity source.
    """

    name = 'iedb'
    description = ('IEDB: curated B-cell/antibody receptor sequences and '
                   'assay annotations')
    homepage = 'https://www.iedb.org/'
    collections = ('bcr', 'tcr', 'bcell', 'bcr_to_bcell')
    collection_help = {
        'bcr': 'BCR (antibody) receptor sequences, from the AIRR bulk export',
        'tcr': 'TCR receptor sequences, from the same AIRR bulk export',
        'bcell': 'B-cell assay records: antigen, epitope and qualitative outcome',
        'bcr_to_bcell': 'join table linking BCR receptor groups to B-cell assay records',
    }
    airr_collections = _BULK_TABLES
    #: The raw mirror (the bulk ZIP and the tables extracted from it) is
    #: deleted once bcr.tsv, tcr.tsv and the rest are written.
    keep_raw = False

    #: TODO: confirm exact license wording against https://www.iedb.org/about
    #: before this is relied on for redistribution terms.
    license = 'Freely available; see https://www.iedb.org/about'
    citation = (
        'Vita R, Mahajan S, Overton JA, Dhanda SK, Martini S, Cantrell JR, '
        'Wheeler DK, Sette A, Peters B. The Immune Epitope Database (IEDB): '
        '2018 update. Nucleic Acids Res. 2019;47(D1):D339-D343. '
        'doi:10.1093/nar/gky1006',
    )

    def harvestSchema(self):
        """
        Build the hand curated snapshot.

        There is no live form to scrape, so this returns the same field
        definitions packaged in schema.yaml, timestamped now. `immvoke
        schema refresh --source iedb` exists mainly so the timestamp and any
        future edits to this method stay reproducible the same way a real
        harvest would.

        Returns:
          SourceSchema: the snapshot.
        """
        from datetime import datetime

        from immvoke.Schema import Collection, Field, SourceSchema
        from immvoke.Version import __version__

        air_species = Field(name='air_species', values=_AIR_SPECIES)
        collections = {
            'bcr': Collection(name='bcr', fields=(air_species,)),
            'tcr': Collection(name='tcr', fields=(air_species,)),
            'bcell': Collection(name='bcell', fields=(
                Field(name='qualitative_measure',
                      values=_QUALITATIVE_MEASURES),
            )),
            'bcr_to_bcell': Collection(name='bcr_to_bcell'),
        }

        return SourceSchema(
            source=self.name,
            harvested=datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ'),
            harvested_by=('immvoke %s (hand curated, not scraped: IEDB '
                          'exposes no HTML form)' % __version__),
            source_urls={'airr_bulk': BULK_URL, 'api': API_BASE},
            collections=collections)

    def searchUnits(self, query):
        """
        Resolve a query to the one DataUnit its table represents.

        A table is not an index of many files; the table itself is the unit,
        so there is always exactly one.

        Arguments:
          query (Query): the validated request.

        Returns:
          list: a single DataUnit.
        """
        unit_id, url = _TABLE_SOURCE[query.collection]

        return [DataUnit(unit_id=unit_id, collection=query.collection, url=url,
                         metadata=dict(query.filters))]

    def fetchUnit(self, unit, outdir, resume=True, progress=True):
        """
        Download one table.

        Arguments:
          unit (DataUnit): what to fetch.
          outdir (Path): the mirror root.
          resume (bool): continue an interrupted transfer if possible.
          progress (bool): show a progress bar.

        Returns:
          DownloadResult: the outcome.
        """
        if unit.collection in _BULK_TABLES:
            return self._fetchFromBulkZip(unit, outdir, resume=resume,
                                          progress=progress)

        return self._fetchFromApi(unit, outdir)

    def _fetchFromBulkZip(self, unit, outdir, resume=True, progress=True):
        """
        Fetch the shared AIRR bulk ZIP and extract one table from it.

        The ZIP is cached at a fixed path under outdir: a rerun with the file
        already present reuses it (HttpClient.fetch skips a download whose
        destination already exists), so fetching bcr and then tcr in the same
        output directory downloads the ZIP only once.

        Arguments:
          unit (DataUnit): the bcr or tcr table.
          outdir (Path): the mirror root.
          resume (bool): continue an interrupted ZIP transfer if possible.
          progress (bool): show a progress bar for the ZIP transfer.

        Returns:
          DownloadResult: the outcome, describing the extracted table file.
        """
        zip_dest = Path(outdir) / BULK_ZIP
        self.client.fetch(unit.url, zip_dest, resume=resume, progress=progress)

        dest = self.resolveOutputPath(unit, outdir)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_dest) as archive:
            member = _findZipMember(archive, unit.unit_id)
            with archive.open(member) as source, open(dest, 'wb') as target:
                target.write(source.read())

        return DownloadResult(unit=unit, path=dest, sha256=hashFile(dest),
                              size_bytes=dest.stat().st_size)

    def _fetchFromApi(self, unit, outdir):
        """
        Page the PostgREST API for one table and write it as JSON.

        Arguments:
          unit (DataUnit): the bcell or bcr_to_bcell table.
          outdir (Path): the mirror root.

        Returns:
          DownloadResult: the outcome.
        """
        records = []
        for batch in pageByRange(self.client, unit.url, page_size=API_PAGE_SIZE):
            records.extend(batch)

        dest = self.resolveOutputPath(unit, outdir)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(records))

        return DownloadResult(unit=unit, path=dest, sha256=hashFile(dest),
                              size_bytes=dest.stat().st_size)

    def readUnit(self, path, unit, chunksize=50000):
        """
        Open a downloaded table.

        Arguments:
          path (Path): the downloaded file.
          unit (DataUnit): what it is.
          chunksize (int): rows per chunk, for the AIRR TSV tables. The
            JSON tables are modest enough to read as a single chunk.

        Returns:
          tuple: (metadata, iterator of raw record chunks). IEDB tables carry
          no per-unit metadata line; for bcr and tcr, metadata holds
          `kept_cells`, the receptors normalizeChunk keeps (see
          _speciesCells), or None to keep every receptor when no
          --air_species was given. bcr and tcr are read as text, blank cells
          as '' rather than NaN.
        """
        if unit.collection in _BULK_TABLES:
            name = (unit.metadata or {}).get('air_species')
            kept = (None if not name or name == '*'
                    else _speciesCells(path, chunksize, _speciesCodes(name)))
            metadata = {'kept_cells': kept}
            return metadata, pandas.read_csv(path, sep='\t', dtype=str,
                                             na_filter=False, chunksize=chunksize)

        with open(path) as handle:
            records = json.load(handle)

        return {}, iter([pandas.DataFrame.from_records(records)])

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Normalize a chunk, applying provenance and any filter the fetch step
        could not apply server-side.

        bcr and tcr are already AIRR rearrangement records, one per chain,
        and pass through unchanged; so do bcell and bcr_to_bcell's own
        columns, filtered where the fetch step could not filter.

        Arguments:
          metadata (dict): from readUnit; for bcr and tcr, `kept_cells`.
          chunk (pandas.DataFrame): raw records, IEDB's own columns.
          unit (DataUnit): what they came from.
          offset (int): unused; kept for interface parity with SourceBase.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the normalized chunk, with provenance columns
          added.
        """
        if unit.collection in _BULK_TABLES:
            report['rows_in'] += len(chunk)
            kept = metadata['kept_cells']
            frame = chunk if kept is None else chunk[chunk['cell_id'].isin(kept)]
            report['rows_out'] += len(frame)
            frame = coerceAirrTypes(frame.copy())
        else:
            frame = self._normalizeAssayChunk(chunk, unit, report)

        # apply(..., axis=1) on an empty frame cannot infer a row shape, so it
        # is skipped rather than trusted to return an empty result.
        if len(frame):
            row_hash = frame.apply(_rowHash, axis=1)
        else:
            row_hash = pandas.Series([], dtype=str, index=frame.index)

        frame = frame.copy()
        frame['immvoke_source'] = self.name
        frame['immvoke_collection'] = unit.collection
        frame['immvoke_unit_id'] = unit.unit_id
        frame['immvoke_row_hash'] = row_hash

        return frame

    def _normalizeAssayChunk(self, chunk, unit, report):
        """
        Pass an assay or join table chunk's columns through unchanged,
        applying any filter the fetch step could not apply server-side.

        Arguments:
          chunk (pandas.DataFrame): raw records, IEDB's own columns.
          unit (DataUnit): what they came from.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the filtered chunk.
        """
        report['rows_in'] += len(chunk)

        frame = self._applyClientFilter(chunk, unit)

        report['rows_out'] += len(frame)

        return frame

    def _applyClientFilter(self, frame, unit):
        """
        Apply qualitative_measure filtering that bcell_search cannot do itself.

        Arguments:
          frame (pandas.DataFrame): raw records for one chunk.
          unit (DataUnit): the table, carrying the resolved filters.

        Returns:
          pandas.DataFrame: the filtered chunk, or the chunk unchanged for
          tables and filters this does not apply to.
        """
        wanted = (unit.metadata or {}).get('qualitative_measure')
        if not wanted or wanted == '*' or 'qualitative_measure' not in frame.columns:
            return frame

        return frame[frame['qualitative_measure'] == wanted]
