"""
IEDB - B-cell and T-cell receptor sequences with their epitopes and assays

IEDB (https://www.iedb.org/) publishes receptor sequences as a single bulk
ZIP export; its PostgREST API (https://query-api.iedb.org/api/v1) is used
only to look up the antigen organisms of those receptors' epitopes. A table
is not an index into many separately downloadable files: the table itself
*is* the data unit, so each table produces exactly one DataUnit.

IEDB has no HTML search form to scrape, so harvestSchema returns a hand
curated snapshot rather than one built from a live page.

The receptor tables (bcr, tcr) come from IEDB's own AIRR export,
airr_full.zip: one AIRR rearrangement TSV per table, a record per chain
linked by cell_id, built from IEDB's calculated (IMGT standardized) fields
and carrying the epitope and assay context alongside. They are already AIRR,
so normalizeChunk passes them through, adding only provenance columns. Every
receptor species is kept unless --air_species names one (see
RECEPTOR_SPECIES and _speciesCells), and every antigen unless
--antigen_species names its organism (see _antigenEpitopes). Each record
already carries its epitope, assay ids and outcome, so IEDB's separate
B-cell assay tables (bcell_search, bcr_to_bcell) are not offered.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import difflib
import logging
import zipfile
from datetime import UTC
from pathlib import Path

import pandas

# Immvoke imports
from immvoke.Convert import coerceAirrTypes
from immvoke.Exceptions import IedbParseError, SchemaError
from immvoke.Http import hashFile
from immvoke.Sources.Base import DataUnit, DownloadResult, SourceBase
from immvoke.Sources.Specificity.Paginate import pageByRange
from immvoke.Sources.Specificity.RowHash import rowHashes

log = logging.getLogger(__name__)

#: Endpoints.
API_BASE = 'https://query-api.iedb.org/api/v1'
BULK_URL = 'https://www.iedb.org/downloader.php?file_name=doc/airr_full.zip'
#: Where the bulk ZIP is cached under the output directory.
BULK_ZIP = 'airr_full.zip'

#: Rows requested per page from the PostgREST API.
API_PAGE_SIZE = 2000

#: table -> its member of the AIRR bulk ZIP.
_TABLE_MEMBER = {
    'bcr': 'bcr_full_airr.tsv',
    'tcr': 'tcr_full_airr.tsv',
}

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

#: --outcome's vocabulary. The AIRR export writes iedb_assay_outcome as
#: positive or negative, in lower case, so --outcome matches it ignoring
#: case.
_RECEPTOR_OUTCOMES = ('Positive', 'Negative')

#: --antigen_species takes any antigen organism name IEDB uses, so the schema
#: carries only these as examples (see Field.open_values).
_ANTIGEN_SPECIES_EXAMPLES = ('Influenza A virus', 'SARS-CoV2',
                             'Human immunodeficiency virus 1', 'Homo sapiens')

#: Epitope ids asked for per epitope_search request when resolving
#: --antigen_species; keeps each request URL a few kilobytes long.
EPITOPE_BATCH = 300


def _baseName(name):
    """
    An IEDB organism name without its bracketed part, lower case:
    'Homo sapiens (human)' -> 'homo sapiens', and a strain such as
    'Influenza A virus (A/Puerto Rico/8/1934(H1N1))' -> 'influenza a virus'.
    """
    return name.split(' (', 1)[0].strip().lower()


def _resolveAntigen(name, records):
    """
    Find the IEDB organism an --antigen_species name stands for.

    IEDB's organism names and ids come as two lists per epitope, not always
    in the same order, so a name is paired with its id only where an epitope
    has a single source organism. The name matches ignoring case, and
    without IEDB's bracketed common name ('Homo sapiens' for 'Homo sapiens
    (human)') as long as that leaves a single organism; the full name always
    wins, so 'Influenza A virus' is the species, not one of its strains.

    Arguments:
      name (str): the --antigen_species value.
      records (list): epitope_search records, with source_organism_names
        and source_organism_iris.

    Returns:
      tuple: (IEDB's own name, its organism id).

    Raises:
      SchemaError: if no organism of the table's epitopes has that name,
        naming the closest ones.
    """
    exact, base = {}, {}
    for record in records:
        names = record.get('source_organism_names') or []
        iris = record.get('source_organism_iris') or []
        if len(names) == 1 and len(iris) == 1:
            exact.setdefault(names[0].lower(), (names[0], iris[0]))
            base.setdefault(_baseName(names[0]), set()).add((names[0], iris[0]))

    wanted = name.strip().lower()
    if wanted in exact:
        return exact[wanted]
    if len({iri for _, iri in base.get(wanted, ())}) == 1:
        return sorted(base[wanted])[0]

    close = difflib.get_close_matches(wanted, list(exact), n=5, cutoff=0.6)
    raise SchemaError(
        "no antigen organism named '%s' among this table's epitopes%s"
        % (name, '; closest IEDB names: %s' % ', '.join(exact[x][0] for x in close)
           if close else ''))


def _antigenEpitopes(client, path, chunksize, name):
    """
    The epitopes (reactivity_ref) of a receptor AIRR table whose antigen
    comes from the organism `name` or one under it.

    The AIRR export records each record's epitope but not its antigen's
    organism, so this asks IEDB's epitope_search for the organisms of the
    table's own epitopes, in batches. Its source_organism_iri_search lists
    each epitope's organism and every organism above it, so 'Influenza A
    virus' takes in every flu strain without immvoke knowing the taxonomy.

    Arguments:
      client (HttpClient): the shared HTTP client.
      path (Path): the downloaded AIRR TSV.
      chunksize (int): rows per chunk while scanning it.
      name (str): the --antigen_species value.

    Returns:
      frozenset: the reactivity_ref values to keep.
    """
    ids = set()
    for chunk in pandas.read_csv(path, sep='\t', dtype=str, na_filter=False,
                                 usecols=['reactivity_ref'], chunksize=chunksize):
        ids.update(chunk['reactivity_ref'])
    numbers = sorted(x.rsplit(':', 1)[-1] for x in ids if x)

    records = []
    for start in range(0, len(numbers), EPITOPE_BATCH):
        params = {'select': 'structure_id,source_organism_names,'
                            'source_organism_iris,source_organism_iri_search',
                  'structure_id': 'in.(%s)'
                                  % ','.join(numbers[start:start + EPITOPE_BATCH])}
        for batch in pageByRange(client, API_BASE + '/epitope_search',
                                 page_size=API_PAGE_SIZE, params=params):
            records.extend(batch)

    own, iri = _resolveAntigen(name, records)
    kept = frozenset('IEDB_EPITOPE:%s' % record['structure_id'] for record in records
                     if iri in (record.get('source_organism_iri_search') or []))
    log.info("%s: antigen '%s' (%s): %d of %d epitopes kept", Path(path).name,
             own, iri, len(kept), len(numbers))

    return kept

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
    description = ('IEDB: curated B-cell and T-cell receptor sequences with '
                   'their epitopes and assays')
    homepage = 'https://www.iedb.org/'
    collections = ('bcr', 'tcr')
    collection_help = {
        'bcr': 'BCR (antibody) receptor sequences, from the AIRR bulk export',
        'tcr': 'TCR receptor sequences, from the same AIRR bulk export',
    }
    airr_collections = frozenset(collections)
    #: The raw mirror (the bulk ZIP and the tables extracted from it) is
    #: deleted once iedb_bcr.tsv and iedb_tcr.tsv are written.
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
        outcome = Field(name='outcome', values=_RECEPTOR_OUTCOMES)
        antigen_species = Field(name='antigen_species',
                                values=_ANTIGEN_SPECIES_EXAMPLES, open_values=True)
        receptor_fields = (air_species, antigen_species, outcome)
        collections = {
            'bcr': Collection(name='bcr', fields=receptor_fields),
            'tcr': Collection(name='tcr', fields=receptor_fields),
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
        return [DataUnit(unit_id=_TABLE_MEMBER[query.collection],
                         collection=query.collection, url=BULK_URL,
                         metadata=dict(query.filters))]

    def fetchUnit(self, unit, outdir, resume=True, progress=True):
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

    def readUnit(self, path, unit, chunksize=50000):
        """
        Open a downloaded table.

        Arguments:
          path (Path): the downloaded AIRR TSV.
          unit (DataUnit): what it is.
          chunksize (int): rows per chunk.

        Returns:
          tuple: (metadata, iterator of raw record chunks). metadata holds
          `kept_cells`, the receptors normalizeChunk keeps (see
          _speciesCells), and `kept_epitopes`, the epitopes it keeps (see
          _antigenEpitopes); each is None to keep everything when its
          filter (--air_species, --antigen_species) was not given. The table
          is read as text, blank cells as '' rather than NaN.
        """
        filters = unit.metadata or {}
        species = filters.get('air_species')
        antigen = filters.get('antigen_species')
        metadata = {
            'kept_cells': (None if not species or species == '*'
                           else _speciesCells(path, chunksize,
                                              _speciesCodes(species))),
            'kept_epitopes': (None if not antigen or antigen == '*'
                              else _antigenEpitopes(self.client, path,
                                                    chunksize, antigen)),
        }

        return metadata, pandas.read_csv(path, sep='\t', dtype=str,
                                         na_filter=False, chunksize=chunksize)

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Normalize a chunk, applying provenance and the filters.

        bcr and tcr are already AIRR rearrangement records, one per chain,
        and pass through unchanged, less the rows the filters drop.

        Arguments:
          metadata (dict): from readUnit: `kept_cells` and `kept_epitopes`.
          chunk (pandas.DataFrame): raw records, IEDB's own columns.
          unit (DataUnit): what they came from.
          offset (int): unused; kept for interface parity with SourceBase.
          report (dict): counters to accumulate into.

        Returns:
          pandas.DataFrame: the normalized chunk, with provenance columns
          added.
        """
        report['rows_in'] += len(chunk)
        frame = chunk
        if metadata['kept_cells'] is not None:
            frame = frame[frame['cell_id'].isin(metadata['kept_cells'])]
        if metadata['kept_epitopes'] is not None:
            frame = frame[frame['reactivity_ref'].isin(metadata['kept_epitopes'])]
        outcome = (unit.metadata or {}).get('outcome')
        if outcome and outcome != '*':
            frame = frame[frame['iedb_assay_outcome'].str.lower() == outcome.lower()]
        report['rows_out'] += len(frame)
        frame = coerceAirrTypes(frame.copy())

        row_hash = rowHashes(frame)

        frame = frame.copy()
        frame['immvoke_source'] = self.name
        frame['immvoke_collection'] = unit.collection
        frame['immvoke_unit_id'] = unit.unit_id
        frame['immvoke_row_hash'] = row_hash

        return frame
