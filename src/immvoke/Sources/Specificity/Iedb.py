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
RECEPTOR_SPECIES and _speciesCells), every antigen unless
--antigen_species names its organism (see _antigenEpitopes), and every
epitope structure unless --antigen_structure_type names one (see
_structureEpitopes). Each record
already carries its epitope, assay ids and outcome, so IEDB's separate
B-cell assay tables (bcell_search, bcr_to_bcell) are not offered.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import difflib
import logging
import os
import time
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
#: The bulk ZIP's file name in the cache (see cacheDir).
BULK_ZIP = 'airr_full.zip'
#: Seconds a cached bulk ZIP is reused before it is fetched again. IEDB
#: refreshes the export weekly and its download sends no Last-Modified or
#: ETag to check against, so freshness goes by age.
CACHE_MAX_AGE = 24 * 60 * 60

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

#: --antigen_structure_type's vocabulary: IEDB's epitope structure_type, the
#: same choice as its home page search's Epitope Structure. 'Discontinuous
#: peptide' also keeps 'Discontinuous peptide on multi chain', the way a
#: species keeps its strains.
_ANTIGEN_STRUCTURE_TYPES = ('Linear peptide', 'Discontinuous peptide',
                            'Discontinuous peptide on multi chain', 'Non-peptidic')

#: Epitope ids asked for per epitope_search request (see _epitopeRecords);
#: keeps each request URL a few kilobytes long.
EPITOPE_BATCH = 300

#: The columns added to bcr and tcr from each record's epitope: its
#: antigen, the organism that antigen comes from, and the epitope's
#: structure type (see _epitopeAntigens), placed after iedb_epitope_name.
ANTIGEN_COL = 'iedb_antigen'
ANTIGEN_SPECIES_COL = 'iedb_antigen_source_species'
ANTIGEN_STRUCTURE_COL = 'iedb_antigen_structure_type'
#: Joins an epitope's antigens when IEDB curates several.
ANTIGEN_SEP = ' | '


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


def _epitopeRecords(client, path, chunksize):
    """
    Look up the epitopes of a receptor AIRR table on IEDB's epitope_search.

    The AIRR export records each record's epitope (reactivity_ref) but not
    the epitope's antigen, that antigen's organism or the epitope's
    structure type, so this asks the API about the table's own epitopes, in
    batches of EPITOPE_BATCH.

    Arguments:
      client (HttpClient): the shared HTTP client.
      path (Path): the downloaded AIRR TSV.
      chunksize (int): rows per chunk while scanning it.

    Returns:
      list: epitope_search records: structure_id, the source organism names,
      ids and ids with every organism above them (source_organism_iri_search),
      curated_source_antigens and structure_type.
    """
    ids = set()
    for chunk in pandas.read_csv(path, sep='\t', dtype=str, na_filter=False,
                                 usecols=['reactivity_ref'], chunksize=chunksize):
        ids.update(chunk['reactivity_ref'])
    numbers = sorted(x.rsplit(':', 1)[-1] for x in ids if x)

    records = []
    for start in range(0, len(numbers), EPITOPE_BATCH):
        params = {'select': 'structure_id,source_organism_names,'
                            'source_organism_iris,source_organism_iri_search,'
                            'curated_source_antigens,structure_type',
                  'structure_id': 'in.(%s)'
                                  % ','.join(numbers[start:start + EPITOPE_BATCH])}
        for batch in pageByRange(client, API_BASE + '/epitope_search',
                                 page_size=API_PAGE_SIZE, params=params):
            records.extend(batch)
    log.info('%s: looked up %d epitopes on IEDB', Path(path).name, len(numbers))

    return records


def _antigenEpitopes(records, name):
    """
    The epitopes (reactivity_ref) whose antigen comes from the organism
    `name` or one under it.

    IEDB's source_organism_iri_search lists each epitope's organism and
    every organism above it, so 'Influenza A virus' takes in every flu
    strain without immvoke knowing the taxonomy.

    Arguments:
      records (list): epitope_search records, from _epitopeRecords.
      name (str): the --antigen_species value.

    Returns:
      frozenset: the reactivity_ref values to keep.
    """
    own, iri = _resolveAntigen(name, records)
    kept = frozenset('IEDB_EPITOPE:%s' % record['structure_id'] for record in records
                     if iri in (record.get('source_organism_iri_search') or []))
    log.info("antigen '%s' (%s): %d of %d epitopes kept", own, iri, len(kept),
             len(records))

    return kept


def _structureEpitopes(records, name):
    """
    The epitopes (reactivity_ref) of the structure type `name`, or of one
    under it ('Discontinuous peptide' takes in 'Discontinuous peptide on
    multi chain'). The schema has already checked the name against
    _ANTIGEN_STRUCTURE_TYPES; IEDB's own value is compared ignoring case, in
    case its capitalization drifts.

    Arguments:
      records (list): epitope_search records, from _epitopeRecords.
      name (str): the --antigen_structure_type value.

    Returns:
      frozenset: the reactivity_ref values to keep.
    """
    wanted = name.strip().lower()
    kept = frozenset(
        'IEDB_EPITOPE:%s' % record['structure_id'] for record in records
        if (record.get('structure_type') or '').lower() in (wanted, wanted + ' on multi chain'))
    log.info("antigen structure type '%s': %d of %d epitopes kept", name, len(kept),
             len(records))

    return kept


def _epitopeAntigens(records):
    """
    Each epitope's antigen and that antigen's source species, from IEDB's
    curated_source_antigens, and its structure_type, for the iedb_antigen,
    iedb_antigen_source_species and iedb_antigen_structure_type columns.

    Arguments:
      records (list): epitope_search records, from _epitopeRecords.

    Returns:
      dict: reactivity_ref -> (antigen names, their source species,
      structure type), the first two ' | ' joined in IEDB's order when an
      epitope has several antigens, each blank when IEDB names none.
    """
    antigens = {}
    for record in records:
        curated = record.get('curated_source_antigens') or []
        names = list(dict.fromkeys(x.get('name') or '' for x in curated))
        species = list(dict.fromkeys(x.get('source_organism_name') or ''
                                     for x in curated))
        antigens['IEDB_EPITOPE:%s' % record['structure_id']] = (
            ANTIGEN_SEP.join(x for x in names if x),
            ANTIGEN_SEP.join(x for x in species if x),
            record.get('structure_type') or '')

    return antigens

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


def cacheDir():
    """
    Where IEDB's bulk ZIP is cached between runs: $XDG_CACHE_HOME/immvoke/iedb,
    ~/.cache/immvoke/iedb when that is not set.

    Returns:
      Path: the cache folder (not created here).
    """
    root = os.environ.get('XDG_CACHE_HOME') or Path.home() / '.cache'

    return Path(root) / 'immvoke' / 'iedb'


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
    #: The tables extracted from the bulk ZIP are deleted once iedb_bcr.tsv
    #: and iedb_tcr.tsv are written; the ZIP itself stays in the cache.
    keep_raw = False
    filter_help = {
        'air_species': ('SPECIES', "the receptor's species, from each chain's "
                                   'locus_species; a species includes its '
                                   'strains'),
        'antigen_species': ('ORGANISM', 'the organism its epitope comes from: '
                                        'any name IEDB uses, ignoring case; an '
                                        'organism includes those under it'),
        'antigen_structure_type': ('TYPE', "its epitope's structure, IEDB's "
                                           'structure_type; Discontinuous '
                                           'peptide includes its multi chain '
                                           'kind'),
        'outcome': ('OUTCOME', 'the assay outcome, iedb_assay_outcome'),
    }
    search_columns = dict.fromkeys(collections, (
        ('iedb_receptor_group', 'receptor groups', 'distinct'),
        ('reactivity_ref', 'epitopes', 'distinct'),
        ('locus', 'locus', 'top'),
        ('locus_species', 'receptor species', 'top'),
        (ANTIGEN_SPECIES_COL, 'antigen species', 'top'),
        (ANTIGEN_STRUCTURE_COL, 'antigen structure type', 'top'),
        ('iedb_assay_outcome', 'outcome', 'top'),
    ))

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
        antigen_structure_type = Field(name='antigen_structure_type',
                                       values=_ANTIGEN_STRUCTURE_TYPES)
        receptor_fields = (air_species, antigen_species, antigen_structure_type, outcome)
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

    def labelValue(self, column, value):
        """
        Name a locus_species code as its species (see RECEPTOR_SPECIES), and
        a blank cell as not recorded.

        Arguments:
          column (str): the column.
          value (str): one of its values.

        Returns:
          str: the label to show.
        """
        if value == '':
            return '(not recorded)'
        if column == 'locus_species':
            return RECEPTOR_SPECIES.get(value, (value, None))[0]

        return value

    def fetchUnit(self, unit, outdir, resume=True, progress=True):
        """
        Fetch the shared AIRR bulk ZIP and extract one table from it.

        The ZIP is kept in the user's cache (see cacheDir), not under outdir,
        so download and search share it, and fetching bcr and then tcr
        downloads it only once. It is reused for CACHE_MAX_AGE, then fetched
        again.

        Arguments:
          unit (DataUnit): the bcr or tcr table.
          outdir (Path): the mirror root, where the table is extracted to.
          resume (bool): continue an interrupted ZIP transfer if possible.
          progress (bool): show a progress bar for the ZIP transfer.

        Returns:
          DownloadResult: the outcome, describing the extracted table file.
        """
        zip_dest = cacheDir() / BULK_ZIP
        if zip_dest.exists() and time.time() - zip_dest.stat().st_mtime > CACHE_MAX_AGE:
            log.info('%s is over a day old; fetching the current export', zip_dest)
            zip_dest.unlink()
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
          _antigenEpitopes and _structureEpitopes, both applying when both
          are given), each None to keep everything when its filter
          (--air_species; --antigen_species, --antigen_structure_type) was
          not given; and `antigens`, each epitope's antigen, source species
          and structure type (see _epitopeAntigens). The table is read as
          text, blank cells as '' rather than NaN.
        """
        filters = unit.metadata or {}
        species = filters.get('air_species')
        antigen = filters.get('antigen_species')
        structure = filters.get('antigen_structure_type')
        records = _epitopeRecords(self.client, path, chunksize)
        kept_epitopes = None
        if antigen and antigen != '*':
            kept_epitopes = _antigenEpitopes(records, antigen)
        if structure and structure != '*':
            by_structure = _structureEpitopes(records, structure)
            kept_epitopes = (by_structure if kept_epitopes is None
                             else kept_epitopes & by_structure)
        metadata = {
            'kept_cells': (None if not species or species == '*'
                           else _speciesCells(path, chunksize,
                                              _speciesCodes(species))),
            'kept_epitopes': kept_epitopes,
            'antigens': _epitopeAntigens(records),
        }

        return metadata, pandas.read_csv(path, sep='\t', dtype=str,
                                         na_filter=False, chunksize=chunksize)

    def normalizeChunk(self, metadata, chunk, unit, offset, report):
        """
        Normalize a chunk, applying provenance and the filters.

        bcr and tcr are already AIRR rearrangement records, one per chain,
        and pass through unchanged, less the rows the filters drop, plus
        ANTIGEN_COL, ANTIGEN_SPECIES_COL and ANTIGEN_STRUCTURE_COL after
        iedb_epitope_name. Those
        come from IEDB's API rather than its export, so the row hash leaves
        them out.

        Arguments:
          metadata (dict): from readUnit: `kept_cells`, `kept_epitopes` and
            `antigens`.
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
        found = frame['reactivity_ref'].map(metadata['antigens'])
        at = frame.columns.get_loc('iedb_epitope_name') + 1
        frame.insert(at, ANTIGEN_COL, [x[0] if isinstance(x, tuple) else ''
                                       for x in found])
        frame.insert(at + 1, ANTIGEN_SPECIES_COL, [x[1] if isinstance(x, tuple) else ''
                                                   for x in found])
        frame.insert(at + 2, ANTIGEN_STRUCTURE_COL, [x[2] if isinstance(x, tuple) else ''
                                                     for x in found])
        frame['immvoke_source'] = self.name
        frame['immvoke_collection'] = unit.collection
        frame['immvoke_unit_id'] = unit.unit_id
        frame['immvoke_row_hash'] = row_hash

        return frame
