"""
immvoke commandline interface

Download data from online immune repertoire databases and format it for
Immcantation.
"""

# Info
__author__ = 'Susanna Marquez'

# Imports
import logging
import shutil
import sys
import time
from argparse import ArgumentParser, Namespace
from pathlib import Path

# Immvoke imports
from immvoke import Catalog, Convert, Provenance, Reference
from immvoke.Airrflow import buildSamplesheet
from immvoke.Commandline import CommonHelpFormatter, setupLogging
from immvoke.Exceptions import ImmvokeError
from immvoke.Http import HttpClient
from immvoke.Schema import loadSchema, saveSchema
from immvoke.Sources import ALIASES, REGISTRY, canonicalName, getSource
from immvoke.Sources.Specificity import REGISTRY as SPECIFICITY_REGISTRY
from immvoke.Sources.Specificity import Report, getSpecificitySource
from immvoke.Sources.Specificity.AnnotateCli import addAnnotateParser, handleAnnotate
from immvoke.Sources.Specificity.Provenance import writeSpecificityMetadata
from immvoke.Sources.Specificity.Search import (
    DEFAULT_TOP,
    SpecificityHelpFormatter,
    addFilterGroup,
    searchTables,
)
from immvoke.Version import __date__, __version__

log = logging.getLogger('immvoke')

#: Output formats the download and convert subcommands can produce.
FORMATS = ('raw', 'airr', 'fasta')

#: Pseudo-collection for germline sources: fetch every species immvoke supports
#: into one reference_base. Not offered for OAS (its collections are paired and
#: unpaired, not species) nor for search (two species would merge two hit lists).
ALL_SPECIES = 'all'

#: Above this many values, a filter flag's help lists only a sample instead of
#: everything, and points at `schema show` for the rest. Enumerated fields are
#: categorical (species, disease, ...), so in practice this rarely bites; it
#: exists so one field having many values can't blow up every --help screen.
VALUE_LIST_CAP = 20

#: Every source `immvoke schema` works on: both registries and the aliases.
#: Specificity databases keep a snapshot like any other source, so `schema
#: show` lists their filter values -- which is what a rejected filter value
#: tells the user to run (see Schema.validateFilters).
SCHEMA_SOURCES = sorted(set(REGISTRY) | set(SPECIFICITY_REGISTRY) | set(ALIASES))


def getAnySource(name, client):
    """
    Instantiate a source from whichever registry holds it.

    `immvoke schema` is the one command that spans both: every source keeps
    a snapshot, whether it is mounted under `immvoke <source>` or
    `immvoke specificity <db>`.

    Arguments:
      name (str): the source name or alias.
      client (HttpClient): the shared HTTP client.

    Returns:
      SourceBase: the source.
    """
    if canonicalName(name) in REGISTRY:
        return getSource(name, client)

    return getSpecificitySource(name, client)


def loadSchemaQuietly(name):
    """
    Load a packaged snapshot, returning None instead of raising.

    Parser construction must work on a checkout that has no snapshot yet,
    otherwise `immvoke schema refresh` could never be run to create one.

    Arguments:
      name (str): the source name.

    Returns:
      SourceSchema: the snapshot, or None.
    """
    try:
        return loadSchema(name)
    except Exception:
        return None


def addFilterArgs(parser, schema, source, collection):
    """
    Generate one commandline flag per searchable field.

    Flags come from the stored snapshot, never from a hardcoded list, so the day
    the source adds a field its flag appears with no code change. This reads only
    packaged data and performs no network access: it runs at documentation build
    time as well as at runtime.

    Arguments:
      parser (ArgumentParser): the subparser to add to.
      schema (SourceSchema): the snapshot, or None if none is installed.
      source (str): the source name, for the `schema show` pointer on overflow.
      collection (str): which collection's fields to add.
    """
    if schema is None or collection not in schema.collections:
        return

    for item in schema.getCollection(collection).fields:
        if item.pseudo_values:
            summary = 'filter on whether %s is recorded' % item.name
        elif item.open_values:
            summary = ('any name %s uses, checked when fetching, e.g. %s'
                       % (source, ', '.join(item.values)))
        elif len(item.values) <= VALUE_LIST_CAP:
            summary = '%d values: %s' % (len(item.values), ', '.join(item.values))
        else:
            # Overflow only: today's fields (species, disease, ...) all stay well
            # under the cap. Once `immvoke build` (the interactive command
            # builder, see plan phase 6) exists, point there instead.
            shown = ', '.join(item.values[:VALUE_LIST_CAP])
            summary = ('%d values, e.g. %s, ... run `immvoke schema show '
                       '--source %s --collection %s --field %s` for the full list'
                       % (len(item.values), shown, source, collection, item.name))

        parser.add_argument(item.flag, dest='filter_%s' % item.name,
                            metavar='VALUE', default=None, help=summary)


def collectFilters(args):
    """
    Gather the generated filter flags the user actually supplied.

    Arguments:
      args (Namespace): parsed arguments.

    Returns:
      dict: field name to value.
    """
    return {k[len('filter_'):]: v for k, v in vars(args).items()
            if k.startswith('filter_') and v is not None}


def getArgParser():
    """
    Build the top level argument parser.

    Defined as a function returning an ArgumentParser so that
    sphinxcontrib-autoprogram can document the tool, per Immcantation's
    CONTRIBUTING.md.

    Returns:
      argparse.ArgumentParser: the top level parser.
    """
    parser = ArgumentParser(prog='immvoke', description=__doc__,
                            formatter_class=CommonHelpFormatter)
    # NB: %(prog)s is expanded by argparse, so it must not be part of the string
    # being %-formatted here.
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s %s' % (__version__, __date__))

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-v', '--verbose', action='store_true',
                       help='report debug level status messages')
    group.add_argument('-q', '--quiet', action='store_true',
                       help='report errors only')

    commands = parser.add_subparsers(title='subcommands', dest='command',
                                     metavar='')

    sources = commands.add_parser(
        'sources', help='list available sources',
        description='List every data source immvoke knows how to fetch '
                    'from, along with a one-line description and its homepage.',
        formatter_class=CommonHelpFormatter)
    sources.add_subparsers(dest='action', metavar='').add_parser(
        'list', help='list the sources immvoke knows about',
        description='List every data source immvoke knows how to fetch '
                    'from, along with a one-line description and its homepage.',
        formatter_class=CommonHelpFormatter)

    _addSchemaParser(commands)
    _addReferenceParser(commands)
    for name, source in sorted(REGISTRY.items()):
        _addSourceParser(commands, name, source)
    _addSpecificityGroupParser(commands)
    addAnnotateParser(commands)

    return parser


def _addReferenceParser(commands):
    """Add the reference subcommand: validate and build from a reference folder."""
    reference = commands.add_parser(
        'reference', help='validate a germline reference folder and build '
                          'IgBLAST databases from it',
        description='Check that a folder of germline FASTAs is in a format '
                    'airrflow can use, and build the IgBLAST databases from it. '
                    'Files are recognised by name, in any directory layout: the '
                    'species and chain, with an optional source prefix and an '
                    'optional aa marker for translated V, as in human_IGHV.fasta '
                    'or imgt_human_IGHV.fasta.',
        formatter_class=CommonHelpFormatter)
    actions = reference.add_subparsers(dest='action', metavar='ACTION',
                                       required=True)

    build = actions.add_parser(
        'build', help='validate a reference folder and build IgBLAST databases',
        description='Validate the reference folder and build the IgBLAST '
                    'databases from it. With --check, only validate and report, '
                    'building nothing (and needing no makeblastdb).',
        formatter_class=CommonHelpFormatter)
    build.add_argument('folder', type=Path,
                       help='a reference_base tree or a flat folder of germline '
                            'FASTAs named <species>_<CHAIN>.fasta')
    build.add_argument('--out', type=Path, default=None,
                       help='directory to write igblast_base into; required '
                            'unless --check')
    build.add_argument('--check', action='store_true',
                       help='validate the folder and report what would build, '
                            'without building anything')
    build.add_argument('--species', nargs='+', choices=list(Reference.SPECIES),
                       default=None,
                       help='limit to these species; default is every species '
                            'found in the folder')
    build.add_argument('--map', dest='map_file', type=Path, default=None,
                       metavar='MANIFEST', help='a manifest declaring the species and chain of files whose names do not say, one per line: <file> <species> <CHAIN> [aa]. It overrides the naming rule, so it can also correct a file the rule misreads')

    diff = actions.add_parser(
        'diff', help='compare two germline reference folders allele by allele',
        description='Compare two reference folders allele by allele and report '
                    'what is identical, added, removed or changed. Files are '
                    'matched by name in any layout, and sequences are compared '
                    'without gaps, so a gapped and an ungapped copy of the same '
                    'allele are not a false difference. Exits non-zero when the '
                    'two references differ.',
        formatter_class=CommonHelpFormatter)
    diff.add_argument('reference_a', type=Path,
                      help='the baseline reference folder')
    diff.add_argument('reference_b', type=Path,
                      help='the reference folder to compare against it')
    diff.add_argument('--species', nargs='+', choices=list(Reference.SPECIES),
                      default=None,
                      help='limit to these species; default is every species '
                           'found in either folder')
    diff.add_argument('--map', dest='map_file', type=Path, default=None,
                      metavar='MANIFEST', help='a manifest declaring the species and chain of files whose names do not say, one per line: <file> <species> <CHAIN> [aa]. It overrides the naming rule, so it can also correct a file the rule misreads')

    show = actions.add_parser(
        'show', help='report what a reference folder is and where it came from',
        description='Read a reference folder\'s provenance sidecars -- '
                    'IMGT.yaml, AIRRC.yaml and immvoke_build.yaml -- and '
                    'report the release and sets it was built from, what was '
                    'built, and what the folder holds. Accepts a reference_base, '
                    'a directory containing one, or an igblast_base, which '
                    'carries copies of the same sidecars.',
        formatter_class=CommonHelpFormatter)
    show.add_argument('folder', type=Path,
                      help='the reference folder to describe')
    show.add_argument('--map', dest='map_file', type=Path, default=None,
                      metavar='MANIFEST', help='a manifest declaring the species and chain of files whose names do not say, one per line: <file> <species> <CHAIN> [aa]. It overrides the naming rule, so it can also correct a file the rule misreads')


def _addSchemaParser(commands):
    """Add the schema subcommand tree."""
    schema = commands.add_parser(
        'schema', help='inspect and refresh snapshots',
        description='Inspect the schema snapshot checked into the package, '
                    'or re-harvest it from the remote source when the '
                    'upstream API changes (new organisms, fields, values, ...).',
        formatter_class=CommonHelpFormatter)
    actions = schema.add_subparsers(dest='action', metavar='')

    show = actions.add_parser(
        'show', help='print a stored snapshot',
        description='Print the schema snapshot stored for a source: its '
                    'collections and how many fields each has, one '
                    'collection\'s fields and how many values each accepts, '
                    'or every value a single field accepts.',
        formatter_class=CommonHelpFormatter)
    show.add_argument('--source', required=True, choices=SCHEMA_SOURCES,
                      help='which source (or alias) to read the snapshot of, '
                           'a repertoire source or a specificity database')
    show.add_argument('--collection', default=None,
                      help='list this collection\'s fields; without it, print '
                           'one summary line per collection')
    show.add_argument('--field', default=None,
                      help='print every value this field accepts, one per line; '
                           'needs --collection')

    refresh = actions.add_parser(
        'refresh', help='re-harvest a snapshot',
        description='Contact a remote source, re-harvest its schema and '
                    'catalogs, and write the result over the packaged '
                    'snapshot (or alongside it, with --out). Previously '
                    'fetched detail-page enrichment is carried forward '
                    'unless --refresh-details asks to redo it.',
        formatter_class=CommonHelpFormatter)
    refresh.add_argument('--source', required=True, choices=SCHEMA_SOURCES,
                         help='which source to contact and re-harvest')
    refresh.add_argument('--out', default=None, type=Path,
                         help='directory to write into; without it the packaged '
                              'snapshot is rewritten in place')
    refresh.add_argument('--collection', action='append', default=None,
                         help='limit to one collection; repeatable')
    refresh.add_argument('--refresh-details', default='auto',
                         choices=['auto', 'all', 'none'],
                         help='fetch per unit detail pages for fields the search '
                              'results omit; auto fetches only units not already '
                              'read, all re-reads every unit')
    refresh.add_argument('--detail-limit', type=int, default=None,
                         help='stop after this many detail pages')


def _addSourceParser(commands, name, source):
    """Add one source's subcommand tree, with a level per collection."""
    schema = loadSchemaQuietly(name)

    parser = commands.add_parser(
        name, aliases=list(source.aliases), help=source.description,
        description='%s\n\nHomepage: %s' % (source.description, source.homepage),
        formatter_class=CommonHelpFormatter)
    actions = parser.add_subparsers(dest='action', metavar='ACTION',
                                    required=True)

    action_help = {
        'search': ('list matching data units',
                   'Search a collection for data units matching the given '
                   'filters and print (or save with --out) a summary of what '
                   'matched. Nothing is downloaded.'),
        'download': ('download matching data units',
                    'Search a collection for data units matching the given '
                    'filters, download each one, and optionally convert it '
                    'to AIRR and/or FASTA, writing a samplesheet for each '
                    'converted format.'),
    }
    for action, (helptext, description) in action_help.items():
        action_parser = actions.add_parser(action, help=helptext,
                                           description=description,
                                           formatter_class=CommonHelpFormatter)
        # The collection is a subcommand rather than a flag because the two
        # collections have genuinely different field sets. argparse cannot vary
        # options by a flag's value, so a flag would force a union and make
        # --help describe fields that do not apply.
        #
        # NB: a named metavar, not the '' used elsewhere, because argparse names
        # the missing argument by its metavar and an empty one produces a
        # required-argument error that says nothing.
        collections = action_parser.add_subparsers(dest='collection',
                                                   metavar='COLLECTION',
                                                   required=True)
        names = list(source.collections)
        if action == 'download' and source.output == 'reference' and len(names) > 1:
            names.append(ALL_SPECIES)

        for collection in names:
            # Passing help is what makes argparse list the collection at all.
            collection_help = source.collection_help.get(
                collection,
                'every species immvoke supports for this source (%s), into '
                'one reference_base. Not every species the source publishes'
                % ', '.join(source.collections))
            leaf = collections.add_parser(
                collection, help=collection_help,
                description='%s the %s collection (%s), optionally narrowed '
                            'down with the filter flags below.'
                            % (action.capitalize(), collection, collection_help),
                formatter_class=CommonHelpFormatter)
            addFilterArgs(leaf, schema, name,
                           source.collections[0] if collection == ALL_SPECIES
                           else collection)
            leaf.add_argument('--limit', type=int, default=None,
                              help='stop after this many units')
            if action == 'search':
                leaf.add_argument('-o', '--out', type=Path, default=None,
                                  help='write the hits to a TSV file')
            else:
                leaf.add_argument('--outdir', type=Path, required=True,
                                  help='directory to write into')
                leaf.add_argument('--dry-run', action='store_true',
                                  help='report what would be fetched, then stop')
                leaf.add_argument('--no-resume', action='store_true',
                                  help='re-download in full rather than '
                                       'continuing a partly fetched file')
                if source.output == 'reference':
                    # Reference sources build a germline reference_base rather
                    # than converting to AIRR, so they take the igblast options
                    # instead of the format and AIRR-strictness ones.
                    leaf.add_argument('--igblast', action='store_true',
                                      help='also build the IgBLAST databases '
                                           'from the reference; needs makeblastdb '
                                           'on PATH')
                    leaf.add_argument('--igblast-out', type=Path, default=None,
                                      help='where to write igblast_base; '
                                           'defaults to <outdir>/igblast_base')
                    leaf.add_argument('--from', dest='from_ref', type=Path,
                                      default=None, metavar='REFERENCE',
                                      help='re-download the versions pinned in a '
                                           'reference_base (its IMGT.yaml / '
                                           'AIRRC.yaml), or one of those files, '
                                           'instead of the latest')
                    leaf.add_argument('--compare', type=Path, default=None,
                                      metavar='REFERENCE',
                                      help='after building, compare the result '
                                           'allele by allele against this '
                                           'reference folder and report; a '
                                           'difference is a non-zero exit')
                    leaf.add_argument('--resolve-doi', action='store_true',
                                      help='resolve each OGRDB set\'s Zenodo DOI '
                                           'into AIRRC.yaml (scrapes the OGRDB '
                                           'web UI; ignored by imgt)')
                else:
                    leaf.add_argument('--format', action='append', dest='formats',
                                      choices=FORMATS,
                                      help='what to write, repeatable to write '
                                           'several; raw mirrors the source files '
                                           'untouched and is always written because '
                                           'the others are converted from it, so '
                                           'omitting this writes raw alone')
                    leaf.add_argument('--strict-airr', action='store_true',
                                      help='drop columns the AIRR schema does '
                                           'not define')


#: Request delay for the IgBLAST support mirror. The default is cautious because
#: IMGT and OGRDB are small academic servers; NCBI's file server and GitHub's raw
#: host are bulk services, and the mirror alone makes 100+ requests, so the
#: default would spend most of a build asleep rather than transferring.
MIRROR_DELAY = 0.05


def _sharedFilterTable(schema, source):
    """
    A table whose filters every table of the database shares, so 'all' can
    offer them too (IEDB's bcr and tcr filter on the same fields), or None
    when the tables filter differently or not at all.

    Arguments:
      schema (SourceSchema): the database's snapshot, or None.
      source (type): the specificity source class.

    Returns:
      str: the first table, or None.
    """
    if schema is None or not all(t in schema.collections for t in source.collections):
        return None
    fields = {schema.getCollection(t).fields for t in source.collections}
    if len(fields) != 1 or not next(iter(fields)):
        return None

    return source.collections[0]


def _addSpecificityParser(commands, name, source):
    """
    Add one specificity database's subcommand tree, with a level per table.

    A trimmed copy of `_addSourceParser`: no --format/samplesheet flags, since
    those are repertoire specific, and one extra synthetic 'all' table that
    loops every table the database offers.
    """
    schema = loadSchemaQuietly(name)

    parser = commands.add_parser(
        name, help=source.description,
        description='%s\n\nHomepage: %s' % (source.description, source.homepage),
        formatter_class=CommonHelpFormatter)
    actions = parser.add_subparsers(dest='action', metavar='ACTION',
                                    required=True)

    action_help = {
        'search': ('count what a table holds',
                   'Count what a table holds, optionally narrowed by the '
                   'filters: its rows and the most common values of its key '
                   'columns. The table is fetched into a temporary folder and '
                   'nothing else is written, unless --out asks for the '
                   'counts.'),
        'download': ('download a table',
                    'Fetch a table, optionally narrowed by the filters, and '
                    'write it as one normalized TSV: a real AIRR '
                    'rearrangement TSV for a table listed in the database\'s '
                    'airr_collections, a plain TSV of the source\'s own '
                    'columns otherwise.'),
    }
    for action, (helptext, description) in action_help.items():
        action_parser = actions.add_parser(action, help=helptext,
                                           description=description,
                                           formatter_class=CommonHelpFormatter)
        tables = action_parser.add_subparsers(dest='table', metavar='TABLE',
                                              required=True)

        entries = [(table, source.collection_help.get(table), table)
                   for table in source.collections]
        entries.append(('all', 'every table', None))

        for table_name, table_help, filter_table in entries:
            if filter_table is None:
                summary = 'every table this database offers'
            else:
                summary = 'the %s table: %s' % (table_name, table_help)
            leaf = tables.add_parser(
                table_name, help=table_help,
                description='%s %s.' % (action.capitalize(), summary),
                formatter_class=SpecificityHelpFormatter)
            if filter_table is not None:
                addFilterGroup(leaf, schema, source, filter_table)
            elif _sharedFilterTable(schema, source) is not None:
                addFilterGroup(leaf, schema, source, _sharedFilterTable(schema, source))
            if action == 'search':
                leaf.add_argument('--top', type=int, default=DEFAULT_TOP,
                                  help='show at most this many values per column')
                leaf.add_argument('-o', '--out', type=Path, default=None,
                                  help='also write every count to this TSV file')
            else:
                if filter_table is None:
                    leaf.add_argument('--limit', type=int, default=None,
                                      help='stop after this many tables')
                leaf.add_argument('--outdir', type=Path, required=True,
                                  help='directory to write into, one '
                                       'folder per database: '
                                       '<db>/<db>_<table>.tsv is the '
                                       'normalized table')
                leaf.add_argument('--dry_run', action='store_true',
                                  help='report what would be fetched, then stop')
                leaf.add_argument('--no_resume', action='store_true',
                                  help='re-download in full rather than '
                                       'continuing a partly fetched file')
                leaf.add_argument('--strict_airr', action='store_true',
                                  help='for a table written as AIRR, drop '
                                       'columns the AIRR schema does not '
                                       'define; no effect on other tables')


def _addSpecificityGroupParser(commands):
    """
    Add the `specificity` subcommand: `list`, one per registered database, and
    a synthetic `all` that fans out to every registered database.
    """
    group = commands.add_parser(
        'specificity', help='search and download specificity databases',
        description='Search and download specificity databases (reference '
                    'tables of epitopes, assays and receptor sequences), as '
                    'opposed to repertoire sequencing sources.',
        formatter_class=CommonHelpFormatter)
    dbs = group.add_subparsers(dest='db', metavar='DB', required=True)

    dbs.add_parser(
        'list', help='list the specificity databases immvoke knows about',
        description='List every specificity database immvoke knows how to '
                    'fetch from, along with a one-line description and its '
                    'homepage.',
        formatter_class=CommonHelpFormatter)

    for name, source in sorted(SPECIFICITY_REGISTRY.items()):
        _addSpecificityParser(dbs, name, source)

    all_group = dbs.add_parser(
        'all', help='every registered specificity database',
        description='Act on every registered specificity database at once.',
        formatter_class=CommonHelpFormatter)
    all_actions = all_group.add_subparsers(dest='action', metavar='ACTION',
                                           required=True)
    all_download = all_actions.add_parser(
        'download', help='download every table of every database',
        description='Download every table of every registered specificity '
                    'database. Filters cannot be applied here, since '
                    'different databases have unrelated fields; target one '
                    'database directly to filter it.',
        formatter_class=CommonHelpFormatter)
    all_download.add_argument('--outdir', type=Path, required=True,
                              help='directory to write into; one '
                                   'subdirectory per database')
    all_download.add_argument('--limit', type=int, default=None,
                              help='stop after this many tables, per database')
    all_download.add_argument('--dry_run', action='store_true',
                              help='report what would be fetched, then stop')
    all_download.add_argument('--strict_airr', action='store_true',
                              help='for a table written as AIRR, drop columns '
                                   'the AIRR schema does not define; no effect '
                                   'on other tables')


def makeClient(args):
    """Build the shared HTTP client."""
    return HttpClient()


def handleSources(args):
    """List registered sources."""
    for name, source in sorted(REGISTRY.items()):
        print('%-10s %s' % (name, source.description))
        if source.aliases:
            print('%-10s alias: %s' % ('', ', '.join(source.aliases)))
        print('%-10s %s' % ('', source.homepage))
        if source.license:
            print('%-10s license: %s' % ('', source.license))
        for paper in source.citation:
            print('%-10s cite: %s' % ('', paper))

    return 0


def handleSchemaShow(args):
    """Print a stored snapshot."""
    args.source = canonicalName(args.source)
    schema = loadSchema(args.source)

    if args.collection is None:
        print('source:    %s' % schema.source)
        print('harvested: %s by %s' % (schema.harvested, schema.harvested_by))
        for name in schema.collection_names:
            collection = schema.getCollection(name)
            print('  %s: %d fields (%s)'
                  % (name, len(collection.fields),
                     ', '.join(collection.field_names)))
        return 0

    collection = schema.getCollection(args.collection)
    if args.field is None:
        width = max(len(x.name) for x in collection.fields) if collection.fields else 0
        for item in collection.fields:
            if item.pseudo_values:
                kind = 'presence only'
            elif item.open_values:
                # The source checks the value itself, so what the snapshot
                # holds are examples; calling them the vocabulary would read
                # as though nothing else were accepted.
                kind = 'any value (%d examples)' % len(item.values)
            else:
                kind = '%d values' % len(item.values)
            print('%s %s' % (item.name.ljust(width), kind))
        return 0

    item = collection.getField(args.field)
    if item is None:
        raise ImmvokeError("no field '%s' in %s %s"
                             % (args.field, args.source, args.collection))
    if item.open_values:
        # On stderr, so that the values themselves still pipe cleanly.
        log.info("%s accepts any value %s checks; these are examples",
                 args.field, args.source)
    for value in item.values:
        print(value)

    return 0


def handleSchemaRefresh(args):
    """Re-harvest a snapshot and its catalogs."""
    args.source = canonicalName(args.source)
    client = makeClient(args)
    source = getAnySource(args.source, client)

    log.info('harvesting %s search schema', args.source)
    schema = source.harvestSchema()

    out = args.out
    if out is None:
        from importlib import resources
        out = Path(str(resources.files('immvoke').joinpath(
            'data/schemas', args.source)))

    written, changed = saveSchema(schema, out)
    log.info('%s %s', 'wrote' if changed else 'unchanged, left alone:', written)

    # Only a source whose units come from a catalog of its own has one to
    # refresh. For a germline set or a specificity database the snapshot is
    # the whole of it, so the harvest ends here rather than reaching for a
    # harvestCatalog that source never defined.
    if not hasattr(source, 'harvestCatalog'):
        return 0

    for collection in args.collection or list(source.collections):
        log.info('harvesting %s %s catalog', args.source, collection)
        rows = source.harvestCatalog(collection, schema=schema)

        path = out / ('%s_catalog.tsv' % collection)
        rows = Catalog.mergeEnrichment(Catalog.loadCatalog(path), rows)

        if args.refresh_details != 'none':
            force = args.refresh_details == 'all'
            pending = rows if force else [x for x in rows if Catalog.needsDetail(x)]
            if pending:
                log.info('enriching %d %s units from detail pages',
                         len(pending), collection)
                source.enrichCatalog(rows, limit=args.detail_limit, force=force)

        Catalog.saveCatalog(rows, path)
        log.info('wrote %s (%d units)', path, len(rows))

    return 0


def handleSearch(args):
    """List data units matching a query."""
    client = makeClient(args)
    source = getSource(args.source, client)

    query = source.validateQuery(args.collection, collectFilters(args))
    if args.limit is not None:
        query = type(query)(collection=query.collection, filters=query.filters,
                            limit=args.limit)

    units = source.searchUnits(query)
    total = sum(x.n_sequences or 0 for x in units)
    log.info('%d data units, %s sequences', len(units), format(total, ','))

    rows = [{'unit_id': x.unit_id, 'collection': x.collection,
             'n_unique_sequences': x.n_sequences or '', 'url': x.url,
             **{k: v for k, v in x.metadata.items()
                if k in Catalog.CATALOG_COLUMNS}}
            for x in units]

    if args.out is not None:
        Catalog.saveCatalog(rows, args.out)
        log.info('wrote %s', args.out)
    else:
        for unit in units:
            print('%-64s %10s' % (unit.unit_id, unit.n_sequences or ''))

    return 0


def loadMap(args):
    """
    Read the --map manifest, if one was given.

    Arguments:
      args (Namespace): parsed arguments.

    Returns:
      dict: the manifest, or None.
    """
    if getattr(args, 'map_file', None) is None:
        return None

    return Reference.loadReferenceMap(args.map_file)


def handleReferenceDiff(args):
    """Compare two reference folders and report; non-zero exit if they differ."""
    for folder in (args.reference_a, args.reference_b):
        if not folder.is_dir():
            raise ImmvokeError('no such reference folder: %s' % folder)

    diff = Reference.diffReference(args.reference_a, args.reference_b,
                                   species=args.species, mapping=loadMap(args))
    print(diff.summary())

    return 0 if diff.same else 1


def handleReferenceShow(args):
    """Report what a reference folder is and where it came from."""
    print(Reference.describeReference(args.folder, mapping=loadMap(args)))

    return 0


def handleReference(args):
    """Validate a reference folder and, unless --check, build its IgBLAST base."""
    if args.action == 'diff':
        return handleReferenceDiff(args)
    if args.action == 'show':
        return handleReferenceShow(args)

    if not args.folder.is_dir():
        raise ImmvokeError('no such reference folder: %s' % args.folder)

    plan = Reference.planReference(args.folder, species=args.species,
                                   mapping=loadMap(args))
    print(plan.summary())

    if not plan.ok:
        raise ImmvokeError('no databases can be built from %s; check the file '
                             'names against <species>_<CHAIN>.fasta' % args.folder)

    if args.check:
        return 0

    if args.out is None:
        raise ImmvokeError('--out is required to build; pass --check to only '
                             'validate the folder')

    report = Reference.buildFromPlan(plan, args.out, HttpClient(delay=MIRROR_DELAY))
    for path in Reference.writeBuildMetadata(
            args.out, report, args.folder, Provenance.timestamp()[:10],
            'immvoke %s' % __version__):
        log.info('wrote %s', path)
    log.info('wrote %s', args.out)

    return 0


def applyPins(source, from_ref, species):
    """
    Pin a reference source to the versions recorded in a reference_base.

    Reads the IMGT.yaml / AIRRC.yaml at --from and applies whichever pins the
    source can act on, so an imgt source takes the IMGT release, an ogrdb source
    the set versions, and the blend both. A --from whose pins none of the source
    can use is a mistake worth stopping on rather than silently fetching latest.

    Only the species being downloaded is read. A reference_base can hold several,
    downloaded weeks apart from different releases, and pinning a mouse download
    to human's release would quietly build something other than what was asked
    for.

    Arguments:
      source (ReferenceSource): the source about to download.
      from_ref (Path): a reference_base or an IMGT.yaml/AIRRC.yaml file.
      species (str): the species being downloaded.
    """
    from immvoke.Sources.Germline import ReferenceSource

    pins = Reference.loadReferencePins(from_ref)
    imgt, airrc = pins.get('imgt') or {}, pins.get('airrc') or {}
    can_release = type(source).pinRelease is not ReferenceSource.pinRelease

    applied = []
    entry = (imgt.get('species') or {}).get(species) or {}
    if entry.get('release') and can_release:
        source.pinRelease(entry['release'])
        applied.append('IMGT release %s' % entry['release'])

    sets = [item for item in airrc.get('sets') or []
            if item.get('species') == species]
    if sets and hasattr(source, 'pinSets'):
        source.pinSets(sets)
        applied.append('%d OGRDB set version(s)' % len(sets))

    if not applied:
        raise ImmvokeError('the reference at %s records no %s versions that %s '
                             'can re-download' % (from_ref, species, source.name))
    log.info('re-downloading pinned: %s', '; '.join(applied))


def handleReferenceDownload(args, source):
    """Download one species, or all, and build an airrflow reference_base."""
    species = (list(source.collections) if args.collection == ALL_SPECIES
               else [args.collection])
    if args.resolve_doi and hasattr(source, 'enableDoi'):
        source.enableDoi()

    outdir = Path(args.outdir)
    reference_dir = outdir / 'reference_base'
    provenance = []
    # Each species is pinned, fetched and built on its own; the provenance
    # sidecars merge across them, so several land in one reference_base.
    for name in species:
        if args.from_ref is not None:
            applyPins(source, args.from_ref, name)
        query = source.validateQuery(name, collectFilters(args))
        if args.limit is not None:
            query = type(query)(collection=query.collection,
                                filters=query.filters, limit=args.limit)
        units = source.searchUnits(query)
        log.info('%d germline files for %s %s', len(units), args.source, name)
        if args.dry_run:
            for unit in units:
                print('%-48s %s' % (unit.unit_id, unit.url))
            continue
        entries = []
        for unit in units:
            result = source.fetchUnit(unit, outdir / 'raw',
                                      resume=not args.no_resume)
            entries.append((unit, result.path))
            provenance.append(Provenance.buildUnitRecord(unit, result, outdir, {}))
        source.buildReference(entries, reference_dir).logSummary()
        for path in source.writeReferenceMetadata(reference_dir, units):
            log.info('wrote %s', path)

    if args.dry_run:
        log.info('dry run: nothing downloaded')
        return 0

    log.info('wrote %s', reference_dir)

    formats = ['reference']
    if args.igblast:
        igblast_out = args.igblast_out or (outdir / 'igblast_base')
        report = Reference.buildIgblastBase(reference_dir, igblast_out,
                                            HttpClient(delay=MIRROR_DELAY),
                                            species=species)
        Reference.writeBuildMetadata(igblast_out, report, reference_dir,
                                     Provenance.timestamp()[:10],
                                     'immvoke %s' % __version__)
        log.info('wrote %s', igblast_out)
        formats.append('igblast')

    record = Provenance.writeDownloadMetadata(
        outdir, args.source, args.collection, collectFilters(args), args.limit,
        formats, provenance, schema=source.schema, license=source.license,
        citation=source.citation)
    log.info('wrote %s', record)

    if args.compare is not None:
        if not args.compare.is_dir():
            raise ImmvokeError('no such reference folder: %s' % args.compare)
        diff = Reference.diffReference(args.compare, reference_dir,
                                       species=species)
        print(diff.summary())
        if not diff.same:
            log.warning('the built reference differs from %s', args.compare)
            return 1

    return 0


def handleDownload(args):
    """Download and optionally convert matching data units."""
    client = makeClient(args)
    source = getSource(args.source, client)

    # Reference sources (germline sets) build a reference_base instead of
    # converting repertoires to AIRR and writing a samplesheet.
    if source.output == 'reference':
        return handleReferenceDownload(args, source)

    query = source.validateQuery(args.collection, collectFilters(args))
    if args.limit is not None:
        query = type(query)(collection=query.collection, filters=query.filters,
                            limit=args.limit)

    units = source.searchUnits(query)
    formats = args.formats or ['raw']
    total = sum(x.n_sequences or 0 for x in units)

    log.info('%d data units, %s sequences, formats: %s',
             len(units), format(total, ','), ', '.join(formats))

    if args.dry_run:
        for unit in units:
            print('%-64s %10s' % (unit.unit_id, unit.n_sequences or ''))
        log.info('dry run: nothing downloaded')
        return 0

    outdir = Path(args.outdir)
    raw_dir = outdir / 'raw'
    # The raw mirror is always written, whether or not it was requested, because
    # converting reads from it. So 'raw' always has a bucket here even when the
    # user asked only for airr or fasta.
    written = {x: [] for x in set(formats) | {'raw'}}
    loci = {}
    provenance = []

    for unit in units:
        result = source.fetchUnit(unit, raw_dir, resume=not args.no_resume)
        written['raw'].append((unit, result.path))
        outputs = {}

        stem = unit.unit_id.replace('/', '_').replace('.csv.gz', '')
        if 'airr' in formats:
            _, chunks, report = source.convertUnit(result.path, unit)
            dest = outdir / 'airr' / ('%s.tsv' % stem)
            validation = Convert.writeAirr(chunks, dest, strict=args.strict_airr)
            Convert.writeValidationReport(validation, dest)
            log.info('%s: %d rows, %d invalid, %d rows in',
                     dest.name, validation['rows_checked'],
                     validation['rows_invalid'], report['rows_in'])
            loci[unit.unit_id] = report['loci']
            written['airr'].append((unit, dest))
            outputs['airr'] = dest

        if 'fasta' in formats:
            _, chunks, report = source.convertUnit(result.path, unit)
            dest = outdir / 'fasta' / ('%s.fasta' % stem)
            Convert.writeFasta(chunks, dest)
            loci.setdefault(unit.unit_id, report['loci'])
            written['fasta'].append((unit, dest))
            outputs['fasta'] = dest

        provenance.append(
            Provenance.buildUnitRecord(unit, result, outdir, outputs))

    # A samplesheet is a derived artifact of a data format, so one is written per
    # converted format rather than one ambiguous sheet naming a single file.
    for fmt in ('airr', 'fasta'):
        if written.get(fmt):
            sheet = outdir / ('samplesheet_airrflow_%s.tsv' % fmt)
            buildSamplesheet(written[fmt], sheet, args.collection, outdir,
                             loci=loci)
            log.info('wrote %s', sheet)

    # Written for every run, including raw-only ones: the raw mirror is the part
    # of the output that cannot be regenerated from anything else here.
    record = Provenance.writeDownloadMetadata(
        outdir, args.source, args.collection, collectFilters(args), args.limit,
        formats, provenance, schema=source.schema, license=source.license,
        citation=source.citation)
    log.info('wrote %s', record)

    return 0


def _resolveSpecificityTables(source, args):
    """
    Return the tables a specificity action should act on.

    Arguments:
      source (SourceBase): the specificity source.
      args (Namespace): parsed arguments, carrying `table` and, for
        `download all`, `limit`.

    Returns:
      list: table names. `--limit` caps how many are used when `table` is
      'all'; a single named table, and search, take no --limit.
    """
    if args.table != 'all':
        return [args.table]

    tables = list(source.collections)
    limit = getattr(args, 'limit', None)

    return tables[:limit] if limit is not None else tables


def _writeSpecificityTsv(chunks, out):
    """
    Write normalized specificity records as a TSV, streaming chunk by chunk.

    Unlike `Convert.writeAirr`, this performs no AIRR schema validation:
    specificity columns are a source's own, not AIRR rearrangement fields.

    Arguments:
      chunks (iterable): DataFrames of normalized records.
      out (Path): output path.

    Returns:
      int: rows written.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    header_written = False
    with open(out, 'w', newline='') as handle:
        for frame in chunks:
            frame.to_csv(handle, sep='\t', index=False, header=not header_written)
            header_written = True
            written += len(frame)

    if not header_written:
        out.write_text('')

    return written


def handleSpecificityList(args):
    """List registered specificity databases."""
    for name, source in sorted(SPECIFICITY_REGISTRY.items()):
        print('%-10s %s' % (name, source.description))
        print('%-10s %s' % ('', source.homepage))
        if source.license:
            print('%-10s license: %s' % ('', source.license))
        for paper in source.citation:
            print('%-10s cite: %s' % ('', paper))

    return 0


def handleSpecificitySearch(args):
    """Count what a table (or every table) holds; see Specificity.Search."""
    source = getSpecificitySource(args.db, makeClient(args))

    return searchTables(source, _resolveSpecificityTables(source, args),
                        collectFilters(args), top=args.top, out=args.out)


def handleSpecificityDownload(args):
    """
    Download one or more tables: a raw mirror plus one normalized TSV each,
    all in the database's own folder, `<outdir>/<db>/` -- so several
    databases can share one outdir without mixing their files or their
    download_metadata.yml.
    """
    client = makeClient(args)
    source = getSpecificitySource(args.db, client)
    filters = collectFilters(args)
    tables = _resolveSpecificityTables(source, args)

    outdir = Path(args.outdir) / args.db
    raw_dir = outdir / 'raw'

    records = []
    formats_written = set()
    for table in tables:
        query = source.validateQuery(table, filters)
        units = source.searchUnits(query)

        if args.dry_run:
            for unit in units:
                print('%-16s %s' % (table, unit.unit_id))
            continue

        for unit in units:
            start = time.monotonic()
            Report.started(source, table, filters)
            result = source.fetchUnit(unit, raw_dir, resume=not args.no_resume)
            _, chunks, report = source.convertUnit(result.path, unit)
            tally = Report.Tally(source, table)
            chunks = tally.wrap(chunks)

            # Named after the database as well as the table, so the file
            # says where it came from once copied out of <db>/.
            dest = outdir / ('%s_%s.tsv' % (args.db, table))
            if table in source.airr_collections:
                validation = Convert.writeAirr(chunks, dest,
                                               strict=args.strict_airr)
                Convert.writeValidationReport(validation, dest)
                rows = validation['rows_checked']
                output_format = 'airr'
            else:
                validation = None
                rows = _writeSpecificityTsv(chunks, dest)
                output_format = 'tsv'
            formats_written.add(output_format)
            Report.completed(dest, time.monotonic() - start, report['rows_in'], rows,
                             tally, validation)

            records.append(Provenance.buildUnitRecord(
                unit, result, outdir, {output_format: dest}))

    if args.dry_run:
        log.info('dry run: nothing downloaded')
        return 0

    if not source.keep_raw:
        shutil.rmtree(raw_dir, ignore_errors=True)
        log.debug('removed %s', raw_dir)

    # A specificity download represents the database's current state, not one
    # more addition to a growing set the way OAS downloads are, so this always
    # replaces the metadata file rather than accumulating old runs into it --
    # see Specificity.Provenance for why that is a separate writer.
    record = writeSpecificityMetadata(
        outdir, args.db, args.table, filters, getattr(args, 'limit', None), records,
        schema=source.schema, license=source.license, citation=source.citation,
        formats=formats_written or ('tsv',))
    Report.finished(record)

    return 0


def handleSpecificityDownloadAll(args):
    """
    Download every table of every registered specificity database, each
    into its own `<outdir>/<db>/` (see handleSpecificityDownload).
    """
    outdir = Path(args.outdir)

    for name in sorted(SPECIFICITY_REGISTRY):
        log.info('downloading specificity source %s', name)
        handleSpecificityDownload(Namespace(
            db=name, table='all', outdir=outdir, limit=args.limit,
            dry_run=args.dry_run, no_resume=False,
            strict_airr=args.strict_airr))

    return 0


def main():
    """
    Parse the commandline and dispatch to the selected subcommand.

    Returns:
      int: process exit status.
    """
    parser = getArgParser()
    args = parser.parse_args()

    setupLogging(verbose=args.verbose, quiet=args.quiet)

    if args.command is None:
        parser.print_help(sys.stderr)
        return 1

    try:
        if args.command == 'sources':
            return handleSources(args)

        if args.command == 'schema':
            if args.action == 'show':
                return handleSchemaShow(args)
            if args.action == 'refresh':
                return handleSchemaRefresh(args)
            parser.parse_args([args.command, '--help'])

        if args.command == 'reference':
            return handleReference(args)

        source_name = canonicalName(args.command) if args.command else None
        if source_name in REGISTRY:
            # The action and collection levels are required subparsers, so
            # argparse has already rejected a commandline missing either. The
            # command may be an alias (e.g. 'airrc'); resolve it to the canonical
            # source so schema and provenance use one name.
            args.source = source_name
            if args.action == 'search':
                return handleSearch(args)
            if args.action == 'download':
                return handleDownload(args)

        if args.command == 'specificity':
            if args.db == 'list':
                return handleSpecificityList(args)
            if args.db == 'all':
                return handleSpecificityDownloadAll(args)
            # The action and table levels are required subparsers, so argparse
            # has already rejected a commandline missing either.
            if args.action == 'search':
                return handleSpecificitySearch(args)
            if args.action == 'download':
                return handleSpecificityDownload(args)

        if args.command == 'annotate':
            return handleAnnotate(args)

        parser.print_help(sys.stderr)
        return 1
    except ImmvokeError as error:
        log.error('%s', error)
        return 1


if __name__ == '__main__':
    sys.exit(main())
