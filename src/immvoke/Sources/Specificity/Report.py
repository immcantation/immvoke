"""
What `immvoke specificity <db> download` and `immvoke annotate` tell the
user as they run

Both say when they start and what they are about to do, report what they
are doing in between, and end with how long it took and what came of it:
a download with the rows kept and the counts the source lists in
search_columns, an annotation with how many query rows matched something.
It is all status, so it goes through logging like the rest of immvoke's
status output and --quiet silences it.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def humanDuration(seconds):
    """
    A duration for people: '0.4 s', '8 s', '2 min 5 s', '3 h 10 min'.

    Arguments:
      seconds (float): the duration.

    Returns:
      str: the duration, to a tenth under a second (where rounding to the
      second would say '0 s' for work that did happen), to the second under
      an hour, to the minute above.
    """
    if seconds < 1:
        return '%.1f s' % seconds
    seconds = int(round(seconds))
    if seconds < 60:
        return '%d s' % seconds
    if seconds < 3600:
        return '%d min %d s' % divmod(seconds, 60)

    return '%d h %d min' % (seconds // 3600, seconds % 3600 // 60)


def humanPath(path):
    """
    A path for people: under the home folder as '~/...', else as given.

    Arguments:
      path (Path): the path.

    Returns:
      str: the path to show.
    """
    path = Path(path)
    try:
        return '~/%s' % path.resolve().relative_to(Path.home())
    except ValueError:
        return str(path)


class Tally:
    """
    Count a table's distinct values as its chunks stream past on their way
    to the file, for the columns the source lists as 'distinct' in
    search_columns, so the summary costs no second pass over the table.

    Arguments:
      source (SourceBase): the specificity source.
      table (str): the table being written.
    """

    def __init__(self, source, table):
        self.columns = [(column, label) for column, label, how
                        in source.search_columns.get(table, ()) if how == 'distinct']
        self.seen = {column: set() for column, _ in self.columns}

    def wrap(self, chunks):
        """
        Pass the chunks through, counting as they go.

        Arguments:
          chunks (iterable): DataFrames of normalized records.

        Yields:
          pandas.DataFrame: each chunk, unchanged.
        """
        for chunk in chunks:
            for column, _ in self.columns:
                if column in chunk.columns:
                    self.seen[column].update(chunk[column].tolist())
            yield chunk

    def distinct(self):
        """
        Returns:
          list: (label, count) per column, blanks not counted.
        """
        return [(label, len(self.seen[column] - {'', None}))
                for column, label in self.columns]


def started(source, table, filters):
    """
    Say a table's download has started, and with which filters.

    Arguments:
      source (SourceBase): the specificity source.
      table (str): the table.
      filters (dict): the filter flags given, field to value.
    """
    given = {k: v for k, v in filters.items() if v not in (None, '*')}
    log.info('Download started: %s %s (%s)', source.name, table,
             source.collection_help.get(table, ''))
    log.info('  filters: %s', ', '.join('%s = %s' % x for x in given.items())
             if given else 'none, every record kept')


def completed(dest, seconds, rows_in, rows, tally, validation=None):
    """
    Say a table's download has finished, with what it holds.

    Arguments:
      dest (Path): the table written.
      seconds (float): how long the table took.
      rows_in (int): rows in the source's table, before filtering.
      rows (int): rows written.
      tally (Tally): the distinct counts, from the rows written.
      validation (dict): from Convert.writeAirr, or None for a table that is
        not AIRR.
    """
    figures = [('rows', '%s of %s kept' % (format(rows, ','), format(rows_in, ',')))]
    figures += [(label, format(count, ',')) for label, count in tally.distinct()]
    if validation is not None:
        invalid = validation['rows_invalid']
        verdict = ('every row valid'
                   if not invalid and validation['header_valid'] is not False
                   else '%s of %s rows invalid' % (format(invalid, ','),
                                                   format(validation['rows_checked'], ',')))
        figures.append(('AIRR validation', '%s (%s)' % (
            verdict, dest.with_name(dest.stem + '.airr_validation.log').name)))
    pad = max(len(label) for label, _ in figures) + 2

    log.info('Download completed in %s: %s', humanDuration(seconds), dest)
    for label, value in figures:
        log.info('  %s%s', label.ljust(pad), value)


#: sequence_region -> what it compares, for annotationStarted.
_REGIONS = {
    'CDR3H': 'the heavy chain CDR3 (cdr3_aa)',
    'VDJ': 'the whole V(D)J amino acid sequence (sequence_aa)',
    'CDRs': 'CDR1, CDR2 and CDR3 pooled',
}

#: partition_by -> what a pair must share before it is compared at all.
_PARTITIONS = {
    None: 'any pair of rows is compared',
    'V': 'only rows with the same V gene are compared',
    'VJ': 'only rows with the same V and J gene are compared',
    'VJL': 'only rows with the same V gene, J gene and CDR3 length are compared',
}


def annotationStarted(query, reference, options):
    """
    Say an annotation has started, and what it is about to compare.

    Arguments:
      query (dict): the query's metadata record: `path` and `rows`.
      reference (dict): the reference's, which may also name its `source`
        database.
      options (dict): the run's options, as annotateSpecificity records
        them.
    """
    def named(record, fallback):
        # As the caller wrote it: resolving a relative path to an absolute
        # one only makes a familiar argument harder to recognize.
        return record['path'] or fallback

    log.info('Annotation started: %s rows of %s', format(query['rows'], ','),
             named(query, 'a table passed from Python'))
    log.info('  against %s rows of %s', format(reference['rows'], ','),
             named(reference, 'a reference passed from Python'))

    if options['paired']:
        log.info('  comparing: heavy and light chain of the same antibody, '
                 'same V gene, J gene and exact CDR3')
        return

    threshold = options['identity_threshold']
    log.info('  comparing: %s', _REGIONS.get(options['sequence_region'],
                                             options['sequence_region']))
    log.info('  a hit is: %s', 'an exact match' if threshold is None else
             '%g%% identity or better, %s' % (threshold * 100,
                                              options['distance_metric']))
    log.info('  %s', _PARTITIONS.get(options['partition_by'],
                                     options['partition_by']))


def annotationCompleted(seconds, results, paired):
    """
    Say an annotation has finished, and what it found.

    Arguments:
      seconds (float): how long it took.
      results (dict): the run's counts, as annotateSpecificity records them.
      paired (bool): whether this was a paired run.
    """
    log.info('Annotation completed in %s', humanDuration(seconds))
    if paired:
        figures = [('query antibodies', format(results['query_antibodies'], ','))]
        figures += [(label, format(count, ','))
                    for label, count in results['classification'].items()]
    else:
        figures = [
            ('query rows matched',
             '%s of %s' % (format(results['query_rows_with_hits'], ','),
                           format(results['query_rows_compared'], ','))),
            ('hits within the threshold', format(results['hits'], ',')),
        ]
    if 'full_matrix_rows' in results:
        figures.append(('pairs in the full matrix',
                        format(results['full_matrix_rows'], ',')))
    pad = max(len(label) for label, _ in figures) + 2
    for label, value in figures:
        log.info('  %s%s', label.ljust(pad), value)


def finished(record):
    """
    Say where the download's provenance was written.

    Arguments:
      record (Path): the download_metadata.yml written.
    """
    log.info('Provenance (source, license, citation, command): %s', record)
