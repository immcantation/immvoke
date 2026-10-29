"""
What `immvoke specificity <db> download` tells the user as it runs

A download says when it starts (the table and the filters given), lets the
source say what it is doing in between (IEDB: whether its export comes from
the cache, the epitopes it looks up), and ends with what it wrote: how long
it took, how many rows were kept, the distinct counts the source lists in
search_columns (receptor groups, epitopes, ...), the AIRR validation result
and the files. It is all status, so it goes through logging like the rest of
immvoke's status output and --quiet silences it.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def humanDuration(seconds):
    """
    A duration for people: '8 s', '2 min 5 s', '3 h 10 min'.

    Arguments:
      seconds (float): the duration.

    Returns:
      str: the duration, to the second under an hour, to the minute above.
    """
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


def finished(record):
    """
    Say where the download's provenance was written.

    Arguments:
      record (Path): the download_metadata.yml written.
    """
    log.info('Provenance (source, license, citation, command): %s', record)
