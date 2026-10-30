"""
`immvoke specificity <db> search`, and the filter flags it shares with download

search counts what a table holds once the filters are applied: its rows,
the distinct values of key columns (receptor groups, epitopes, ...) and the
most common values of others (locus, species, outcome, ...), as each source
lists them in search_columns. It fetches the table the way download does,
but into a temporary folder that is removed afterwards (IEDB's bulk ZIP
lands in its cache, so a download straight after reuses it), and writes
nothing else unless --out asks for the counts as a TSV.

The filter flags of both actions come from the source's schema snapshot,
grouped under 'filters', each saying what it filters on (the source's
filter_help) and laying its values out in columns.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import collections
import csv
import logging
import tempfile
import textwrap
from pathlib import Path

# Immvoke imports
from immvoke.Commandline import CommonHelpFormatter

log = logging.getLogger(__name__)

#: Most common values shown per column unless --top says otherwise.
DEFAULT_TOP = 10

#: Marks a help line holding a value list, for SpecificityHelpFormatter to
#: lay out once it knows the width; VALUE_SEP separates values.
VALUE_LIST = '\x00'
VALUE_SEP = '\x1f'

#: A short vocabulary is written inline rather than as an indented list.
INLINE_VALUES = 4


class SpecificityHelpFormatter(CommonHelpFormatter):
    """
    CommonHelpFormatter for the specificity table commands: it keeps a help
    string's own line breaks, lays a filter's values out so that no value is
    split across lines, and leaves '(default: None)' and '(default: False)'
    off a flag that is simply off unless given.
    """

    def _split_lines(self, text, width):
        lines = []
        for line in text.splitlines():
            if line.startswith(VALUE_LIST):
                lines.extend(_flow(line[len(VALUE_LIST):].split(VALUE_SEP), width))
            else:
                lines.extend(textwrap.wrap(line, width) or [''])

        return lines

    def _get_help_string(self, action):
        if action.default is None or action.default is False:
            return action.help

        return super()._get_help_string(action)


def _flow(values, width):
    """
    Pack values onto lines, comma separated, breaking only between values:
    many species names have spaces, so ordinary wrapping would split them.

    Arguments:
      values (list): the values, in order.
      width (int): characters available.

    Returns:
      list: the lines, each indented by two spaces.
    """
    lines, line = [], ''
    for i, value in enumerate(values):
        piece = value + (',' if i < len(values) - 1 else '')
        if line and len(line) + 1 + len(piece) > width - 2:
            lines.append('  ' + line)
            line = piece
        else:
            line = '%s %s' % (line, piece) if line else piece
    if line:
        lines.append('  ' + line)

    return lines


def addFilterGroup(parser, schema, source, table):
    """
    Add a table's filter flags to its subcommand, as a 'filters' group.

    Each flag's placeholder and description come from the source's
    filter_help, and its values from the schema snapshot: listed inline when
    there are a few, in columns when there are more, as examples when the
    source checks the value itself (Field.open_values). The source's
    filter_switches, if any, follow as flags that take no value.

    Arguments:
      parser (ArgumentParser): the table's subparser.
      schema (SourceSchema): the source's snapshot, or None if none is
        installed.
      source (type): the source class, for its filter_help.
      table (str): the table whose fields to add.
    """
    if schema is None or table not in schema.collections:
        return
    fields = schema.getCollection(table).fields
    if not fields:
        return

    group = parser.add_argument_group(
        'filters', 'Each narrows the table down; leave one out to keep every value.')
    for item in fields:
        metavar, what = source.filter_help.get(item.name, ('VALUE', ''))
        lines = [what] if what else []
        if item.pseudo_values:
            lines.append('defined or undefined')
        elif item.open_values:
            lines.append('e.g. %s' % ', '.join('"%s"' % x for x in item.values))
        elif len(item.values) <= INLINE_VALUES:
            lines.append('one of: %s' % ', '.join(item.values))
        else:
            lines.append('%d values:' % len(item.values))
            lines.append(VALUE_LIST + VALUE_SEP.join(item.values))

        group.add_argument(item.flag, dest='filter_%s' % item.name,
                           metavar=metavar, default=None, help='\n'.join(lines))
    for name, what in getattr(source, 'filter_switches', {}).items():
        group.add_argument('--%s' % name, dest='filter_%s' % name, action='store_true',
                           default=None, help=what)


def searchTable(source, table, filters, workdir):
    """
    Count what one table holds once the filters are applied.

    Arguments:
      source (SourceBase): the specificity source.
      table (str): the table.
      filters (dict): the filter flags given, field to value.
      workdir (Path): where the table is fetched to.

    Returns:
      dict: `table`, `filters` (those given), `rows_in` (before filtering),
      `rows` (after), `distinct` (label -> count) and `top` (label ->
      Counter of labelled values), from the source's search_columns.
    """
    query = source.validateQuery(table, filters)
    columns = source.search_columns.get(table, ())
    seen = {column: set() for column, _, how in columns if how == 'distinct'}
    counts = {column: collections.Counter() for column, _, how in columns
              if how == 'top'}

    rows_in = rows = 0
    for unit in source.searchUnits(query):
        result = source.fetchUnit(unit, workdir)
        _, chunks, report = source.convertUnit(result.path, unit)
        for chunk in chunks:
            rows += len(chunk)
            for column in seen.keys() & set(chunk.columns):
                seen[column].update(chunk[column].tolist())
            for column in counts.keys() & set(chunk.columns):
                counts[column].update(chunk[column].tolist())
        rows_in += report['rows_in']

    top = {}
    for column, label, how in columns:
        if how == 'top':
            labelled = collections.Counter()
            for value, count in counts[column].items():
                text = '' if value is None else str(value)
                labelled[source.labelValue(column, text)] += count
            top[label] = labelled

    return {'table': table,
            'filters': {k: v for k, v in query.filters.items() if v != '*'},
            'rows_in': rows_in, 'rows': rows,
            'distinct': {label: len(seen[column] - {''})
                         for column, label, how in columns if how == 'distinct'},
            'top': top}


def formatSummary(source, summary, top=DEFAULT_TOP):
    """
    Render a searchTable summary for the terminal.

    Arguments:
      source (SourceBase): the source, for its name and table description.
      summary (dict): from searchTable.
      top (int): most common values shown per column.

    Returns:
      str: the summary, ending in a newline.
    """
    table = summary['table']
    lines = ['%s %s: %s' % (source.name, table, source.collection_help.get(table, ''))]
    if summary['filters']:
        lines.append('filters: %s' % ', '.join('%s = %s' % item
                                                for item in summary['filters'].items()))
    lines.append('')

    figures = [('rows', '%s of %s' % (format(summary['rows'], ','),
                                      format(summary['rows_in'], ',')))]
    figures += [(label, format(count, ','))
                for label, count in summary['distinct'].items()]
    pad = max(len(label) for label, _ in figures) + 2
    lines += ['  %s%s' % (label.ljust(pad), value) for label, value in figures]

    for label, counter in summary['top'].items():
        if not counter:
            continue
        shown = counter.most_common(top)
        more = len(counter) - len(shown)
        lines += ['', '  %s%s' % (label, ' (top %d of %d)' % (len(shown), len(counter))
                                  if more else '')]
        digits = len(format(shown[0][1], ','))
        lines += ['    %s  %s' % (format(count, ',').rjust(digits), value)
                  for value, count in shown]

    return '\n'.join(lines) + '\n'


def writeSummaries(summaries, out):
    """
    Write every count of the summaries as a TSV: table, field, value, rows.

    Arguments:
      summaries (list): from searchTable.
      out (Path): the TSV to write.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['table', 'field', 'value', 'rows'])
        for summary in summaries:
            table = summary['table']
            writer.writerow([table, 'rows', '', summary['rows']])
            for label, count in summary['distinct'].items():
                writer.writerow([table, label, '', count])
            for label, counter in summary['top'].items():
                for value, count in counter.most_common():
                    writer.writerow([table, label, value, count])


def searchTables(source, tables, filters, top=DEFAULT_TOP, out=None):
    """
    Search one or more tables and print what each holds.

    Arguments:
      source (SourceBase): the specificity source.
      tables (list): the tables to search.
      filters (dict): the filter flags given, field to value.
      top (int): most common values shown per column.
      out (Path): also write every count here as a TSV, or None.

    Returns:
      int: process exit status.
    """
    summaries = []
    with tempfile.TemporaryDirectory(prefix='immvoke-search-') as workdir:
        for table in tables:
            summary = searchTable(source, table, filters, Path(workdir))
            summaries.append(summary)
            print(formatSummary(source, summary, top=top))

    if out is not None:
        writeSummaries(summaries, out)
        log.info('wrote %s', out)

    return 0
