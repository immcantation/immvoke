"""
Antigen specificity annotation: query-vs-reference receptor comparison

For each row of a query AIRR table, find every row of a fixed reference AIRR
table it agrees with, and report its best hits -- one annotation per query
row that matched something (its `sequence_id`, `n_hits`, `identity`,
`hit_ids` and each hit's `antigen`, `epitope` and
`antigen_source_species`), plus, on request, every
hit within the threshold (full_dist_matrix). This is the last pipeline stage: download a
specificity database's table, which lands as an AIRR rearrangement TSV in
the database's own folder (`immvoke specificity iedb download bcr`), then
annotate a user's own AIRR table against it as the reference, into an output
folder of its own.

annotateSpecificity() is the entry point `immvoke annotate` wraps: it picks
the comparison from a few knobs (sequence_region, partition_by,
paired, identity_threshold, distance_metric) -- see its section at the end of
this module -- and runs it through Symscan.runSymscanAnnotation, symscan's
symmetric-deletion index. Rows beyond symscan's limits are scored directly
instead (Distances.py, _scoreSkippedRows), so no row misses a hit.

runPairedExactAnnotation builds on an exact match for the case where a
query row is one chain of an antibody rather than a self-contained unit: it
runs it once per chain (heavy, light) and reconciles the two per antibody,
classifying each as agreeing with the same reference antibody on both
chains, one chain only, or neither.

V/J-gene agreement goes through Genes.py (`gene_cols`), which
parses gene calls rather than comparing them as strings; `exact_cols`
compares plain strings. Substitution-matrix-weighted (BLOSUM) annotation is
not included.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
import sys
from pathlib import Path

import numpy
import pandas
import yaml

# Immvoke imports
from immvoke.Exceptions import AnnotationError
from immvoke.Provenance import (
    DOWNLOAD_METADATA,
    METADATA_VERSION,
    commandLine,
    timestamp,
)
from immvoke.Version import __version__

log = logging.getLogger(__name__)

#: Reference column read as each hit's identifier when the caller does not
#: say otherwise -- the AIRR rearrangement primary key.
DEFAULT_ID_COL = 'sequence_id'
#: An unresolved residue: under exact and Hamming matching, a position where
#: either side has it never counts as a difference (see Symscan.py).
WILDCARD = 'X'

#: runPairedExactAnnotation's fixed column names -- one antibody's heavy and
#: light chain are two rows of the same AIRR table, linked by `cell_id` and
#: told apart by `locus`; there is nothing to make configurable here beyond
#: what the AIRR schema already fixes these fields to mean.
CELL_ID_COL = 'cell_id'
LOCUS_COL = 'locus'
HEAVY_LOCUS = 'IGH'
LIGHT_LOCI = ('IGK', 'IGL')
#: The AIRR count column that picks one heavy and one light chain per query
#: cell in paired annotation, unless the caller names another.
DEFAULT_COUNT_FIELD = 'duplicate_count'

#: What annotateSpecificity writes into its outdir: the result, the record
#: of the run and, with full_dist_matrix, every hit within the threshold --
#: each name prefixed with the reference's database and table
#: (_outputPrefix), e.g. iedb_bcr_annotation.tsv.
ANNOTATION_FILE = 'annotation.tsv'
ANNOTATION_METADATA = 'annotation_metadata.yml'
FULL_MATRIX_FILE = 'annotation_full_matrix.txt'

#: What the result reports about each best hit -> the reference column it
#: comes from (IEDB's, see Iedb.py); a reference without one leaves it blank.
HIT_ANNOTATIONS = (
    ('antigen', 'iedb_antigen'),
    ('epitope', 'iedb_epitope_name'),
    ('antigen_source_species', 'iedb_antigen_source_species'),
)
#: Separates hits in hit_ids and HIT_ANNOTATIONS' columns -- not a comma,
#: which IEDB's discontinuous epitope names ('G58, L59, ...') are full of.
HIT_SEP = ';'
#: Separates one hit's values when it stands for several reference rows:
#: IEDB lists a receptor once per epitope, all under its sequence_id.
VALUE_SEP = ' | '

#: Slack added before truncating a fraction-of-length threshold to a count,
#: so float error can't drop it by one (e.g. (1 - 0.9) * 10 is 0.99999...).
_FRACTION_EPSILON = 1e-9


def _allowedDifferences(fraction, length):
    """
    Resolve a fraction-of-length tolerance to a whole number of differences.

    Arguments:
      fraction (float): the tolerated fraction, e.g. 0.1.
      length (int): the query's length the fraction applies to.

    Returns:
      int: floor(fraction * length), robust to float representation error.
    """
    return int(fraction * length + _FRACTION_EPSILON)


def _mismatchPositions(a, b):
    """
    Every substitution mismatch between `a` and `b`, a WILDCARD on either
    side of a position never counting as one.

    Arguments:
      a (tuple): one row's values, one per comparison column.
      b (tuple): the row it is being compared to, same shape.

    Returns:
      list: (column index, 0-based position) pairs.
    """
    return [(col, pos)
           for col, (x, y) in enumerate(zip(a, b))
           for pos, (cx, cy) in enumerate(zip(x, y))
           if cx != cy and cx != WILDCARD and cy != WILDCARD]


def _formatPositions(positions, n_cols):
    """
    Render mismatch positions as a comma-separated string.

    Arguments:
      positions (list): (column index, position) pairs, from
        _mismatchPositions.
      n_cols (int): how many comparison columns are in play.

    Returns:
      str: bare positions ("3,7") for the common single-column case,
      column-prefixed ("0:3,1:2") when there is more than one -- a bare
      position is ambiguous once more than one column is being compared.
    """
    if n_cols == 1:
        return ','.join(str(pos) for _, pos in positions)

    return ','.join('%d:%d' % (col, pos) for col, pos in positions)


def _rowValues(frame, columns):
    """
    Read one or more columns as a list of per-row tuples.

    Arguments:
      frame (pandas.DataFrame): the table to read from.
      columns (list): column names, in the order values should be compared.

    Returns:
      list: one tuple per row, values in `columns` order, blank standing in
      for a missing value.
    """
    # astype(object): an all-blank column read as float takes a '' fill too.
    return list(frame[columns].astype(object).fillna('').itertuples(index=False,
                                                                    name=None))


def _exactValues(frame, exact_cols):
    """
    _rowValues for the exact-agreement gate, or one empty tuple per row when
    there is no gate.

    Arguments:
      frame (pandas.DataFrame): the table to read from.
      exact_cols (list): column names, or None.

    Returns:
      list: one tuple per row.
    """
    return _rowValues(frame, exact_cols) if exact_cols else [()] * len(frame)


# ── paired heavy+light exact annotation ─────────────────────────────────────

def _hitsByCellId(query_chain, chain_result, cell_id_col):
    """
    Fold one chain's exact match result into a hit set per antibody.

    Arguments:
      query_chain (pandas.DataFrame): the query rows matched for this
        chain, in the same order (its index reset to 0..n-1, matching
        `chain_result`'s `query_row`).
      chain_result (pandas.DataFrame): the exact match result for
        `query_chain`.
      cell_id_col (str): column in `query_chain` identifying the antibody a
        chain row belongs to.

    Returns:
      dict: query antibody id to the set of reference antibody ids (read off
      the reference's `cell_id_col`) whose corresponding chain agreed with
      it.
    """
    cell_ids = query_chain[cell_id_col].reset_index(drop=True)
    hits = {}
    for row in chain_result.itertuples(index=False):
        if not row.hit_ids:
            continue
        cell_id = cell_ids.iloc[row.query_row]
        hits.setdefault(cell_id, set()).update(row.hit_ids.split(';'))

    return hits


def runPairedExactAnnotation(queries, reference, cdr3_col='cdr3_aa',
                             vgene_col='v_call', jgene_col='j_call',
                             threads=1, return_pairs=False):
    """
    Classify each query antibody by whether the *same* reference antibody
    agrees with it on the heavy chain, the light chain, both, or neither.

    An antibody is two rows of an AIRR table -- a heavy chain and, where
    present, a light chain -- sharing one `cell_id` and told apart by
    `locus` (`IGH` for heavy, `IGK`/`IGL` for light). Agreement on a chain
    is the same V-gene and J-gene (any allele, see Genes.py) and
    exact CDR3, run once for the heavy rows and once for the light rows; a
    reference antibody counts toward "both" only if its `cell_id` appears in
    *both* chains' hit sets for the query antibody, not merely in each
    independently -- a query whose heavy chain happens to match one
    reference antibody and whose light chain happens to match an unrelated
    one is `conflicting`, not `both`.

    Arguments:
      queries (pandas.DataFrame): query AIRR rows, both chains, with `locus`
        and `cell_id` columns.
      reference (pandas.DataFrame): reference AIRR rows, same shape.
      cdr3_col (str): the CDR3 amino acid column.
      vgene_col (str): the V-gene call column.
      jgene_col (str): the J-gene call column.
      threads (int): symscan's thread pool; see Symscan.runSymscanAnnotation.
      return_pairs (bool): also return the (query, reference) antibody pairs
        behind the hit ids.

    Returns:
      pandas.DataFrame: one row per query `cell_id`: `classification` (one
      of `both`, `heavy_only`, `light_only`, `conflicting`, `neither`),
      `heavy_hit_ids`, `light_hit_ids` (every reference antibody agreeing on
      that chain alone, semicolon joined, sorted) and `both_hit_ids` (their
      intersection -- non-empty only when `classification` is `both`); its
      `attrs['skipped_query_rows']` counts the query chain rows symscan
      could not compare (see Symscan.py). With `return_pairs`, a tuple: that table, and one row per (query antibody,
      reference antibody) agreeing on either chain -- `cell_id`, `hit_id`
      and `matched_chains` (`both`, `heavy` or `light`).
    """
    query_heavy = queries[queries[LOCUS_COL] == HEAVY_LOCUS]
    query_light = queries[queries[LOCUS_COL].isin(LIGHT_LOCI)]
    reference_heavy = reference[reference[LOCUS_COL] == HEAVY_LOCUS]
    reference_light = reference[reference[LOCUS_COL].isin(LIGHT_LOCI)]

    query_heavy = query_heavy.reset_index(drop=True)
    query_light = query_light.reset_index(drop=True)
    reference_heavy = reference_heavy.reset_index(drop=True)
    reference_light = reference_light.reset_index(drop=True)

    gene_cols = [c for c in (vgene_col, jgene_col) if c is not None] or None
    # Deferred: Symscan builds on this module's helpers.
    from immvoke.Sources.Specificity.Symscan import runSymscanAnnotation

    def _exact(query_chain, reference_chain):
        return runSymscanAnnotation(query_chain, reference_chain, [cdr3_col], 0,
                                    id_col=CELL_ID_COL, gene_cols=gene_cols,
                                    threads=threads)

    heavy_result = _exact(query_heavy, reference_heavy)
    light_result = _exact(query_light, reference_light)

    heavy_hits = _hitsByCellId(query_heavy, heavy_result, CELL_ID_COL)
    light_hits = _hitsByCellId(query_light, light_result, CELL_ID_COL)

    rows, pairs = [], []
    for cell_id in sorted(set(queries[CELL_ID_COL])):
        heavy = heavy_hits.get(cell_id, set())
        light = light_hits.get(cell_id, set())
        both = heavy & light
        for hit in sorted(heavy | light):
            matched = 'both' if hit in both else 'heavy' if hit in heavy else 'light'
            pairs.append({'cell_id': cell_id, 'hit_id': hit, 'matched_chains': matched})

        if both:
            classification = 'both'
        elif heavy and light:
            classification = 'conflicting'
        elif heavy:
            classification = 'heavy_only'
        elif light:
            classification = 'light_only'
        else:
            classification = 'neither'

        rows.append({
            'cell_id': cell_id,
            'classification': classification,
            'heavy_hit_ids': ';'.join(sorted(heavy)),
            'light_hit_ids': ';'.join(sorted(light)),
            'both_hit_ids': ';'.join(sorted(both)),
        })

    result = pandas.DataFrame(rows)
    result.attrs['skipped_query_rows'] = (heavy_result.attrs['skipped_query_rows']
                                          + light_result.attrs['skipped_query_rows'])
    if not return_pairs:
        return result

    return result, pandas.DataFrame(pairs, columns=['cell_id', 'hit_id',
                                                    'matched_chains'])


# ── high level API: annotateSpecificity() ───────────────────────────────────
#
# A file-or-DataFrame-in, DataFrame-and-files-out entry point, for a caller
# who wants "match this receptor table against that reference table". It
# picks the comparison from a small set of knobs (paired, count_field,
# sequence_region, partition_by, identity_threshold, distance_metric) -- the
# flow `immvoke annotate` exposes on the commandline:
#
#     paired  -> one heavy and one light chain per query cell (count_field)
#             -> runPairedExactAnnotation (V, J, exact CDR3 on both chains)
#     else    -> heavy chain only -> sequence_region columns
#             -> partition_by gate (none, V, VJ, VJ + CDR3 length)
#             -> identity_threshold: 1.0 exact, below it Hamming/Levenshtein
#             -> symscan; rows beyond its limits are scored directly instead
#     then    -> each query row's best hits (_bestHits)
#             -> outdir: ANNOTATION_FILE, ANNOTATION_METADATA and, with
#                full_dist_matrix, FULL_MATRIX_FILE, each named after the
#                reference (e.g. iedb_bcr_annotation.tsv)

#: sequence_region -> AIRR column(s) to compare. 'CDRs' pools CDR1-3, so it
#: only runs under Hamming distance (or an exact match), never Levenshtein:
#: edit distance summed over independently aligned columns is not the edit
#: distance of any one alignment, the way Hamming distance pools cleanly.
SEQUENCE_REGION_COLUMNS = {
    'cdr3h': ('cdr3_aa',),
    'vdj': ('sequence_aa',),
    'cdrs': ('cdr1_aa', 'cdr2_aa', 'cdr3_aa'),
}

#: partition_by -> gene call columns two rows must agree on (same gene, any
#: allele; see Genes.py) before their sequence counts as a hit. 'VJL' also
#: gates on exact CDR3 length; see annotateSpecificity.
PARTITION_GENE_COLUMNS = {
    'v': ('v_call',),
    'vj': ('v_call', 'j_call'),
    'vjl': ('v_call', 'j_call'),
}

#: distance_metric values accepted ('levenstein' is a common misspelling).
DISTANCE_METRICS = {'hamming': 'hamming', 'levenshtein': 'levenshtein',
                    'levenstein': 'levenshtein'}

#: Synthetic gate column added to both tables under partition_by='VJL'; it
#: only feeds the exact-column gate and is never written to the result.
_CDR3_LENGTH_COL = '__cdr3_length__'


def _resolveSequenceRegion(sequence_region):
    """
    Arguments:
      sequence_region (str): 'CDR3H', 'VDJ' or 'CDRs', case insensitive.

    Returns:
      list: the AIRR column(s) to compare.

    Raises:
      AnnotationError: if `sequence_region` is not recognized.
    """
    columns = SEQUENCE_REGION_COLUMNS.get(str(sequence_region).strip().lower())
    if columns is None:
        raise AnnotationError("unknown sequence_region %r; choose one of: CDR3H, "
                              "VDJ, CDRs" % (sequence_region,))

    return list(columns)


def _resolvePartitionBy(partition_by):
    """
    Arguments:
      partition_by (str): None, 'V', 'VJ' or 'VJL', case insensitive.

    Returns:
      tuple: (gate columns (list), whether to also gate on CDR3 length).

    Raises:
      AnnotationError: if `partition_by` is not recognized.
    """
    if partition_by is None:
        return [], False

    key = str(partition_by).strip().lower()
    columns = PARTITION_GENE_COLUMNS.get(key)
    if columns is None:
        raise AnnotationError("unknown partition_by %r; choose one of: None, V, "
                              "VJ, VJL" % (partition_by,))

    return list(columns), key == 'vjl'


def _resolveDistanceMetric(distance_metric):
    """
    Arguments:
      distance_metric (str): 'hamming' or 'levenshtein', case insensitive.

    Returns:
      str: 'hamming' or 'levenshtein'.

    Raises:
      AnnotationError: if `distance_metric` is not recognized.
    """
    metric = DISTANCE_METRICS.get(str(distance_metric).strip().lower())
    if metric is None:
        raise AnnotationError("unknown distance_metric %r; choose 'hamming' or "
                              "'levenshtein'" % (distance_metric,))

    return metric


def _readAirrTsv(table):
    """
    Read an AIRR rearrangement TSV as plain text columns.

    Arguments:
      table (str or pandas.DataFrame): the file to read, or an already
        loaded table (returned as is).

    Returns:
      pandas.DataFrame: every column as a string, blank cells as '' rather
      than NaN -- the row comparisons here expect a value to always be a
      string.
    """
    if isinstance(table, pandas.DataFrame):
        return table

    return pandas.read_csv(table, sep='\t', dtype=str, na_filter=False)


def _requireColumns(frame, columns, label, why):
    """
    Raise a readable error, instead of a KeyError, for a missing column.

    Arguments:
      frame (pandas.DataFrame): the table to check.
      columns (list): the column names it needs.
      label (str): 'query' or 'reference', for the message.
      why (str): what needs them, for the message.

    Raises:
      AnnotationError: if any of `columns` is missing from `frame`.
    """
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise AnnotationError('%s has no %s column, needed by %s'
                              % (label, ', '.join(repr(c) for c in missing), why))


def _restrictToHeavyChain(frame, label):
    """
    Keep only the heavy chain (locus IGH), if `frame` has a locus column; a
    table without one is kept whole.

    Arguments:
      frame (pandas.DataFrame): the table to filter.
      label (str): 'query' or 'reference', for messages.

    Returns:
      pandas.DataFrame: the filtered table, index reset to 0..n-1.

    Raises:
      AnnotationError: if no row is left -- rather than silently annotating
        nothing.
    """
    if LOCUS_COL in frame.columns:
        filtered = frame[frame[LOCUS_COL] == HEAVY_LOCUS]
        log.info('%s: %d/%d rows are locus %s', label, len(filtered), len(frame),
                 HEAVY_LOCUS)
        if not len(filtered):
            raise AnnotationError(
                '%s has no heavy chain rows (locus %s) to annotate; its loci '
                'are: %s' % (label, HEAVY_LOCUS,
                             ', '.join(sorted(set(frame[LOCUS_COL]))) or 'none'))
        frame = filtered

    return frame.reset_index(drop=True)


def _withCellId(frame, label):
    """
    Drop the rows with a blank cell_id: nothing links such a chain to a
    partner, so paired annotation cannot pair it -- and must not lump every
    such row together as one antibody.

    Arguments:
      frame (pandas.DataFrame): the table, with a cell_id column.
      label (str): 'query' or 'reference', for the log line only.

    Returns:
      pandas.DataFrame: `frame` without those rows.
    """
    blank = frame[CELL_ID_COL].astype(object).fillna('').astype(str).str.strip() == ''
    if blank.any():
        log.warning('%s: %d rows have no cell_id to pair them by; left out of '
                    'paired annotation', label, int(blank.sum()))

    return frame[~blank]


def _oneChainPerCell(frame, count_field, label):
    """
    Keep one heavy chain (IGH) and one light chain (IGK or IGL) per cell_id:
    the one with the highest `count_field`. On a tie, a blank count, or no
    such column at all, the row that comes first in the file is kept. Rows
    of other loci are kept as they are.

    Arguments:
      frame (pandas.DataFrame): the table, with cell_id and locus columns.
      count_field (str): the AIRR count column to rank chains by, e.g.
        'duplicate_count' or 'umi_count'.
      label (str): 'query' or 'reference', for the log line only.

    Returns:
      pandas.DataFrame: `frame` without the extra chains, rows in their
      original order.
    """
    frame = frame.reset_index(drop=True)
    if count_field in frame.columns:
        counts = pandas.to_numeric(frame[count_field], errors='coerce')
        # A blank count ranks below any number.
        counts = counts.fillna(float('-inf'))
    else:
        counts = pandas.Series(0, index=frame.index)

    drop = []
    report = {}
    for name, chain in (('heavy', frame[LOCUS_COL] == HEAVY_LOCUS),
                        ('light', frame[LOCUS_COL].isin(LIGHT_LOCI))):
        rows = frame.index[chain]
        cells = frame.loc[rows, CELL_ID_COL]
        several = cells[cells.duplicated(keep=False)]
        if several.empty:
            continue
        # Highest count first; a stable sort keeps file order among ties.
        ranked = counts[several.index].sort_values(ascending=False, kind='stable')
        best = frame.loc[ranked.index, CELL_ID_COL].drop_duplicates()
        drop.extend(several.index.difference(best.index))
        report[name] = several.nunique()

    if not drop:
        return frame
    if count_field in frame.columns:
        rule = 'the one with the highest %s (the first row on a tie)' % count_field
    else:
        rule = 'the first row (no %s column to rank them by)' % count_field
    log.info('%s: %s; kept %s', label,
             ', '.join('%d cells had several %s chains' % (n, name)
                       for name, n in report.items()), rule)

    return frame.drop(index=drop)


def _requireAirr(frame, label):
    """
    Require `frame` to be an AIRR rearrangement table: every field the AIRR
    schema marks required present as a column (its values may be blank).

    Arguments:
      frame (pandas.DataFrame): the table to check.
      label (str): 'query' or 'reference', for the message.

    Raises:
      AnnotationError: naming the missing required fields.
    """
    from airr.schema import RearrangementSchema, ValidationError

    try:
        RearrangementSchema.validate_header(list(frame.columns))
    except ValidationError as error:
        raise AnnotationError('%s is not an AIRR rearrangement table: %s'
                              % (label, error)) from None


def _runSingleChain(queries, reference, compare_cols, gene_cols, exact_cols,
                    threshold, metric, threads):
    """
    Run one single-chain comparison through symscan, keeping every hit.

    Arguments:
      queries (pandas.DataFrame): the query rows, already filtered.
      reference (pandas.DataFrame): the reference rows, already filtered.
      compare_cols (list): the column(s) to compare.
      gene_cols (list): gene call columns that must name the same gene,
        possibly empty.
      exact_cols (list): columns that must agree exactly, possibly empty.
      threshold (float): identity in (0.0, 1.0]; 1.0 is an exact match.
      metric (str): 'hamming' or 'levenshtein'.
      threads (int): symscan's thread pool; see annotateSpecificity.

    Returns:
      tuple: (every hit within the threshold, one row per (query row, hit
      id), closest first within a query row -- see
      Symscan.runSymscanAnnotation's `return_pairs`; the query rows beyond
      symscan's limits, scored directly instead).
    """
    exact = threshold >= 1.0
    # An exact match is 0 substitutions under either metric -- the same
    # thing, and the only way to run exact agreement on pooled 'CDRs'.
    distance_type = 'hamming' if exact else metric
    tolerance = 0 if exact else 1 - threshold
    gates = dict(id_col=DEFAULT_ID_COL, exact_cols=exact_cols or None,
                 gene_cols=gene_cols or None)

    # Deferred: Symscan builds on this module's helpers.
    from immvoke.Sources.Specificity.Symscan import runSymscanAnnotation

    options = dict(distance_type=distance_type, threshold_is_fraction=not exact,
                   **gates)
    result, pairs = runSymscanAnnotation(queries, reference, compare_cols, tolerance,
                                         threads=threads, return_pairs=True,
                                         warn_skipped=False, **options)
    skipped = result.attrs['skipped_rows']
    found = _scoreSkippedRows(skipped, queries, reference, compare_cols, tolerance,
                              options)
    if len(found):
        pairs = (pandas.concat([pairs, found[pairs.columns]])
                 .sort_values(['query_row', 'mismatches'], kind='stable')
                 .drop_duplicates(['query_row', 'reference_id'])
                 .reset_index(drop=True))

    return pairs, skipped


def _scoreSkippedRows(rows, queries, reference, compare_cols, tolerance, options):
    """
    Find the hits of the query rows beyond symscan's limits (see Symscan.py)
    by scoring their pairs directly (Distances.allPairDistances), so no row
    misses a hit -- slower, but only for those rows.

    Arguments:
      rows (list): those query rows, runSymscanAnnotation's
        `attrs['skipped_rows']`.
      queries, reference, compare_cols, tolerance: as given to symscan.
      options (dict): the rest of symscan's keyword arguments.

    Returns:
      pandas.DataFrame: their hits as pairs, see
      Symscan.runSymscanAnnotation's `return_pairs`; empty when there are
      none.
    """
    rows = numpy.array(rows, dtype=numpy.int64)
    if not len(rows):
        return pandas.DataFrame()
    log.info("%d query rows are beyond symscan's limits (see Symscan.py); scoring "
             'their pairs directly', len(rows))

    # Deferred: Distances builds on this module's helpers.
    from immvoke.Sources.Specificity.Distances import allPairDistances

    subset = queries.iloc[rows].reset_index(drop=True)
    pairs = pandas.concat(list(allPairDistances(subset, reference, compare_cols,
                                                tolerance, **options)))
    found = pairs[pairs['within_threshold']].drop(columns='within_threshold')
    found = found.assign(query_row=rows[found['query_row'].to_numpy(dtype=int)])

    return found.reset_index(drop=True)


def _hitValues(reference):
    """
    What each reference id says for HIT_ANNOTATIONS: its distinct values
    across every row with that id, VALUE_SEP joined in reference order.

    Arguments:
      reference (pandas.DataFrame): the reference rows that were compared.

    Returns:
      dict: result column -> {reference id -> text}; a column the reference
      lacks maps every id to ''.
    """
    ids = reference[DEFAULT_ID_COL].astype(str)
    values = {}
    for name, column in HIT_ANNOTATIONS:
        if column not in reference.columns:
            values[name] = {}
            continue
        text = reference[column].astype(object).fillna('').astype(str)
        values[name] = {
            hit: VALUE_SEP.join(dict.fromkeys(x for x in group if x))
            for hit, group in text.groupby(ids.to_numpy(), sort=False)}

    return values


def _bestHits(pairs, n_queries, reference):
    """
    Summarize each query row by its best hits: the ones with the highest
    identity, all of them on a tie.

    Arguments:
      pairs (pandas.DataFrame): every hit, from _runSingleChain.
      n_queries (int): how many query rows were compared; one without a
        hit still gets a row.
      reference (pandas.DataFrame): the reference rows that were compared,
        for HIT_ANNOTATIONS.

    Returns:
      pandas.DataFrame: one row per query row, including those without a
      hit (annotateSpecificity drops those once it has counted them):
      `query_row`, `n_hits` (how many best hits), `identity` (theirs, to 4
      decimals), `hit_ids` and HIT_ANNOTATIONS' columns, one entry per hit,
      HIT_SEP joined in the same order; 0 hits and blanks for a query row
      without one.
    """
    names = [name for name, _ in HIT_ANNOTATIONS]
    columns = ['query_row', 'n_hits', 'identity', 'hit_ids'] + names
    result = pandas.DataFrame({'query_row': numpy.arange(n_queries)})

    if len(pairs):
        # One query row's identities share its length, so the fewest
        # mismatches is the highest identity, without comparing floats.
        fewest = pairs.groupby('query_row')['mismatches'].transform('min')
        best = pairs[pairs['mismatches'] == fewest]
        values = _hitValues(reference)
        grouped = best.groupby('query_row', sort=True)
        hits = pandas.DataFrame({
            'n_hits': grouped.size(),
            'identity': grouped['identity'].first().round(4),
            'hit_ids': grouped['reference_id'].agg(HIT_SEP.join)})
        for name in names:
            hits[name] = grouped['reference_id'].agg(
                lambda ids, lookup=values[name]: HIT_SEP.join(lookup.get(x, '')
                                                              for x in ids))
        result = result.join(hits, on='query_row')

    result = result.reindex(columns=columns)
    result['n_hits'] = result['n_hits'].fillna(0).astype(int)
    for name in ['hit_ids'] + names:
        result[name] = result[name].astype(object).fillna('').astype(str)

    return result


def _withQueryIds(table, queries):
    """
    Name each query by its own `sequence_id` instead of its position, so
    the output joins straight back to the query file. A query without a
    `sequence_id` column keeps `query_row`.

    Arguments:
      table (pandas.DataFrame): a result or pairs table, with `query_row`
        indexing `queries`.
      queries (pandas.DataFrame): the query rows that were compared, in
        `query_row` order.

    Returns:
      pandas.DataFrame: `table` with `sequence_id` in place of `query_row`.
    """
    if DEFAULT_ID_COL not in queries.columns:
        return table
    ids = queries[DEFAULT_ID_COL].astype(object).fillna('').astype(str).to_numpy()
    named = table.drop(columns='query_row')
    named.insert(0, DEFAULT_ID_COL, ids[table['query_row'].to_numpy(dtype=int)])

    return named


def _tableRecord(table, frame):
    """
    Describe an input table for the annotation metadata.

    Arguments:
      table (str or pandas.DataFrame): the table as the caller passed it.
      frame (pandas.DataFrame): the table as read.

    Returns:
      dict: its `path` (None for a DataFrame) and `rows`.
    """
    path = None if isinstance(table, pandas.DataFrame) else str(table)

    return {'path': path, 'rows': len(frame)}


def _referenceSource(reference):
    """
    Read which database a reference table came from, off the
    download_metadata.yml `immvoke specificity` writes in the database's
    folder, beside the table.

    Arguments:
      reference (str or pandas.DataFrame): the reference as the caller
        passed it.

    Returns:
      dict: `source` (e.g. 'iedb') and `downloaded` (when), or empty when
      the reference is a DataFrame, or no readable record sits beside it.
    """
    if isinstance(reference, pandas.DataFrame):
        return {}
    path = Path(reference).parent / DOWNLOAD_METADATA
    try:
        record = yaml.safe_load(path.read_text())
    except (OSError, yaml.YAMLError):
        return {}
    if not isinstance(record, dict) or 'source' not in record:
        return {}

    runs = record.get('runs') or [{}]
    return {'source': record['source'], 'downloaded': runs[-1].get('finished')}


def _outputPrefix(reference, source):
    """
    Name the reference for the output files, so runs against different
    databases stay apart even in one outdir.

    Arguments:
      reference (str or pandas.DataFrame): the reference as the caller
        passed it.
      source (str): its database, from _referenceSource, or None.

    Returns:
      str: '<db>_<table>' (e.g. 'iedb_bcr') for a table `immvoke
      specificity` downloaded, read off its file name (iedb_bcr.tsv), with
      the database added when the name lacks it (an older bcr.tsv); the file
      name alone for any other reference file; 'reference' for a DataFrame.
    """
    if isinstance(reference, pandas.DataFrame):
        return 'reference'
    table = Path(reference).name.split('.')[0]
    if not source or table.startswith(source + '_'):
        return table

    return '%s_%s' % (source, table)


def _writeOutputs(outdir, prefix, result, matrix, metadata):
    """
    Write the annotation bundle into `outdir`: ANNOTATION_FILE,
    ANNOTATION_METADATA and, when there is a matrix, FULL_MATRIX_FILE, each
    named `<prefix>_<name>`.

    Arguments:
      outdir (str): the folder to write into, or None to skip writing.
      prefix (str): the reference's name, from _outputPrefix.
      result (pandas.DataFrame): the annotation result.
      matrix (pandas.DataFrame): every hit within the threshold, one row per
        (query, hit), or None when full_dist_matrix was not asked for.
      metadata (dict): the run record, see annotateSpecificity.
    """
    if outdir is None:
        return

    # Rendered once up front, so a value YAML cannot represent fails before
    # any file is written rather than leaving a partial bundle behind.
    yaml.safe_dump(metadata)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    names = ['%s_%s' % (prefix, name) for name in (ANNOTATION_FILE, ANNOTATION_METADATA)]
    result.to_csv(out / names[0], sep='\t', index=False)
    if matrix is not None:
        names.append('%s_%s' % (prefix, FULL_MATRIX_FILE))
        # Identity to 4 decimals: 0.4231, not 0.42307692307692313.
        matrix.round({'identity': 4}).to_csv(out / names[-1], sep='\t', index=False)
    metadata['finished'] = timestamp()
    (out / names[1]).write_text(yaml.safe_dump(metadata, sort_keys=False))
    log.info('wrote %s into %s', ', '.join(names), out)


def annotateSpecificity(query, reference, sequence_region='CDR3H',
                        partition_by=None, paired=False, identity_threshold=None,
                        distance_metric='hamming', threads=1,
                        outdir=None, count_field=DEFAULT_COUNT_FIELD,
                        full_dist_matrix=False):
    """
    Annotate a query receptor table against a reference specificity table.

    Arguments:
      query (str or pandas.DataFrame): the query AIRR rearrangement TSV, or
        the table itself.
      reference (str or pandas.DataFrame): the reference AIRR rearrangement
        TSV -- e.g. a specificity database's downloaded table -- or the
        table itself. It must be an AIRR rearrangement table: every field
        the AIRR schema requires present as a column (values may be blank),
        as immvoke's own AIRR writer always produces.
      sequence_region (str): which sequence to compare: 'CDR3H' (the heavy
        chain CDR3, cdr3_aa),
        'VDJ' (sequence_aa, the full V(D)J amino acid sequence) or 'CDRs'
        (cdr1_aa, cdr2_aa and cdr3_aa pooled -- Hamming only). Ignored when
        `paired` is True.
      partition_by (str): only compare row pairs in the same partition: None
        (no restriction), 'V' (same V-gene), 'VJ' (same V- and J-gene) or
        'VJL' (same V-gene, J-gene and CDR3 length -- needs `cdr3_aa` in both
        tables whatever `sequence_region` is). Ignored when `paired` is True.
      paired (bool): classify each query antibody (`cell_id`) by V-gene,
        J-gene and exact CDR3 agreement on both its heavy and light chain
        with the same reference antibody (runPairedExactAnnotation); rows
        with a blank `cell_id` are left out, with a warning. When False,
        only heavy chain (locus IGH) rows are compared; a table without a
        locus column is compared whole.
      count_field (str): under `paired`, a query cell with several heavy
        (or light) chains keeps only the one with the highest value in this
        AIRR count column -- 'duplicate_count' by default, or e.g.
        'umi_count' -- or the first in the file on a tie, a blank count, or
        no such column.
      identity_threshold (float): minimum sequence identity in (0.0, 1.0];
        None or 1.0 requires an exact match, a lower value allows up to
        floor((1 - identity_threshold) * query length) differences. Ignored
        when `paired` is True.
      distance_metric (str): 'hamming' (substitutions only) or 'levenshtein'
        (substitutions, insertions and deletions; a single-column
        `sequence_region` only, i.e. not 'CDRs'). Only matters when
        `identity_threshold` is below 1.0.
      threads (int): symscan's thread pool size (-1 for every core); only
        the first symscan run in a process sets it.
      outdir (str): if given, also write into this folder, each file named
        after the reference's database and table (e.g.
        iedb_bcr_annotation.tsv; see _outputPrefix): ANNOTATION_FILE (the
        result) and ANNOTATION_METADATA (the options, the query and
        reference with the reference's source database, the counts, and
        the immvoke version, time and command).
      full_dist_matrix (bool): also write FULL_MATRIX_FILE into `outdir`:
        every hit within the threshold, not only the best ones, one row per
        (query, hit) -- `sequence_id`, `hit_id`, `identity`; under
        `paired`, `cell_id`, `hit_id` and `matched_chains` (`both`, `heavy`
        or `light`).

    Returns:
      pandas.DataFrame: when `paired`, one row per query `cell_id` (see
      runPairedExactAnnotation); otherwise one row per heavy chain query
      row that matched something, a row without a hit left out, with
      `sequence_id` (the query's own; `query_row`, its position among the
      heavy chain rows, for a query without that column), then its best
      hits, the ones
      with the highest identity (all of them on a tie): `n_hits`,
      `identity`, `hit_ids`, and each hit's `antigen`, `epitope` and
      `antigen_source_species` from the reference (HIT_ANNOTATIONS), one
      entry per hit, HIT_SEP joined in hit_ids order, a hit standing for
      several reference rows joining its values with VALUE_SEP.

    Raises:
      AnnotationError: if the reference is not an AIRR rearrangement table,
        an argument is not recognized, a needed column is missing, no heavy
        chain row is left to annotate, an argument combination cannot run
        (e.g. 'CDRs' under 'levenshtein'). Single-chain rows beyond
        symscan's limits are scored directly instead (_scoreSkippedRows);
        under `paired`, such rows (over 14 X's) are skipped with a warning.
    """
    queries = _readAirrTsv(query)
    ref = _readAirrTsv(reference)
    _requireAirr(ref, 'reference')
    source = _referenceSource(reference)
    prefix = _outputPrefix(reference, source.get('source'))

    metadata = {
        'immvoke_metadata_version': METADATA_VERSION,
        'generated_by': 'immvoke %s' % __version__,
        'finished': None,
        'query': _tableRecord(query, queries),
        'reference': dict(_tableRecord(reference, ref), **source),
        # .item(): a numpy scalar (e.g. a threshold from numpy.arange) as the
        # plain Python value YAML can write.
        'options': {name: value.item() if hasattr(value, 'item') else value
                    for name, value in (('sequence_region', sequence_region),
                                        ('partition_by', partition_by),
                                        ('paired', paired),
                                        ('identity_threshold', identity_threshold),
                                        ('distance_metric', distance_metric),
                                        ('count_field', count_field),
                                        ('full_dist_matrix', full_dist_matrix))},
    }
    # Run through the commandline rather than called from Python, where
    # sys.argv is some other program's.
    if Path(sys.argv[0]).stem == 'immvoke':
        metadata['command'] = commandLine()

    if paired:
        if partition_by is not None or identity_threshold is not None:
            log.warning('paired annotation always compares V-gene, J-gene and '
                        'exact CDR3; partition_by and identity_threshold are '
                        'ignored')
        paired_cols = [CELL_ID_COL, LOCUS_COL, 'cdr3_aa', 'v_call', 'j_call']
        _requireColumns(queries, paired_cols, 'query', 'paired annotation')
        _requireColumns(ref, paired_cols, 'reference', 'paired annotation')
        queries = _withCellId(queries, 'query')
        ref = _withCellId(ref, 'reference')
        chain_loci = (HEAVY_LOCUS,) + LIGHT_LOCI
        for frame, label in ((queries, 'query'), (ref, 'reference')):
            if not frame[LOCUS_COL].isin(chain_loci).any():
                raise AnnotationError('%s has no antibody chain rows (locus %s) '
                                      'to pair' % (label, ', '.join(chain_loci)))
        queries = _oneChainPerCell(queries, count_field, 'query')
        result, pairs = runPairedExactAnnotation(queries, ref, threads=threads,
                                                 return_pairs=True)
        classes = {label: int(count) for label, count
                   in result['classification'].value_counts().items()}
        log.info('%d query antibodies: %s', len(result),
                 ', '.join('%s=%d' % item for item in classes.items()))
        metadata['results'] = {
            'query_antibodies': len(result), 'classification': classes,
            'skipped_query_rows': result.attrs['skipped_query_rows'],
            'pairs': len(pairs)}
        _writeOutputs(outdir, prefix, result, pairs if full_dist_matrix else None,
                      metadata)

        return result

    threshold = 1.0 if identity_threshold is None else float(identity_threshold)
    if not 0 < threshold <= 1.0:
        raise AnnotationError('identity_threshold must be in (0.0, 1.0]; got %s'
                              % (identity_threshold,))
    compare_cols = _resolveSequenceRegion(sequence_region)
    gene_cols, gate_on_cdr3_length = _resolvePartitionBy(partition_by)
    metric = _resolveDistanceMetric(distance_metric)
    if threshold < 1.0 and metric == 'levenshtein' and len(compare_cols) != 1:
        raise AnnotationError("distance_metric='levenshtein' compares a single "
                              "column; sequence_region=%r pools %d, which only "
                              "'hamming' supports"
                              % (sequence_region, len(compare_cols)))

    needed = compare_cols + gene_cols + (['cdr3_aa'] if gate_on_cdr3_length else [])
    why = 'sequence_region=%r, partition_by=%r' % (sequence_region, partition_by)
    _requireColumns(queries, needed, 'query', why)
    _requireColumns(ref, needed + [DEFAULT_ID_COL], 'reference', why)

    queries = _restrictToHeavyChain(queries, 'query')
    ref = _restrictToHeavyChain(ref, 'reference')
    exact_cols = []
    if gate_on_cdr3_length:
        queries = queries.assign(
            **{_CDR3_LENGTH_COL: queries['cdr3_aa'].str.len().astype(str)})
        ref = ref.assign(**{_CDR3_LENGTH_COL: ref['cdr3_aa'].str.len().astype(str)})
        exact_cols = [_CDR3_LENGTH_COL]

    pairs, skipped = _runSingleChain(queries, ref, compare_cols, gene_cols,
                                     exact_cols, threshold, metric, threads)
    result = _bestHits(pairs, len(queries), ref)
    compared = len(result)
    # A query row with no hit says nothing the metadata's counts do not, and
    # on a real repertoire it is almost every row, so the table holds only
    # the rows that matched something.
    result = result[result['n_hits'].gt(0)].reset_index(drop=True)
    log.info('%d query rows, %d with at least one hit, %d hits within the '
             'threshold in all', compared, len(result), len(pairs))
    metadata['results'] = {
        'query_rows_compared': compared, 'query_rows_with_hits': len(result),
        'query_rows_scored_directly': len(skipped), 'hits': len(pairs)}
    result = _withQueryIds(result, queries)
    matrix = None
    if full_dist_matrix:
        matrix = _withQueryIds(pairs[['query_row', 'reference_id', 'identity']]
                               .rename(columns={'reference_id': 'hit_id'}), queries)
    _writeOutputs(outdir, prefix, result, matrix, metadata)

    return result
