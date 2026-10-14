"""
Antigen specificity annotation: query-vs-reference receptor comparison

For each row of a query AIRR table, find every row of a fixed reference AIRR
table it agrees with -- one annotation per query row (`query_row`,
`n_hits_total`, `is_unique_hit`, `hit_ids`). This is the last pipeline stage:
download a specificity database's table, which lands as an AIRR
rearrangement TSV (`immvoke specificity bcr`), then annotate a user's own
AIRR table against it as the reference.

annotateSpecificity() is the entry point `immvoke annotate` wraps: it picks
the comparison from a few knobs (sequence_region, partition_by,
paired, identity_threshold, distance_metric) -- see its section at the end of
this module -- and runs it through Symscan.runSymscanAnnotation, symscan's
symmetric-deletion index. There is no other matching engine: rows too large
for symscan are skipped with a warning rather than compared another way.

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
from pathlib import Path

import pandas

# Immvoke imports
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity.Genes import inferLocus

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
    return list(frame[columns].fillna('').itertuples(index=False, name=None))


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
                             threads=1):
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

    Returns:
      pandas.DataFrame: one row per query `cell_id`: `classification` (one
      of `both`, `heavy_only`, `light_only`, `conflicting`, `neither`),
      `heavy_hit_ids`, `light_hit_ids` (every reference antibody agreeing on
      that chain alone, semicolon joined, sorted) and `both_hit_ids` (their
      intersection -- non-empty only when `classification` is `both`).
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

    rows = []
    for cell_id in sorted(set(queries[CELL_ID_COL])):
        heavy = heavy_hits.get(cell_id, set())
        light = light_hits.get(cell_id, set())
        both = heavy & light

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

    return pandas.DataFrame(rows)


# ── high level API: annotateSpecificity() ───────────────────────────────────
#
# A file-or-DataFrame-in, DataFrame-and-file-out entry point, for a caller
# who wants "match this receptor table against that reference table". It
# picks the comparison from a small set of knobs (sequence_region,
# partition_by, identity_threshold, distance_metric) -- the flow
# `immvoke annotate` exposes on the commandline:
#
#     paired  -> runPairedExactAnnotation (V, J, exact CDR3 on both chains)
#     else    -> heavy chain only -> sequence_region columns
#             -> partition_by gate (none, V, VJ, VJ + CDR3 length)
#             -> identity_threshold: 1.0 exact, below it Hamming/Levenshtein
#             -> symscan; rows beyond its limits are skipped, with a warning

#: sequence_region -> AIRR column(s) to compare. 'CDRs' pools CDR1-3, so it
#: only runs under Hamming distance (or an exact match), never Levenshtein:
#: edit distance summed over independently aligned columns is not the edit
#: distance of any one alignment, the way Hamming distance pools cleanly.
SEQUENCE_REGION_COLUMNS = {
    'cdr3': ('cdr3_aa',),
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
      sequence_region (str): 'CDR3', 'VDJ' or 'CDRs', case insensitive.

    Returns:
      list: the AIRR column(s) to compare.

    Raises:
      AnnotationError: if `sequence_region` is not recognized.
    """
    columns = SEQUENCE_REGION_COLUMNS.get(str(sequence_region).strip().lower())
    if columns is None:
        raise AnnotationError("unknown sequence_region %r; choose one of: CDR3, "
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


def _fillLocus(frame, label):
    """
    Fill a blank `locus` from the row's own gene calls (Genes.inferLocus):
    the V call's locus, else the J call's, else left blank.

    Arguments:
      frame (pandas.DataFrame): the table to fill.
      label (str): 'query' or 'reference', for the log line only.

    Returns:
      pandas.DataFrame: `frame`, with blank loci filled where its gene
      calls name one -- unchanged if it has no locus column or no blank.
    """
    if LOCUS_COL not in frame.columns:
        return frame
    blank = frame[LOCUS_COL].fillna('') == ''
    if not blank.any():
        return frame

    v_calls = frame['v_call'] if 'v_call' in frame.columns else pandas.Series('', frame.index)
    j_calls = frame['j_call'] if 'j_call' in frame.columns else pandas.Series('', frame.index)
    inferred = pandas.Series([inferLocus(v, j) for v, j in zip(v_calls[blank], j_calls[blank])],
                             index=frame.index[blank])
    filled = inferred[inferred != '']
    if len(filled):
        frame = frame.copy()
        frame.loc[filled.index, LOCUS_COL] = filled
    log.info('%s: %d rows had a blank locus; filled from their gene calls: %s; '
             '%d left blank', label, int(blank.sum()),
             ', '.join('%s=%d' % kv for kv in sorted(filled.value_counts().items()))
             or 'none', int(blank.sum()) - len(filled))

    return frame


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
    Run one single-chain comparison through symscan.

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
      pandas.DataFrame: one row per query row, see
      Symscan.runSymscanAnnotation.
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

    return runSymscanAnnotation(queries, reference, compare_cols, tolerance,
                                distance_type=distance_type,
                                threshold_is_fraction=not exact,
                                threads=threads, **gates)


def _writeIfRequested(result, output):
    """
    Write `result` to `output` as a TSV, if a path was given.

    Arguments:
      result (pandas.DataFrame): the annotation result.
      output (str): where to write it, or None to skip writing.
    """
    if output is None:
        return

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, sep='\t', index=False)
    log.info('wrote %s', out)


def annotateSpecificity(query, reference, sequence_region='CDR3',
                        partition_by=None, paired=False, identity_threshold=None,
                        distance_metric='hamming', threads=1,
                        output=None):
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
      sequence_region (str): which sequence to compare: 'CDR3' (cdr3_aa),
        'VDJ' (sequence_aa, the full V(D)J amino acid sequence) or 'CDRs'
        (cdr1_aa, cdr2_aa and cdr3_aa pooled -- Hamming only). Ignored when
        `paired` is True.
      partition_by (str): only compare row pairs in the same partition: None
        (no restriction), 'V' (same V-gene), 'VJ' (same V- and J-gene) or
        'VJL' (same V-gene, J-gene and CDR3 length -- needs `cdr3_aa` in both
        tables whatever `sequence_region` is). Ignored when `paired` is True.
      paired (bool): classify each query antibody (`cell_id`) by V-gene,
        J-gene and exact CDR3 agreement on both its heavy and light chain
        with the same reference antibody (runPairedExactAnnotation). When
        False, only heavy chain (locus IGH) rows are compared; a table
        without a locus column is compared whole.
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
      output (str): if given, also write the result here as a TSV.

    Returns:
      pandas.DataFrame: when `paired`, one row per query `cell_id` (see
      runPairedExactAnnotation); otherwise one row per heavy chain query
      row, every one of them -- 0 hits if unmatched -- with `query_row`,
      `n_hits_total`, `min_mismatches`, `is_unique_hit`, `hit_ids` (plus
      `mismatch_positions` under Hamming).

    Raises:
      AnnotationError: if the reference is not an AIRR rearrangement table,
        an argument is not recognized, a needed column is missing, no heavy
        chain row is left to annotate, an argument combination cannot run
        (e.g. 'CDRs' under 'levenshtein'). Rows beyond symscan's limits are
        skipped with a warning instead (see Symscan.py).
    """
    queries = _readAirrTsv(query)
    ref = _readAirrTsv(reference)
    _requireAirr(ref, 'reference')
    queries = _fillLocus(queries, 'query')
    ref = _fillLocus(ref, 'reference')

    if paired:
        if partition_by is not None or identity_threshold is not None:
            log.warning('paired annotation always compares V-gene, J-gene and '
                        'exact CDR3; partition_by and identity_threshold are '
                        'ignored')
        paired_cols = [CELL_ID_COL, LOCUS_COL, 'cdr3_aa', 'v_call', 'j_call']
        _requireColumns(queries, paired_cols, 'query', 'paired annotation')
        _requireColumns(ref, paired_cols, 'reference', 'paired annotation')
        chain_loci = (HEAVY_LOCUS,) + LIGHT_LOCI
        for frame, label in ((queries, 'query'), (ref, 'reference')):
            if not frame[LOCUS_COL].isin(chain_loci).any():
                raise AnnotationError('%s has no antibody chain rows (locus %s) '
                                      'to pair' % (label, ', '.join(chain_loci)))
        result = runPairedExactAnnotation(queries, ref, threads=threads)
        log.info('%d query antibodies: %s', len(result),
                 ', '.join('%s=%d' % (label, count) for label, count
                          in result['classification'].value_counts().items()))
        _writeIfRequested(result, output)

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

    result = _runSingleChain(queries, ref, compare_cols, gene_cols, exact_cols,
                             threshold, metric, threads)
    hit = int(result['n_hits_total'].gt(0).sum())
    log.info('%d query rows, %d with at least one hit, %d unique (0 or 1 hit)',
             len(result), hit, int(result['is_unique_hit'].sum()))
    _writeIfRequested(result, output)

    return result
