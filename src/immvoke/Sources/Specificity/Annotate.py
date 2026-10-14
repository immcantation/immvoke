"""
Antigen specificity annotation: query-vs-reference receptor comparison

For each row of a query AIRR table, find every row of a fixed reference AIRR
table it agrees with -- one annotation per query row (`query_row`,
`n_hits_total`, `is_unique_hit`, `hit_ids`). This is the third pipeline
stage: download a specificity database (`immvoke specificity iedb download
bcr`), which already lands as a real AIRR rearrangement TSV, then annotate a
user's own AIRR table against it as the reference.

Two pure-Python ways to decide a hit (an exact match is either at 0):

    runFuzzyAnnotation        one or more columns, up to `max_mismatches`
                              substitutions allowed (pooled across columns
                              if more than one), pure Python
    runLevenshteinAnnotation  a single column, up to `max_edits`
                              substitutions/insertions/deletions combined,
                              pure Python (rapidfuzz)

runPairedExactAnnotation builds on an exact match for the case where a
query row is one chain of an antibody rather than a self-contained unit: it
runs it once per chain (heavy, light) and reconciles the two per antibody,
classifying each as agreeing with the same reference antibody on both
chains, one chain only, or neither.

Symscan.py runs the same criteria through symscan's
symmetric-deletion index instead, returning the same result tables; it is
the engine annotateSpecificity() uses. The pure-Python functions here run
only when asked for (`--backend python`), never as a fallback, and are the
reference behaviour Symscan.py is tested against.

annotateSpecificity() is the high level entry point `immvoke annotate`
wraps: it picks among all of the above from a few knobs (sequence_region,
partition_by, gene_match, paired, identity_threshold, distance_metric,
backend) -- see its section at the end of this module. V/J-gene agreement
there goes through Genes.py (`gene_cols`/`gene_match`), which parses gene
calls rather than comparing them as strings; `exact_cols` compares plain
strings.

runLevenshteinAnnotation has neither the wildcard handling nor the
length-bucketed index runFuzzyAnnotation has: Levenshtein distance lets a query
and a candidate differ in length, so length can't prefilter candidates, and a
wildcard has no clean meaning once an alignment can shift position. It also
takes only a single column, unlike runFuzzyAnnotation -- see its section
header for why edit distance does not pool across columns the way Hamming
distance does.

Substitution-matrix-weighted (BLOSUM) annotation is intentionally not
included in this first pass -- see the sibling project this was ported from
(`analysis/matching/matching.py` in specificity-annotations) if it is needed
later.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
from pathlib import Path

import pandas
from rapidfuzz.distance import Levenshtein

# Immvoke imports
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity.Genes import GENE_MATCH_LEVELS, buildGeneGate

log = logging.getLogger(__name__)

#: Reference column read as each hit's identifier when the caller does not
#: say otherwise -- the AIRR rearrangement primary key.
DEFAULT_ID_COL = 'sequence_id'
#: A position where either side has this character never counts as a
#: mismatch, under Hamming distance (runFuzzyAnnotation, and symscan's own
#: Python comparison of such rows).
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


def _valueDistance(a, b):
    """
    Count of substitution mismatches between `a` and `b`.

    Arguments:
      a (str): the first value.
      b (str): the second value.

    Returns:
      int: mismatch count, or None if the lengths differ -- this only models
      substitutions, not insertions/deletions. A WILDCARD on either side of a
      position is never a mismatch.
    """
    if len(a) != len(b):
        return None

    return sum(ca != cb and ca != WILDCARD and cb != WILDCARD
              for ca, cb in zip(a, b))


def _rowDistance(a, b):
    """
    Total substitution mismatches summed across every comparison column.

    Arguments:
      a (tuple): one row's values, one per comparison column.
      b (tuple): the row it is being compared to, same shape.

    Returns:
      int: the summed distance, or None if any column pair differs in length.
    """
    total = 0
    for x, y in zip(a, b):
        distance = _valueDistance(x, y)
        if distance is None:
            return None
        total += distance

    return total


def _mismatchPositions(a, b):
    """
    Every substitution mismatch between `a` and `b`.

    Arguments:
      a (tuple): one row's values, one per comparison column.
      b (tuple): the row it is being compared to, same shape.

    Returns:
      list: (column index, 0-based position) pairs, same wildcard rule as
      _valueDistance.
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


def _bucketKey(row_values, exact_values):
    """
    Build the index bucket a row belongs to.

    Arguments:
      row_values (tuple): the row's comparison column values.
      exact_values (tuple): the row's exact column values, empty if none.

    Returns:
      tuple: exact column values followed by the length of each comparison
      column value -- substitutions can't change length, so two rows can
      only agree when every comparison column is the same length.
    """
    return exact_values + tuple(len(value) for value in row_values)


def buildAnnotationIndex(reference, compare_cols, id_col=DEFAULT_ID_COL,
                         exact_cols=None):
    """
    Bucket a reference table for fast lookup.

    Arguments:
      reference (pandas.DataFrame): the fixed table being annotated against.
      compare_cols (list): column names to compare.
      id_col (str): column holding each reference row's identifier.
      exact_cols (list): column names that must agree exactly before
        `compare_cols` are even compared, e.g. `['v_call']` to only ever
        compare CDR3s within the same V-gene. None for no such gate.

    Returns:
      dict: bucket key (see _bucketKey) to a list of (row_values, id,
      reference position) tuples.
    """
    values = _rowValues(reference, compare_cols)
    exact_values = (_rowValues(reference, exact_cols) if exact_cols
                    else [()] * len(reference))
    ids = reference[id_col].tolist()

    index = {}
    for position, (row_values, row_exact, row_id) in enumerate(
            zip(values, exact_values, ids)):
        index.setdefault(_bucketKey(row_values, row_exact), []).append(
            (row_values, row_id, position))

    return index


def runFuzzyAnnotation(queries, reference, compare_cols, max_mismatches,
                       id_col=DEFAULT_ID_COL, exact_cols=None,
                       threshold_is_fraction=False, gene_cols=None,
                       gene_match='gene'):
    """
    Annotate each query row against `reference` allowing up to
    `max_mismatches` substitutions across `compare_cols` (pooled if more
    than one) instead of requiring exact agreement.

    `exact_cols` matters beyond speed: an unconstrained CDR3-only fuzzy
    annotation is prone to coincidental hits (especially on short CDR3s,
    where 1-2 substitutions is a large fraction of the sequence) that don't
    share a V-gene with the query and are unlikely to be a real related
    receptor.

    Arguments:
      queries (pandas.DataFrame): the input rows to annotate.
      reference (pandas.DataFrame): the fixed table being annotated against.
      compare_cols (list): column names to compare.
      max_mismatches (float): the mismatch tolerance -- an absolute count, or
        a fraction of the query's length when `threshold_is_fraction` is True.
      id_col (str): column in `reference` holding each row's identifier.
      exact_cols (list): column names that must agree exactly first, e.g.
        `['v_call']`. None for no such gate.
      threshold_is_fraction (bool): treat `max_mismatches` as a fraction of
        the query's summed comparison column length (e.g. 0.15 allows up to
        15% mismatches) instead of a fixed count. Resolved to an absolute
        count once per query row, since every candidate in its bucket
        already shares its exact length.
      gene_cols (list): gene call columns that must agree at `gene_match`
        level first, e.g. `['v_call', 'j_call']`; blank or unrecognized
        calls agree with anything (see Genes.py). None for no such gate.
      gene_match (str): 'allele', 'gene' or 'family'; see Genes.py.

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `min_mismatches`, `is_unique_hit`, `hit_ids`,
      `mismatch_positions` (per id in `hit_ids`, same order; column
      prefixed when `compare_cols` has more than one column).
    """
    index = buildAnnotationIndex(reference, compare_cols, id_col, exact_cols)
    gate = buildGeneGate(queries, reference, gene_cols, gene_match)
    values = _rowValues(queries, compare_cols)
    exact_values = (_rowValues(queries, exact_cols) if exact_cols
                    else [()] * len(values))

    rows = []
    for query_row, (row_values, row_exact) in enumerate(zip(values, exact_values)):
        candidates = index.get(_bucketKey(row_values, row_exact), [])
        max_distance = (_allowedDifferences(max_mismatches,
                                            sum(len(v) for v in row_values))
                        if threshold_is_fraction else max_mismatches)

        best_by_id = {}
        for cand_values, row_id, ref_row in candidates:
            if gate is not None and not gate.agrees(query_row, ref_row):
                continue
            distance = _rowDistance(row_values, cand_values)
            if distance is not None and distance <= max_distance:
                if row_id not in best_by_id or distance < best_by_id[row_id][0]:
                    best_by_id[row_id] = (distance, cand_values)

        hit_ids = sorted(best_by_id, key=lambda x: (best_by_id[x][0], str(x)))
        mismatch_strs = [
            _formatPositions(_mismatchPositions(row_values, best_by_id[x][1]),
                            len(compare_cols))
            for x in hit_ids
        ]
        rows.append({
            'query_row': query_row,
            'n_hits_total': len(hit_ids),
            # float('nan'), not None: keeps the column a uniform float64
            # rather than a mixed-type object column when some query rows
            # have no hit at all.
            'min_mismatches': float(min(d for d, _ in best_by_id.values()))
                              if best_by_id else float('nan'),
            'is_unique_hit': len(hit_ids) <= 1,
            'hit_ids': ';'.join(str(x) for x in hit_ids),
            'mismatch_positions': ';'.join(mismatch_strs),
        })

    return pandas.DataFrame(rows)


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
                             threads=1, backend='symscan', gene_match='gene'):
    """
    Classify each query antibody by whether the *same* reference antibody
    agrees with it on the heavy chain, the light chain, both, or neither.

    An antibody is two rows of an AIRR table -- a heavy chain and, where
    present, a light chain -- sharing one `cell_id` and told apart by
    `locus` (`IGH` for heavy, `IGK`/`IGL` for light). Agreement on a chain
    is the same V-gene and J-gene (at `gene_match` level, see Genes.py) and
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
      backend (str): 'symscan' (Symscan.runSymscanAnnotation at 0
        substitutions) or 'python' (runFuzzyAnnotation at 0 substitutions)
        -- the engine behind each chain's exact match; the result is the
        same either way.
      gene_match (str): 'allele', 'gene' or 'family' -- how V/J-gene calls
        are compared; blank or unrecognized calls agree with anything.

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
    if backend == 'symscan':
        # Deferred: Symscan builds on this module's helpers.
        from immvoke.Sources.Specificity.Symscan import runSymscanAnnotation

        def _exact(query_chain, reference_chain):
            return runSymscanAnnotation(query_chain, reference_chain, [cdr3_col], 0,
                                        id_col=CELL_ID_COL, gene_cols=gene_cols,
                                        gene_match=gene_match, threads=threads)
    elif backend == 'python':
        def _exact(query_chain, reference_chain):
            return runFuzzyAnnotation(query_chain, reference_chain, [cdr3_col], 0,
                                      id_col=CELL_ID_COL, gene_cols=gene_cols,
                                      gene_match=gene_match)
    else:
        raise AnnotationError("unknown backend %r for paired annotation; choose "
                              "'symscan' or 'python'" % (backend,))

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


# ── indel-tolerant annotation (Levenshtein) ─────────────────────────────────
#
# Unlike runFuzzyAnnotation, there is no length-bucketed index and no
# 'X'-wildcard handling here: Levenshtein allows a query and a candidate to
# differ in length, so a per-column length key can't prefilter candidates, and
# a wildcard has no clean meaning once an alignment can shift positions.
# Candidates are only pre-filtered by `exact_cols` (e.g. V-gene, J-gene), then
# compared to every query brute-force within that block. Also unlike
# runFuzzyAnnotation, this only ever takes a single seq_col: edit distance
# summed across multiple columns independently aligned doesn't correspond to
# any one alignment of the concatenated sequence, so it does not generalize
# the way Hamming distance pools cleanly across columns.

def runLevenshteinAnnotation(queries, reference, seq_col, max_edits,
                             id_col=DEFAULT_ID_COL, exact_cols=None,
                             threshold_is_fraction=False, gene_cols=None,
                             gene_match='gene'):
    """
    Annotate each query row against `reference` allowing up to `max_edits`
    substitutions, insertions and deletions combined in `seq_col` -- the
    indel-tolerant counterpart to runFuzzyAnnotation's substitution-only
    criterion, via rapidfuzz's Levenshtein distance.

    Arguments:
      queries (pandas.DataFrame): the input rows to annotate.
      reference (pandas.DataFrame): the fixed table being annotated against.
      seq_col (str): the single column to compare.
      max_edits (float): the edit-distance tolerance (substitutions +
        insertions + deletions combined) -- an absolute count, or a fraction
        of the query's length when `threshold_is_fraction` is True.
      id_col (str): column in `reference` holding each row's identifier.
      exact_cols (list): column names that must agree exactly first, e.g.
        `['v_call', 'j_call']`. None for no such gate.
      threshold_is_fraction (bool): treat `max_edits` as a fraction of the
        query's `seq_col` length (e.g. 0.1 allows up to 10% edits) instead
        of a fixed count, resolved to a whole count once per query row.
      gene_cols (list): see runFuzzyAnnotation.
      gene_match (str): see runFuzzyAnnotation.

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `min_mismatches` (edit distance here, not substitution
      count), `is_unique_hit`, `hit_ids`.
    """
    reference_exact = (_rowValues(reference, exact_cols) if exact_cols
                       else [()] * len(reference))
    reference_seqs = reference[seq_col].fillna('').astype(str).tolist()
    reference_ids = reference[id_col].tolist()

    buckets = {}
    for ref_row, (seq, exact_values, row_id) in enumerate(
            zip(reference_seqs, reference_exact, reference_ids)):
        buckets.setdefault(exact_values, []).append((seq, row_id, ref_row))
    gate = buildGeneGate(queries, reference, gene_cols, gene_match)

    query_exact = (_rowValues(queries, exact_cols) if exact_cols
                  else [()] * len(queries))
    query_seqs = queries[seq_col].fillna('').astype(str).tolist()

    rows = []
    for query_row, (query_seq, row_exact) in enumerate(zip(query_seqs, query_exact)):
        limit = (_allowedDifferences(max_edits, len(query_seq))
                 if threshold_is_fraction else int(max_edits))
        best_by_id = {}
        for cand_seq, row_id, ref_row in buckets.get(row_exact, []):
            if gate is not None and not gate.agrees(query_row, ref_row):
                continue
            distance = Levenshtein.distance(query_seq, cand_seq,
                                            score_cutoff=limit)
            if distance <= limit:
                if row_id not in best_by_id or distance < best_by_id[row_id]:
                    best_by_id[row_id] = distance

        hit_ids = sorted(best_by_id, key=lambda x: (best_by_id[x], str(x)))
        rows.append({
            'query_row': query_row,
            'n_hits_total': len(hit_ids),
            'min_mismatches': float(min(best_by_id.values()))
                              if best_by_id else float('nan'),
            'is_unique_hit': len(hit_ids) <= 1,
            'hit_ids': ';'.join(str(x) for x in hit_ids),
        })

    return pandas.DataFrame(rows)


# ── high level API: annotateSpecificity() ───────────────────────────────────
#
# A file-or-DataFrame-in, DataFrame-and-file-out entry point over the
# primitives above, for a caller who wants "match this receptor table against
# that reference table" without choosing which run*Annotation to call. It
# picks one from a small set of knobs (sequence_region, partition_by,
# identity_threshold, distance_metric, backend) -- the flow `immvoke annotate`
# exposes on the commandline:
#
#     paired  -> runPairedExactAnnotation (exact V, J, CDR3 on both chains)
#     else    -> heavy chain only -> sequence_region columns
#             -> partition_by gate (none, V, VJ, VJ + CDR3 length)
#             -> identity_threshold: 1.0 exact, below it Hamming/Levenshtein
#             -> backend: symscan (a required dependency), or Python only
#                when asked for -- never as a fallback
#                (symscan itself compares in Python what it cannot index,
#                e.g. a distance above Symscan.MAX_DISTANCE)

#: sequence_region -> AIRR column(s) to compare. 'CDRs' pools CDR1-3, so it
#: only runs under Hamming distance (or an exact match), never Levenshtein --
#: see the Levenshtein section note on why edit distance does not pool.
SEQUENCE_REGION_COLUMNS = {
    'cdr3': ('cdr3_aa',),
    'vdj': ('sequence_aa',),
    'cdrs': ('cdr1_aa', 'cdr2_aa', 'cdr3_aa'),
}

#: partition_by -> gene call columns two rows must agree on (at gene_match
#: level, see Genes.py) before their sequence counts as a hit. 'VJL' also
#: gates on exact CDR3 length; see annotateSpecificity.
PARTITION_GENE_COLUMNS = {
    'v': ('v_call',),
    'vj': ('v_call', 'j_call'),
    'vjl': ('v_call', 'j_call'),
}

#: distance_metric values accepted ('levenstein' is a common misspelling).
DISTANCE_METRICS = {'hamming': 'hamming', 'levenshtein': 'levenshtein',
                    'levenstein': 'levenshtein'}

#: backend values accepted. symscan, a required dependency, is the default;
#: python runs only when asked for, never as a fallback.
BACKENDS = ('symscan', 'python')

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


def _resolveBackend(backend):
    """
    Validate `backend`, the engine that will run.

    Arguments:
      backend (str): one of BACKENDS, case insensitive.

    Returns:
      str: 'symscan' or 'python'.

    Raises:
      AnnotationError: if `backend` is not recognized.
    """
    name = str(backend).strip().lower()
    if name not in BACKENDS:
        raise AnnotationError("unknown backend %r; choose one of: %s"
                              % (backend, ', '.join(BACKENDS)))
    log.info('matching backend: %s', name)

    return name


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
                    gene_match, threshold, metric, backend, threads):
    """
    Run one single-chain comparison on `backend`.

    Arguments:
      queries (pandas.DataFrame): the query rows, already filtered.
      reference (pandas.DataFrame): the reference rows, already filtered.
      compare_cols (list): the column(s) to compare.
      gene_cols (list): gene call columns that must agree at `gene_match`
        level, possibly empty.
      exact_cols (list): columns that must agree exactly, possibly empty.
      gene_match (str): 'allele', 'gene' or 'family'.
      threshold (float): identity in (0.0, 1.0]; 1.0 is an exact match.
      metric (str): 'hamming' or 'levenshtein'.
      backend (str): 'symscan' or 'python'.
      threads (int): symscan's thread pool; see annotateSpecificity.

    Returns:
      pandas.DataFrame: one row per query row, see runFuzzyAnnotation.
    """
    exact = threshold >= 1.0
    # An exact match is 0 substitutions under either metric -- the same
    # thing, and the only way to run exact agreement on pooled 'CDRs'.
    distance_type = 'hamming' if exact else metric
    tolerance = 0 if exact else 1 - threshold
    gates = dict(id_col=DEFAULT_ID_COL, exact_cols=exact_cols or None,
                 gene_cols=gene_cols or None, gene_match=gene_match)

    if backend == 'symscan':
        # Deferred: Symscan builds on this module's helpers.
        from immvoke.Sources.Specificity.Symscan import runSymscanAnnotation

        return runSymscanAnnotation(queries, reference, compare_cols, tolerance,
                                    distance_type=distance_type,
                                    threshold_is_fraction=not exact,
                                    threads=threads, **gates)

    if distance_type == 'hamming':
        return runFuzzyAnnotation(queries, reference, compare_cols, tolerance,
                                  threshold_is_fraction=not exact, **gates)

    return runLevenshteinAnnotation(queries, reference, compare_cols[0], tolerance,
                                    threshold_is_fraction=True, **gates)


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
                        distance_metric='hamming', gene_match='gene', threads=1,
                        output=None, backend='symscan'):
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
      gene_match (str): how `partition_by` and `paired` compare v_call/
        j_call: 'allele' (same gene and allele), 'gene' (the default: same
        gene, any allele) or 'family'. A row with several calls agrees if
        any one does; a blank or unrecognized call is unknown and agrees
        with anything. See Genes.py.
      identity_threshold (float): minimum sequence identity in (0.0, 1.0];
        None or 1.0 requires an exact match, a lower value allows up to
        floor((1 - identity_threshold) * query length) differences. Ignored
        when `paired` is True.
      distance_metric (str): 'hamming' (substitutions only) or 'levenshtein'
        (substitutions, insertions and deletions; a single-column
        `sequence_region` only, i.e. not 'CDRs'). Only matters when
        `identity_threshold` is below 1.0.
      threads (int): symscan's thread pool size (-1 for every core); only
        the first symscan run in a process sets it. The python backend is
        single threaded.
      output (str): if given, also write the result here as a TSV.
      backend (str): 'symscan' (the default) or 'python' (only when asked
        for; never used as a fallback). Either gives the same result;
        symscan compares in Python whatever it cannot index (an 'X'
        wildcard, a distance above 14, very long sequences at a high
        tolerance).

    Returns:
      pandas.DataFrame: when `paired`, one row per query `cell_id` (see
      runPairedExactAnnotation); otherwise one row per heavy chain query
      row, every one of them -- 0 hits if unmatched -- with `query_row`,
      `n_hits_total`, `min_mismatches`, `is_unique_hit`, `hit_ids` (plus
      `mismatch_positions` under Hamming).

    Raises:
      AnnotationError: if the reference is not an AIRR rearrangement table,
        an argument is not recognized, a needed column is missing, no heavy
        chain row is left to annotate, or an argument combination cannot
        run (e.g. 'CDRs' under 'levenshtein').
    """
    if gene_match not in GENE_MATCH_LEVELS:
        raise AnnotationError('unknown gene_match %r; choose one of: %s'
                              % (gene_match, ', '.join(GENE_MATCH_LEVELS)))
    backend = _resolveBackend(backend)
    queries = _readAirrTsv(query)
    ref = _readAirrTsv(reference)
    _requireAirr(ref, 'reference')

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
        result = runPairedExactAnnotation(queries, ref, threads=threads,
                                          backend=backend, gene_match=gene_match)
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
                             gene_match, threshold, metric, backend, threads)
    hit = int(result['n_hits_total'].gt(0).sum())
    log.info('%d query rows, %d with at least one hit, %d unique (0 or 1 hit)',
             len(result), hit, int(result['is_unique_hit'].sum()))
    _writeIfRequested(result, output)

    return result
