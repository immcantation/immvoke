"""
symscan annotation: exact, Hamming and Levenshtein

    runSymscanAnnotation  up to `max_distance` substitutions (distance_type
                          'hamming', one or more pooled columns) or
                          substitutions/insertions/deletions combined
                          ('levenshtein', a single column); an exact match is
                          `max_distance` 0

symscan (https://github.com/yutanagano/symscan) builds a Rust
symmetric-deletion index and only ever looks at likely matches, rather than
scoring every query x reference pair -- so it stays fast as the reference
grows. It is the only engine that finds hits: no hit is searched for in
Python. (Distances.py scores every gated pair for a threshold run's pairs
file; the hits still come from here.)

Blanks. A row with a blank value in any compared column has nothing to
compare, so it never matches: a blank query gets no hit, and a blank
reference is never hit.

Gates. symscan has no gate of its own: `exact_cols` (exact string
agreement), each column's length for a pooled multi-column Hamming
comparison, and `gene_cols` (V/J-gene agreement, same gene and any allele,
see Genes.py) are applied to the pairs it returns.

`gene_cols` is also applied before the search, as a partition: the
comparison runs once per gene key rather than once over everything (see
_geneGroups). The index cost grows combinatorially with the tolerance, and
in one index every query collides with references of every gene
combination -- collisions whose distance is computed only for the gate to
discard them. Grouping removes that work instead of paying for it: on a
21,412 row repertoire against IEDB at 0.7 identity, 57 seconds became 3.
The gate still runs afterwards and still decides what is kept, so the
answer does not depend on the partition being clever.

The "X" wildcard. An X is an unresolved residue: under 'hamming' (and so
under exact matching) a position where either side has one never counts as a
difference. symscan has no wildcard and compares X as an ordinary letter,
so the X's are taken out of the comparison instead: references are grouped
by where their X's sit, and those positions are cut out of the group's
references and of every query of the same length before symscan sees them
-- a reference X then never counts, and the search runs at the real limit.
A query's own X's, rarer, are covered by slack (its limit raised by its X
count), and every returned pair involving an X is re-scored under the
wildcard rule, keeping only those within the real limit. Levenshtein has no
wildcard rule: an X there is an ordinary letter.

Limits. symscan's index holds every deletion variant of every string on both
sides -- sum(C(length, k) for k <= distance) per string -- and it aborts the
whole process (an uncatchable Rust panic or failed allocation) once that
outgrows memory, or once the distance exceeds MAX_DISTANCE. Each call is
therefore kept within MAX_VARIANTS by grouping queries by their limit (and,
under 'hamming', by length, since only same-length references can match)
and chunking both the queries and the references of a group to half the
budget each. What cannot fit at all -- a limit above MAX_DISTANCE, X slack
included, or one sequence whose index alone exceeds half the budget (a very
long sequence at a loose threshold) -- is skipped, with a warning naming the
query rows that may therefore miss hits.

symscan parallelizes via rayon's global pool, which reads RAYON_NUM_THREADS
once, on first use -- so `threads` only takes effect on the first call in a
process; later calls reuse that pool.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import itertools
import logging
import math
import os
from functools import cache

import numpy
import pandas
import symscan

# Immvoke imports
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity.Annotate import (
    DEFAULT_ID_COL,
    WILDCARD,
    _allowedDifferences,
    _exactValues,
    _formatPositions,
    _mismatchPositions,
    _rowValues,
)
from immvoke.Sources.Specificity.Genes import buildGeneGate, parseCall

log = logging.getLogger(__name__)

#: Largest max_distance symscan 0.8 accepts without panicking.
MAX_DISTANCE = 14

#: Most deletion variants, both sides together, one symscan call may index.
#: symscan's peak memory measured at roughly 13-40 bytes per variant, so this
#: keeps a call to around 1-2 GB.
MAX_VARIANTS = 40_000_000

#: RAYON_NUM_THREADS the rayon pool was started with, once one has been.
_pool_threads = None


@cache
def _deletionVariants(length, distance):
    """
    Arguments:
      length (int): a string's length.
      distance (int): the distance it is indexed at.

    Returns:
      int: how many deletion variants symscan's index holds for it.
    """
    return sum(math.comb(length, k) for k in range(min(distance, length) + 1))


def _resolveLimits(query_values, max_distance, threshold_is_fraction):
    """
    Resolve the tolerance to an absolute distance per query row.

    Arguments:
      query_values (list): per-query comparison column tuples.
      max_distance (float): an absolute count, or a fraction of the query's
        summed comparison column length when `threshold_is_fraction` is True.
      threshold_is_fraction (bool): see runSymscanAnnotation.

    Returns:
      numpy.ndarray: one int per query row.
    """
    if threshold_is_fraction:
        return numpy.array([_allowedDifferences(max_distance,
                                                sum(len(v) for v in values))
                            for values in query_values], dtype=numpy.int64)

    return numpy.full(len(query_values), int(max_distance), dtype=numpy.int64)


def _startPool(threads):
    """
    Size rayon's global pool before symscan's first call in this process.

    Arguments:
      threads (int): worker threads, -1 for every available core.
    """
    global _pool_threads
    wanted = os.cpu_count() if threads == -1 else threads
    if _pool_threads is None:
        os.environ['RAYON_NUM_THREADS'] = str(wanted)
        _pool_threads = wanted
    elif wanted != _pool_threads:
        log.debug('symscan thread pool already started with %d threads; '
                  'threads=%d ignored', _pool_threads, threads)


def _keyCodes(query_keys, reference_keys):
    """
    Encode each row's gate key as an integer, consistently across both
    sides, so candidate pairs can be filtered with one vectorized compare.

    Arguments:
      query_keys (list): one tuple per query row.
      reference_keys (list): one tuple per reference row.

    Returns:
      tuple: (query codes, reference codes), numpy int arrays.
    """
    joined = ['\x1f'.join(map(str, key)) for key in query_keys + reference_keys]
    codes, _ = pandas.factorize(pandas.Series(joined, dtype=object))

    return codes[:len(query_keys)], codes[len(query_keys):]


def _without(seq, positions):
    """
    Returns:
      str: `seq` with the characters at `positions` (a set) cut out.
    """
    return ''.join(c for i, c in enumerate(seq) if i not in positions)


def _describeRows(rows, queries, shown=5):
    """
    Returns:
      str: the first few query rows, for a log line -- by their own
      `sequence_id` where the query has one, else by position.
    """
    more = ', ...' if len(rows) > shown else ''
    if DEFAULT_ID_COL in queries.columns:
        ids = queries[DEFAULT_ID_COL].iloc[rows[:shown]].astype(str)
        return '%s %s%s' % (DEFAULT_ID_COL, ', '.join(ids), more)

    return 'query_row %s%s' % (', '.join(str(int(r)) for r in rows[:shown]), more)


def _concat(parts):
    """
    Arguments:
      parts (list): (query positions, reference positions, distances) tuples.

    Returns:
      tuple: the three arrays, each concatenated across `parts`.
    """
    if not parts:
        return tuple(numpy.array([], dtype=numpy.int64) for _ in range(3))

    return tuple(numpy.concatenate(x).astype(numpy.int64) for x in zip(*parts))


def _chunks(rows, costs, budget):
    """
    Split rows into consecutive chunks whose summed cost stays within budget.

    Arguments:
      rows (numpy.ndarray): positions, each costing no more than `budget`.
      costs (list): each row's deletion-variant count, same order.
      budget (int): the most one chunk may cost.

    Returns:
      list: numpy int arrays of positions.
    """
    chunks, chunk, total = [], [], 0
    for row, cost in zip(rows.tolist(), costs):
        if chunk and total + cost > budget:
            chunks.append(chunk)
            chunk, total = [], 0
        chunk.append(row)
        total += cost
    if chunk:
        chunks.append(chunk)

    return [numpy.array(c, dtype=numpy.int64) for c in chunks]


def _geneKeys(frame, gene_cols):
    """
    Each row's gene keys: one per combination of the genes it calls, so a
    row calling two V genes and one J gene gets two keys.

    Two rows pass the gene gate when they name a gene in common in every
    gene column (Genes.callsAgree is set intersection), which is exactly
    when their key sets intersect -- so grouping by key never separates a
    pair the gate would keep. A row with a blank or unrecognized call in
    any column gets no key at all, matching the gate's rule that such a
    call agrees with nothing.

    Arguments:
      frame (pandas.DataFrame): the rows.
      gene_cols (list): gene call columns, e.g. ['v_call', 'j_call'].

    Returns:
      list: one tuple of keys per row, each key a tuple of one gene per
      column.
    """
    parsed = []
    for col in gene_cols:
        values = frame[col].astype(object).fillna('').astype(str).tolist()
        # Each distinct call parsed once: a repertoire repeats its calls.
        cache = {value: sorted(parseCall(value)) for value in set(values)}
        parsed.append([cache[value] for value in values])

    return [tuple(itertools.product(*genes)) for genes in zip(*parsed)]


def _geneGroups(queries, reference, gene_cols, query_rows, ref_rows):
    """
    Split the comparison into one group per gene key, so symscan indexes
    only the rows that could pass the gene gate.

    symscan's cost grows combinatorially with the tolerance, and in one
    index every query collides with references of every gene combination --
    collisions the gate then throws away, after their distance has been
    computed. Grouping first removes that work rather than paying for it.
    The gate still runs afterwards: grouping decides what is compared, the
    gate decides what is kept.

    Arguments:
      queries (pandas.DataFrame): the query rows.
      reference (pandas.DataFrame): the reference rows.
      gene_cols (list): gene call columns to group on, or None for no
        grouping -- then the whole comparison is one group.
      query_rows (numpy.ndarray): the query positions to compare.
      ref_rows (numpy.ndarray): the reference positions to compare.

    Yields:
      tuple: (query positions, reference positions) per group, only where
      both sides have rows.
    """
    if not gene_cols:
        yield query_rows, ref_rows
        return

    ref_keys = _geneKeys(reference, gene_cols)
    by_key = {}
    for ref in ref_rows.tolist():
        for key in ref_keys[ref]:
            by_key.setdefault(key, []).append(ref)

    query_keys = _geneKeys(queries, gene_cols)
    wanted = {}
    for query in query_rows.tolist():
        for key in query_keys[query]:
            if key in by_key:
                wanted.setdefault(key, []).append(query)

    for key, members in wanted.items():
        yield (numpy.array(members, dtype=numpy.int64),
               numpy.array(by_key[key], dtype=numpy.int64))


def _searchGroup(query_seqs, ref_seqs, query_lengths, limits, query_x, ref_x,
                 distance_type, query_rows, ref_rows):
    """
    Run one group's comparison through symscan, handling the X wildcard.

    References are grouped by where their X's sit, and those positions are
    cut out of the group's references and of every same-length query, so a
    reference X never counts; a query's own X's are covered by slack
    instead (see the module note).

    Arguments:
      query_seqs (list): every query's compared sequence.
      ref_seqs (list): every reference's.
      query_lengths (numpy.ndarray): len() of each query sequence.
      limits (numpy.ndarray): each query's allowed differences.
      query_x (numpy.ndarray): each query's X count.
      ref_x (numpy.ndarray): each reference's X count.
      distance_type (str): 'hamming' or 'levenshtein'.
      query_rows (numpy.ndarray): the query positions to compare.
      ref_rows (numpy.ndarray): the reference positions to compare.

    Returns:
      tuple: (list of (rows, cols, dists) pair blocks, list of skipped
      (query positions, reference positions) blocks), all in the positions
      of the whole table.
    """
    found, skipped = [], []
    patterns = {}
    for ref in ref_rows[ref_x[ref_rows] > 0].tolist():
        seq = ref_seqs[ref]
        key = (len(seq), tuple(i for i, c in enumerate(seq) if c == WILDCARD))
        patterns.setdefault(key, []).append(ref)

    clean_refs = ref_rows[ref_x[ref_rows] == 0]
    if len(clean_refs) and len(query_rows):
        pairs, blocks = _symscanPairs(query_seqs, ref_seqs, limits + query_x,
                                      distance_type, query_rows, clean_refs)
        found.append(pairs)
        skipped.extend(blocks)
    for (length, wild), refs in patterns.items():
        members = query_rows[query_lengths[query_rows] == length]
        if not len(members):
            continue
        cut = set(wild)
        masked_queries = {q: _without(query_seqs[q], cut) for q in members.tolist()}
        masked_refs = {r: _without(ref_seqs[r], cut) for r in refs}
        pairs, blocks = _symscanPairs(masked_queries, masked_refs, limits + query_x,
                                      distance_type, members,
                                      numpy.array(refs, dtype=numpy.int64))
        found.append(pairs)
        skipped.extend(blocks)

    return found, skipped


def _symscanPairs(query_seqs, ref_seqs, limits, distance_type, query_rows,
                  ref_rows):
    """
    Every (query, reference) pair within its query's limit, via symscan,
    with each call kept within MAX_VARIANTS and MAX_DISTANCE.

    Queries are grouped by limit (and, under Hamming, length), and both
    sides of a group are chunked to half the budget each, so every call
    fits however large the tables. What cannot fit at all is skipped: a
    group whose limit exceeds MAX_DISTANCE, or a single sequence whose own
    index exceeds half the budget.

    Arguments:
      query_seqs (list): every query row's comparison string.
      ref_seqs (list): every reference row's comparison string.
      limits (numpy.ndarray): every query row's resolved distance limit.
      distance_type (str): 'hamming' or 'levenshtein'.
      query_rows (numpy.ndarray): the query positions to compare.
      ref_rows (numpy.ndarray): the reference positions to compare them to.

    Returns:
      tuple: ((query positions, reference positions, distances), skipped)
      -- numpy int arrays indexing the full query and reference tables, and
      a list of (query positions, reference positions) blocks that were not
      compared.
    """
    hamming = distance_type == 'hamming'
    half = MAX_VARIANTS // 2
    query_lengths = numpy.array([len(query_seqs[i]) for i in query_rows], dtype=int)
    ref_lengths = numpy.array([len(ref_seqs[i]) for i in ref_rows], dtype=int)

    # Under Hamming only same-length strings can match, so a group is one
    # (limit, length) and indexes only the references of that length.
    if hamming:
        keys = numpy.stack([limits[query_rows], query_lengths], axis=1)
    else:
        keys = limits[query_rows][:, None]
    groups = {}
    for position, key in enumerate(map(tuple, keys)):
        groups.setdefault(key, []).append(position)

    found, skipped = [], []
    for key, positions in groups.items():
        limit = int(key[0])
        members = query_rows[positions]
        refs = ref_rows[ref_lengths == key[1]] if hamming else ref_rows
        if len(refs) == 0:
            continue
        if limit > MAX_DISTANCE:
            skipped.append((members, refs))
            continue

        # Python ints: a long sequence at a high distance overflows int64.
        query_costs = [_deletionVariants(len(query_seqs[i]), limit) for i in members]
        ref_costs = [_deletionVariants(len(ref_seqs[i]), limit) for i in refs]
        query_fits = numpy.array([cost <= half for cost in query_costs], dtype=bool)
        ref_fits = numpy.array([cost <= half for cost in ref_costs], dtype=bool)
        fit_queries, fit_refs = members[query_fits], refs[ref_fits]
        if not ref_fits.all():
            skipped.append((members, refs[~ref_fits]))
        if not query_fits.all():
            skipped.append((members[~query_fits], fit_refs))

        query_chunks = _chunks(fit_queries, [c for c, ok in zip(query_costs, query_fits)
                                             if ok], half)
        ref_chunks = _chunks(fit_refs, [c for c, ok in zip(ref_costs, ref_fits) if ok],
                             half)
        for ref_chunk in ref_chunks:
            subset_refs = [ref_seqs[i] for i in ref_chunk]
            for chunk in query_chunks:
                r, c, d = symscan.get_neighbors_across(
                    [query_seqs[i] for i in chunk], subset_refs,
                    max_distance=limit, distance_type=distance_type)
                found.append((chunk[r.astype(numpy.int64)],
                              ref_chunk[c.astype(numpy.int64)], d.astype(numpy.int64)))

    return _concat(found), skipped


def runSymscanAnnotation(queries, reference, compare_cols, max_distance,
                         distance_type='hamming', id_col=DEFAULT_ID_COL,
                         exact_cols=None, threshold_is_fraction=False, threads=1,
                         gene_cols=None, return_pairs=False, warn_skipped=True):
    """
    Annotate each query row against `reference` allowing up to
    `max_distance` differences in `compare_cols`, via
    symscan.get_neighbors_across.

    symscan has no gate of its own, so `exact_cols` (and, for a pooled
    multi-column Hamming comparison, each column's length, so the columns
    line up one by one) and `gene_cols` are applied to the pairs symscan
    returns. `gene_cols` additionally partitions the search, which is only
    an optimization: the gate decides what is kept either way (see the
    module note). An 'X' matches
    any residue under 'hamming'; see the module note.

    Arguments:
      queries (pandas.DataFrame): the input rows to annotate.
      reference (pandas.DataFrame): the fixed table being annotated against.
      compare_cols (list): column names to compare -- pooled (concatenated)
        when more than one, which only 'hamming' supports.
      max_distance (float): the tolerance -- an absolute count, or a
        fraction of the query's summed comparison column length when
        `threshold_is_fraction` is True.
      distance_type (str): 'hamming' (substitutions only) or 'levenshtein'
        (substitutions, insertions and deletions combined).
      id_col (str): column in `reference` holding each row's identifier.
      exact_cols (list): column names that must agree exactly, e.g.
        `['v_call']`. None for no such gate.
      threshold_is_fraction (bool): treat `max_distance` as a fraction of
        the query's length instead of a fixed count.
      threads (int): rayon worker threads (-1 for every available core);
        only the first call in a process sets this -- see the module note.
      gene_cols (list): gene call columns that must name the same gene
        (any allele), e.g. `['v_call', 'j_call']`; a blank or unrecognized
        call agrees with nothing (see Genes.py). None for no such gate.
      return_pairs (bool): also return the pairs behind `hit_ids`.
      warn_skipped (bool): warn about the query rows beyond symscan's
        limits; False for a caller that scores them some other way (see
        Annotate._scoreSkippedRows).

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `min_mismatches`, `is_unique_hit`, `hit_ids`, plus
      `mismatch_positions` under 'hamming'. `min_mismatches` is the
      closest hit's distance (edit distance under 'levenshtein'), and
      `hit_ids` the hits' `id_col` values, closest first. Pairs beyond
      symscan's limits are skipped with a warning (see the module note);
      `attrs['skipped_query_rows']` counts the query rows that may
      therefore miss hits, and `attrs['skipped_rows']` lists them.
      With `return_pairs`, a tuple: that table, and one row per (query row,
      hit id) in the same order -- `query_row`, `reference_id`, `mismatches`,
      `identity` (1 - mismatches / the query's compared length, an X
      counting as identical under 'hamming'), plus `mismatch_positions`
      under 'hamming'.

    Raises:
      AnnotationError: if `distance_type` is not recognized, or
        'levenshtein' is given several columns.
    """
    if distance_type not in ('hamming', 'levenshtein'):
        raise AnnotationError("unknown distance_type %r; choose 'hamming' or "
                              "'levenshtein'" % (distance_type,))
    hamming = distance_type == 'hamming'
    if not hamming and len(compare_cols) != 1:
        raise AnnotationError("distance_type='levenshtein' compares a single column; "
                              'got %d' % len(compare_cols))

    query_values = _rowValues(queries, compare_cols)
    ref_values = _rowValues(reference, compare_cols)
    query_seqs = [''.join(v) for v in query_values]
    ref_seqs = [''.join(v) for v in ref_values]
    limits = _resolveLimits(query_values, max_distance, threshold_is_fraction)
    query_exact = _exactValues(queries, exact_cols)
    ref_exact = _exactValues(reference, exact_cols)

    _startPool(threads)
    # A blank value has nothing to compare: a row with one in any compared
    # column never matches, on either side.
    query_blank = numpy.array([any(v == '' for v in values) for values in query_values],
                              dtype=bool)
    ref_blank = numpy.array([any(v == '' for v in values) for values in ref_values],
                            dtype=bool)
    if query_blank.any() or ref_blank.any():
        # Only the side that has them, so a clean query says nothing.
        blank = ', '.join('%s %s rows' % (format(int(n.sum()), ','), label)
                          for label, n in (('query', query_blank),
                                           ('reference', ref_blank)) if n.any())
        log.info('  %s have a blank %s; they are not compared', blank,
                 '/'.join(compare_cols))
    all_queries = numpy.flatnonzero(~query_blank).astype(numpy.int64)
    all_refs = numpy.flatnonzero(~ref_blank).astype(numpy.int64)

    # The X wildcard (Hamming only); see the module note. References are
    # grouped by where their X's sit, and those positions are cut out of the
    # group's references and of every same-length query, so a reference X
    # never counts; a query's own X's are covered by slack instead.
    if hamming:
        query_x = numpy.array([seq.count(WILDCARD) for seq in query_seqs],
                              dtype=numpy.int64)
        ref_x = numpy.array([seq.count(WILDCARD) for seq in ref_seqs], dtype=numpy.int64)
    else:
        query_x = numpy.zeros(len(query_seqs), dtype=numpy.int64)
        ref_x = numpy.zeros(len(ref_seqs), dtype=numpy.int64)
    # One search per gene key rather than one over everything (see
    # _geneGroups): the gene gate below still decides what is kept, but
    # pairs of different genes are no longer compared only to be discarded.
    query_lengths = numpy.array([len(seq) for seq in query_seqs], dtype=int)
    found, skipped = [], []
    for group_queries, group_refs in _geneGroups(queries, reference, gene_cols,
                                                 all_queries, all_refs):
        pairs, blocks = _searchGroup(query_seqs, ref_seqs, query_lengths, limits,
                                     query_x, ref_x, distance_type,
                                     group_queries, group_refs)
        found.extend(pairs)
        skipped.extend(blocks)
    rows, cols, dists = _concat(found)
    skipped_queries = []
    if skipped:
        skipped_queries = numpy.unique(numpy.concatenate([q for q, _ in skipped]))
        skipped_refs = numpy.unique(numpy.concatenate([r for _, r in skipped]))
    if skipped and warn_skipped:
        log.warning(
            '%d query rows could not be compared with %d reference rows, beyond '
            "symscan's limits (over %d differences, counting one per %r, or a "
            'sequence too long for its index at this threshold); those pairs '
            'were skipped, so these query rows may miss hits: %s',
            len(skipped_queries), len(skipped_refs), MAX_DISTANCE, WILDCARD,
            _describeRows(skipped_queries, queries))

    # Gate symscan's pairs: exact columns always, per-column lengths when
    # several Hamming columns were pooled (the concatenated length alone
    # could match while the individual columns don't line up).
    query_keys, ref_keys = query_exact, ref_exact
    if hamming and len(compare_cols) > 1:
        query_keys = [k + tuple(map(len, v)) for k, v in zip(query_keys, query_values)]
        ref_keys = [k + tuple(map(len, v)) for k, v in zip(ref_keys, ref_values)]
    if exact_cols or (hamming and len(compare_cols) > 1):
        query_codes, ref_codes = _keyCodes(query_keys, ref_keys)
        keep = query_codes[rows] == ref_codes[cols]
        rows, cols, dists = rows[keep], cols[keep], dists[keep]

    # Re-score every pair involving an X under the wildcard rule, and drop
    # those the slack let in beyond the real limit.
    wild = numpy.flatnonzero((query_x[rows] > 0) | (ref_x[cols] > 0))
    if len(wild):
        dists = dists.copy()
        dists[wild] = [len(_mismatchPositions(query_values[rows[i]], ref_values[cols[i]]))
                       for i in wild]
        keep = dists <= limits[rows]
        rows, cols, dists = rows[keep], cols[keep], dists[keep]

    gate = buildGeneGate(queries, reference, gene_cols)
    if gate is not None:
        keep = gate.keep(rows, cols)
        rows, cols, dists = rows[keep], cols[keep], dists[keep]

    ref_ids = reference[id_col].to_numpy()
    pairs = pandas.DataFrame({'query_row': rows, 'col': cols, 'dist': dists,
                              'id': ref_ids[cols].astype(str)})
    # One reference id can appear on several rows; keep its closest one
    # (lowest reference position on a tie).
    pairs = (pairs.sort_values(['query_row', 'dist', 'id', 'col'])
             .drop_duplicates(['query_row', 'id']))

    grouped = pairs.groupby('query_row', sort=True)
    hits = pandas.DataFrame({'n_hits_total': grouped.size(),
                             'min_mismatches': grouped['dist'].min().astype(float),
                             'hit_ids': grouped['id'].agg(';'.join)})
    columns = ['query_row', 'n_hits_total', 'min_mismatches', 'is_unique_hit',
               'hit_ids']
    if hamming:
        n_cols = len(compare_cols)
        pairs['positions'] = [
            _formatPositions(_mismatchPositions(query_values[q], ref_values[c]), n_cols)
            for q, c in zip(pairs['query_row'].tolist(), pairs['col'].tolist())]
        hits['mismatch_positions'] = pairs.groupby('query_row', sort=True)[
            'positions'].agg(';'.join)
        columns.append('mismatch_positions')

    result = pandas.DataFrame({'query_row': numpy.arange(len(queries))}).join(
        hits, on='query_row')
    result['n_hits_total'] = result['n_hits_total'].fillna(0).astype(int)
    # An all-missing column comes out float: cast it to object before filling
    # with '' (pandas will refuse a str fill on a float column), then to str.
    text_cols = ['hit_ids'] + (['mismatch_positions'] if hamming else [])
    for col in text_cols:
        result[col] = result[col].astype(object).fillna('').astype(str)
    result['is_unique_hit'] = result['n_hits_total'] <= 1
    result = result[columns]
    result.attrs['skipped_query_rows'] = len(skipped_queries)
    result.attrs['skipped_rows'] = [int(row) for row in skipped_queries]
    if not return_pairs:
        return result

    # The pairs behind hit_ids, in the same order. Identity is against the
    # query's compared length, the length identity_threshold scales with.
    query_rows = pairs['query_row'].to_numpy()
    lengths = numpy.array([len(seq) for seq in query_seqs], dtype=numpy.int64)
    pair_table = pandas.DataFrame({
        'query_row': query_rows, 'reference_id': pairs['id'].to_numpy(),
        'mismatches': pairs['dist'].to_numpy(),
        'identity': 1 - pairs['dist'].to_numpy() / lengths[query_rows]})
    if hamming:
        pair_table['mismatch_positions'] = pairs['positions'].to_numpy()

    return result, pair_table
