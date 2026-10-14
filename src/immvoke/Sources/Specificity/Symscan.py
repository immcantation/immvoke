"""
symscan-accelerated annotation: exact, Hamming and Levenshtein

    runSymscanAnnotation       up to `max_distance` substitutions
                               (distance_type 'hamming', one or more pooled
                               columns) or substitutions/insertions/deletions
                               combined ('levenshtein', a single column); an
                               exact match is `max_distance` 0

symscan (https://github.com/yutanagano/symscan) builds a Rust
symmetric-deletion index and only ever looks at likely matches, rather than
scoring every query x reference pair -- so it stays fast as the reference
grows. It returns the same result table as its pure-Python counterpart in
Annotate.py, so the two are interchangeable.

symscan has no gate of its own: `exact_cols` (exact string agreement) and
`gene_cols` (V/J-gene agreement at `gene_match` level, see Genes.py) are
applied to the pairs it returns -- distance first, gene check second.

symscan cannot compare every row, so the rows it cannot are compared in
Python instead, under the same rules, and merged with its own pairs:

  * It has no concept of the "X" wildcard. Under 'hamming' (and so under
    exact matching too), query rows that carry one are compared against the
    whole reference, and clean query rows against the reference's own
    X-containing rows. Levenshtein has no wildcard handling under either
    backend, so it never splits.
  * Its index holds every deletion variant of every string on both sides --
    sum(C(length, k) for k <= distance) per string -- and it aborts the whole
    process (an uncatchable Rust panic or failed allocation) once that
    outgrows memory, or once the distance exceeds MAX_DISTANCE. Each call is
    therefore kept within MAX_VARIANTS by grouping queries by their resolved
    distance (and, under 'hamming', by length, since only same-length
    references can match) and chunking the queries in a group; a group whose
    reference side alone exceeds the budget, or whose distance exceeds
    MAX_DISTANCE, is compared in Python.

symscan parallelizes via rayon's global pool, which reads RAYON_NUM_THREADS
once, on first use -- so `threads` only takes effect on the first call in a
process; later calls reuse that pool.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
import math
import os
from functools import cache

import numpy
import pandas
import symscan
from rapidfuzz.distance import Levenshtein

# Immvoke imports
from immvoke.Exceptions import AnnotationError
from immvoke.Sources.Specificity.Annotate import (
    DEFAULT_ID_COL,
    WILDCARD,
    _allowedDifferences,
    _bucketKey,
    _exactValues,
    _formatPositions,
    _mismatchPositions,
    _rowDistance,
    _rowValues,
)
from immvoke.Sources.Specificity.Genes import buildGeneGate

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


def _hasWildcard(values):
    """
    Arguments:
      values (list): per-row comparison column tuples.

    Returns:
      numpy.ndarray: one bool per row, True if any of its values holds
      WILDCARD.
    """
    return numpy.array([any(WILDCARD in v for v in row) for row in values],
                       dtype=bool)


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


def _symscanPairs(query_seqs, ref_seqs, limits, distance_type, query_rows,
                  ref_rows):
    """
    Every (query, reference) pair within its query's limit, via symscan,
    with each call kept within MAX_VARIANTS and MAX_DISTANCE.

    Arguments:
      query_seqs (list): every query row's comparison string.
      ref_seqs (list): every reference row's comparison string.
      limits (numpy.ndarray): every query row's resolved distance limit.
      distance_type (str): 'hamming' or 'levenshtein'.
      query_rows (numpy.ndarray): the query positions to compare.
      ref_rows (numpy.ndarray): the reference positions to compare them to.

    Returns:
      tuple: ((query positions, reference positions, distances), leftover)
      -- numpy int arrays indexing the full query and reference tables, and
      a list of (query positions, reference positions) groups symscan could
      not index, for the caller to compare another way.
    """
    hamming = distance_type == 'hamming'
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

    found, leftover = [], []
    for key, positions in groups.items():
        limit = int(key[0])
        members = query_rows[positions]
        refs = ref_rows[ref_lengths == key[1]] if hamming else ref_rows
        if len(refs) == 0:
            continue

        ref_cost = sum(_deletionVariants(len(ref_seqs[i]), limit) for i in refs)
        if limit > MAX_DISTANCE or ref_cost > MAX_VARIANTS:
            leftover.append((members, refs))
            continue

        # Chunk the queries so each call, reference side included, fits.
        subset_refs = [ref_seqs[i] for i in refs]
        chunk, chunk_cost, chunks, too_large = [], ref_cost, [], []
        for row in members:
            cost = _deletionVariants(len(query_seqs[row]), limit)
            if ref_cost + cost > MAX_VARIANTS:
                too_large.append(row)
                continue
            if chunk and chunk_cost + cost > MAX_VARIANTS:
                chunks.append(chunk)
                chunk, chunk_cost = [], ref_cost
            chunk.append(row)
            chunk_cost += cost
        if chunk:
            chunks.append(chunk)
        if too_large:
            leftover.append((numpy.array(too_large, dtype=numpy.int64), refs))

        for chunk in chunks:
            chunk = numpy.array(chunk, dtype=numpy.int64)
            r, c, d = symscan.get_neighbors_across(
                [query_seqs[i] for i in chunk], subset_refs,
                max_distance=limit, distance_type=distance_type)
            found.append((chunk[r.astype(numpy.int64)], refs[c.astype(numpy.int64)],
                          d.astype(numpy.int64)))

    return _concat(found), leftover


def _hammingPairs(query_values, query_exact, ref_values, ref_exact, limits,
                  query_rows, ref_rows):
    """
    Every (query, reference) Hamming pair within its query's limit, compared
    in Python under Annotate.py's wildcard rule -- for the rows symscan
    cannot compare.

    Arguments:
      query_values (list): every query row's comparison column tuple.
      query_exact (list): every query row's exact column tuple.
      ref_values (list): every reference row's comparison column tuple.
      ref_exact (list): every reference row's exact column tuple.
      limits (numpy.ndarray): every query row's resolved distance limit.
      query_rows (numpy.ndarray): the query positions to compare.
      ref_rows (numpy.ndarray): the reference positions to compare them to.

    Returns:
      tuple: (query positions, reference positions, distances), numpy int
      arrays indexing the full query and reference tables. Already gated on
      exact columns and per-column length, since the bucket key is both.
    """
    index = {}
    for col in ref_rows:
        index.setdefault(_bucketKey(ref_values[col], ref_exact[col]), []).append(col)

    rows, cols, dists = [], [], []
    for row in query_rows:
        for col in index.get(_bucketKey(query_values[row], query_exact[row]), []):
            distance = _rowDistance(query_values[row], ref_values[col])
            if distance is not None and distance <= limits[row]:
                rows.append(row)
                cols.append(col)
                dists.append(distance)

    return _concat([(rows, cols, dists)])


def _levenshteinPairs(query_seqs, query_exact, ref_seqs, ref_exact, limits,
                      query_rows, ref_rows):
    """
    Every (query, reference) Levenshtein pair within its query's limit,
    compared in Python (rapidfuzz) -- for the groups symscan cannot index.

    Arguments:
      query_seqs (list): every query row's comparison string.
      query_exact (list): every query row's exact column tuple.
      ref_seqs (list): every reference row's comparison string.
      ref_exact (list): every reference row's exact column tuple.
      limits (numpy.ndarray): every query row's resolved distance limit.
      query_rows (numpy.ndarray): the query positions to compare.
      ref_rows (numpy.ndarray): the reference positions to compare them to.

    Returns:
      tuple: (query positions, reference positions, distances), numpy int
      arrays indexing the full query and reference tables, already gated on
      exact columns.
    """
    index = {}
    for col in ref_rows:
        index.setdefault(ref_exact[col], []).append(col)

    rows, cols, dists = [], [], []
    for row in query_rows:
        limit = int(limits[row])
        for col in index.get(query_exact[row], []):
            distance = Levenshtein.distance(query_seqs[row], ref_seqs[col],
                                            score_cutoff=limit)
            if distance <= limit:
                rows.append(row)
                cols.append(col)
                dists.append(distance)

    return _concat([(rows, cols, dists)])


def runSymscanAnnotation(queries, reference, compare_cols, max_distance,
                         distance_type='hamming', id_col=DEFAULT_ID_COL,
                         exact_cols=None, threshold_is_fraction=False, threads=1,
                         gene_cols=None, gene_match='gene'):
    """
    Annotate each query row against `reference` allowing up to
    `max_distance` differences in `compare_cols`, via
    symscan.get_neighbors_across.

    symscan has no gate of its own, so `exact_cols` (and, for a pooled
    multi-column Hamming comparison, each column's length -- the same
    per-column length rule Annotate.runFuzzyAnnotation buckets on) and
    `gene_cols` are applied to the pairs symscan returns rather than before
    it runs: the distance comes first, the gene check second. Rows
    symscan cannot compare are compared in Python; see the module note.

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
      gene_cols (list): gene call columns that must agree at `gene_match`
        level, e.g. `['v_call', 'j_call']`; blank or unrecognized calls
        agree with anything (see Genes.py). None for no such gate.
      gene_match (str): 'allele', 'gene' or 'family'; see Genes.py.

    Returns:
      pandas.DataFrame: one row per query row, columns `query_row`,
      `n_hits_total`, `min_mismatches`, `is_unique_hit`, `hit_ids`, plus
      `mismatch_positions` under 'hamming' -- identical to
      Annotate.runFuzzyAnnotation ('hamming') or runLevenshteinAnnotation
      ('levenshtein').

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

    # Only Hamming has a wildcard rule to honour; see the module note.
    if hamming:
        query_x = _hasWildcard(query_values)
        ref_x = _hasWildcard(ref_values)
    else:
        query_x = numpy.zeros(len(query_values), dtype=bool)
        ref_x = numpy.zeros(len(ref_values), dtype=bool)
    clean_queries, x_queries = numpy.flatnonzero(~query_x), numpy.flatnonzero(query_x)
    clean_refs, x_refs = numpy.flatnonzero(~ref_x), numpy.flatnonzero(ref_x)
    if len(x_queries) or len(x_refs):
        log.info("'X' wildcard in %d query and %d reference rows; comparing "
                 'those in Python, the rest via symscan',
                 len(x_queries), len(x_refs))

    _startPool(threads)
    (rows, cols, dists), leftover = _symscanPairs(query_seqs, ref_seqs, limits,
                                                  distance_type, clean_queries,
                                                  clean_refs)

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

    # Everything symscan did not compare, in Python (already gated): the
    # groups beyond its limits, then X-containing query rows vs the whole
    # reference and clean query rows vs the reference's X-containing rows.
    found = [(rows, cols, dists)]
    if leftover:
        log.warning('%d query rows need more than symscan can index (a '
                    'distance above %d, or too many deletion variants for '
                    'their length); comparing those in Python',
                    sum(len(members) for members, _ in leftover), MAX_DISTANCE)
    for members, refs in leftover:
        if hamming:
            found.append(_hammingPairs(query_values, query_exact, ref_values,
                                       ref_exact, limits, members, refs))
        else:
            found.append(_levenshteinPairs(query_seqs, query_exact, ref_seqs,
                                           ref_exact, limits, members, refs))
    if len(x_queries) or len(x_refs):
        found.append(_hammingPairs(query_values, query_exact, ref_values, ref_exact,
                                   limits, x_queries, numpy.arange(len(ref_values))))
        found.append(_hammingPairs(query_values, query_exact, ref_values, ref_exact,
                                   limits, clean_queries, x_refs))
    rows, cols, dists = _concat(found)
    gate = buildGeneGate(queries, reference, gene_cols, gene_match)
    if gate is not None:
        keep = gate.keep(rows, cols)
        rows, cols, dists = rows[keep], cols[keep], dists[keep]

    ref_ids = reference[id_col].to_numpy()
    pairs = pandas.DataFrame({'query_row': rows, 'col': cols, 'dist': dists,
                              'id': ref_ids[cols].astype(str)})
    # One reference id can appear on several rows; keep its closest one
    # (lowest reference position on a tie, as runFuzzyAnnotation does).
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
    # astype(str): an all-missing column would otherwise stay object dtype.
    result['hit_ids'] = result['hit_ids'].fillna('').astype(str)
    if hamming:
        result['mismatch_positions'] = result['mismatch_positions'].fillna('').astype(str)
    result['is_unique_hit'] = result['n_hits_total'] <= 1

    return result[columns]
