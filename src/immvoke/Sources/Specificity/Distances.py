"""
Every gated pair's distance, whatever the threshold

symscan only ever returns the pairs within a threshold. The pairs file of a
Hamming or Levenshtein run lists every (query, reference) pair that passes
the run's gates instead -- the gene gate (`gene_cols`), the exact columns
(`exact_cols`) and, under Hamming, each compared column's length -- with
its distance, within the threshold or not. allPairDistances scores those
pairs directly with numpy:

    Hamming      a block of queries against the same-length references,
                 compared position by position; an X on either side never
                 counts (as in Symscan.py)
    Levenshtein  the edit distance table, filled one query residue at a
                 time for a whole block of pairs at once; an X is an
                 ordinary letter (as in Symscan.py)

Candidates come from indexes rather than from testing every pair:
references are indexed by their exact-column key (with the per-column
lengths under Hamming) and by every gene they name, so a query only meets
the references that pass its gates. Queries go in blocks of consecutive
rows, so the pairs come out in `query_row` order and never all sit in
memory at once -- without a gene gate, a repertoire against IEDB is on the
order of 10^8 pairs.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging

import numpy
import pandas

# Immvoke imports
from immvoke.Sources.Specificity.Annotate import (
    DEFAULT_ID_COL,
    WILDCARD,
    _exactValues,
    _formatPositions,
    _rowValues,
)
from immvoke.Sources.Specificity.Genes import parseCall
from immvoke.Sources.Specificity.Symscan import _resolveLimits

log = logging.getLogger(__name__)

#: Consecutive query rows scored, and yielded, together.
QUERY_BLOCK = 2000

#: Most residues (pairs x sequence length) one vectorized step holds, which
#: bounds its memory to some tens of MB.
RESIDUE_BUDGET = 5_000_000

_WILDCARD_CODE = ord(WILDCARD)
_NO_ROWS = numpy.array([], dtype=numpy.int64)


def _codes(seqs):
    """
    Arguments:
      seqs (list): equal-length strings.

    Returns:
      numpy.ndarray: one row of byte codes per string.
    """
    joined = ''.join(seqs).encode('latin-1', 'replace')

    return numpy.frombuffer(joined, dtype=numpy.uint8).reshape(len(seqs), -1)


def _byLength(rows, seqs):
    """
    Returns:
      dict: sequence length to the `rows` (a numpy array) with that length.
    """
    lengths = numpy.array([len(seqs[i]) for i in rows], dtype=numpy.int64)

    return {int(length): rows[lengths == length] for length in numpy.unique(lengths)}


def _levenshtein(a, b):
    """
    The edit distance of every pair at once.

    Arguments:
      a (numpy.ndarray): query codes, one row per pair.
      b (numpy.ndarray): reference codes, one row per pair.

    Returns:
      numpy.ndarray: one distance per pair.
    """
    pairs, n = b.shape
    steps = numpy.arange(1, n + 1, dtype=numpy.int32)
    prev = numpy.tile(numpy.arange(n + 1, dtype=numpy.int32), (pairs, 1))
    for i in range(1, a.shape[1] + 1):
        # Substitution (or match) and deletion first; insertion runs along
        # the row: cur[j] = min(best[j], cur[j-1] + 1) with cur[0] = i,
        # which unrolls to min(i + j, j + min over k <= j of best[k] - k).
        best = numpy.minimum(prev[:, :-1] + (b != a[:, i - 1:i]), prev[:, 1:] + 1)
        cur = numpy.empty_like(prev)
        cur[:, 0] = i
        cur[:, 1:] = numpy.minimum(
            i + steps, numpy.minimum.accumulate(best - steps, axis=1) + steps)
        prev = cur

    return prev[:, n]


def _levenshteinPairs(query_seqs, ref_seqs, query_rows, ref_rows):
    """
    Score every query of `query_rows` against every reference of `ref_rows`
    under Levenshtein, a block of same-length queries against same-length
    references at a time.

    Yields:
      tuple: (query positions, reference positions, distances), numpy
      arrays.
    """
    query_groups = _byLength(query_rows, query_seqs)
    for ref_length, refs in _byLength(ref_rows, ref_seqs).items():
        ref_codes = _codes([ref_seqs[i] for i in refs])
        for query_length, queries in query_groups.items():
            query_codes = _codes([query_seqs[i] for i in queries])
            step = max(1, RESIDUE_BUDGET // (len(refs) * (query_length + ref_length)))
            for start in range(0, len(queries), step):
                block = query_codes[start:start + step]
                distances = _levenshtein(numpy.repeat(block, len(refs), axis=0),
                                         numpy.tile(ref_codes, (len(block), 1)))
                yield (numpy.repeat(queries[start:start + step], len(refs)),
                       numpy.tile(refs, len(block)), distances)


def _hammingPairs(query_seqs, ref_seqs, query_rows, ref_rows, col_lengths):
    """
    Score every query of `query_rows` against every reference of `ref_rows`
    under Hamming -- all one length, the gate saw to that.

    Arguments:
      col_lengths (tuple): the length of each compared column, pooled into
        the sequences in this order, for the mismatch positions.

    Yields:
      tuple: (query positions, reference positions, distances, mismatch
      positions as strings, see Annotate._formatPositions).
    """
    starts = numpy.cumsum((0,) + tuple(col_lengths))[:-1]
    ref_codes = _codes([ref_seqs[i] for i in ref_rows])
    query_codes = _codes([query_seqs[i] for i in query_rows])
    step = max(1, RESIDUE_BUDGET // (len(ref_rows) * ref_codes.shape[1]))
    for start in range(0, len(query_rows), step):
        a = query_codes[start:start + step, None, :]
        b = ref_codes[None, :, :]
        differ = (a != b) & (a != _WILDCARD_CODE) & (b != _WILDCARD_CODE)
        differ = differ.reshape(-1, ref_codes.shape[1])
        positions = []
        for row in differ:
            where = numpy.flatnonzero(row)
            cols = numpy.searchsorted(starts, where, side='right') - 1
            positions.append(_formatPositions(
                list(zip(cols.tolist(), (where - starts[cols]).tolist())),
                len(col_lengths)))
        count = len(a)
        yield (numpy.repeat(query_rows[start:start + step], len(ref_rows)),
               numpy.tile(ref_rows, count), differ.sum(axis=1), positions)


def allPairDistances(queries, reference, compare_cols, max_distance,
                     distance_type='hamming', id_col=DEFAULT_ID_COL, exact_cols=None,
                     threshold_is_fraction=False, gene_cols=None):
    """
    Score every (query, reference) pair that passes the gates, whatever its
    distance.

    Arguments: as Symscan.runSymscanAnnotation -- `max_distance` and
      `threshold_is_fraction` only decide `within_threshold` here.

    Yields:
      pandas.DataFrame: the pairs of one block of consecutive query rows,
      in `query_row` order, closest first: `query_row`, `reference_id` (the
      reference's `id_col`, once per query row -- its closest row),
      `mismatches` (the distance), `identity` (1 - mismatches / the query's
      compared length), `within_threshold` (whether symscan's threshold
      admits it), plus `mismatch_positions` under 'hamming'. At least one
      table, empty when no pair passes.
    """
    hamming = distance_type == 'hamming'
    query_values = _rowValues(queries, compare_cols)
    ref_values = _rowValues(reference, compare_cols)
    query_seqs = [''.join(v) for v in query_values]
    ref_seqs = [''.join(v) for v in ref_values]
    limits = _resolveLimits(query_values, max_distance, threshold_is_fraction)

    # Gate keys, as Symscan applies them: exact columns, and under Hamming
    # each compared column's length (only same-length strings compare).
    query_keys = _exactValues(queries, exact_cols)
    ref_keys = _exactValues(reference, exact_cols)
    if hamming:
        query_keys = [k + tuple(map(len, v)) for k, v in zip(query_keys, query_values)]
        ref_keys = [k + tuple(map(len, v)) for k, v in zip(ref_keys, ref_values)]
    query_blank = [any(v == '' for v in values) for values in query_values]
    by_key = {}
    for row, (key, values) in enumerate(zip(ref_keys, ref_values)):
        if not any(v == '' for v in values):
            by_key.setdefault(key, []).append(row)
    by_key = {key: numpy.array(rows, dtype=numpy.int64) for key, rows in by_key.items()}

    # Gene gate: each reference row under every gene it names, per column.
    def parsed(frame, col):
        values = frame[col].astype(object).fillna('').astype(str).tolist()
        calls = {value: parseCall(value) for value in set(values)}
        return [calls[value] for value in values]

    gene_cols = gene_cols or []
    query_genes = list(zip(*(parsed(queries, col) for col in gene_cols))) or \
        [()] * len(queries)
    gene_index = []
    for col in gene_cols:
        index = {}
        for row, call in enumerate(parsed(reference, col)):
            for gene in call:
                index.setdefault(gene, []).append(row)
        gene_index.append({gene: numpy.array(rows, dtype=numpy.int64)
                           for gene, rows in index.items()})

    cache = {}

    def candidates(key, genes):
        if (key, genes) not in cache:
            rows = by_key.get(key, _NO_ROWS)
            for index, call in zip(gene_index, genes):
                # A blank or unrecognized call agrees with nothing.
                agreeing = [index[g] for g in call if g in index]
                rows = (numpy.intersect1d(rows, numpy.concatenate(agreeing))
                        if agreeing else _NO_ROWS)
            cache[key, genes] = rows
        return cache[key, genes]

    ref_ids = reference[id_col].astype(str).to_numpy()
    lengths = numpy.array([len(seq) for seq in query_seqs], dtype=numpy.int64)
    columns = ['query_row', 'reference_id', 'mismatches', 'identity',
               'within_threshold'] + (['mismatch_positions'] if hamming else [])
    yielded = False
    for start in range(0, len(queries), QUERY_BLOCK):
        groups = {}
        for row in range(start, min(start + QUERY_BLOCK, len(queries))):
            if not query_blank[row]:
                groups.setdefault((query_keys[row], query_genes[row]), []).append(row)

        parts = []
        for (key, genes), rows in groups.items():
            refs = candidates(key, genes)
            if not len(refs):
                continue
            rows = numpy.array(rows, dtype=numpy.int64)
            if hamming:
                col_lengths = key[len(key) - len(compare_cols):]
                parts.extend(_hammingPairs(query_seqs, ref_seqs, rows, refs, col_lengths))
            else:
                parts.extend(part + (None,) for part in
                             _levenshteinPairs(query_seqs, ref_seqs, rows, refs))

        if parts:
            pairs = pandas.DataFrame({
                'query_row': numpy.concatenate([p[0] for p in parts]),
                'col': numpy.concatenate([p[1] for p in parts]),
                'dist': numpy.concatenate([p[2] for p in parts])})
            pairs['id'] = ref_ids[pairs['col'].to_numpy()]
            if hamming:
                pairs['positions'] = [x for p in parts for x in p[3]]
            # One reference id can sit on several rows; keep its closest one
            # (lowest reference position on a tie), as Symscan does.
            pairs = (pairs.sort_values(['query_row', 'dist', 'id', 'col'])
                     .drop_duplicates(['query_row', 'id']))
        else:
            pairs = pandas.DataFrame({'query_row': _NO_ROWS, 'col': _NO_ROWS,
                                      'dist': _NO_ROWS, 'id': [], 'positions': []})
        rows, dist = pairs['query_row'].to_numpy(), pairs['dist'].to_numpy()
        table = pandas.DataFrame({
            'query_row': rows, 'reference_id': pairs['id'].to_numpy(),
            'mismatches': dist, 'identity': 1 - dist / lengths[rows],
            'within_threshold': dist <= limits[rows]})
        if hamming:
            table['mismatch_positions'] = pairs['positions'].to_numpy()
        if len(table) or not yielded:
            yielded = True
            yield table[columns]

    if not yielded:
        yield pandas.DataFrame({col: [] for col in columns})
