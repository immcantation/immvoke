"""
A brute-force reference matcher, for tests only

Compares every query row with every reference row directly -- no index, no
grouping, no slack -- and builds the same result and pairs tables as
Symscan.runSymscanAnnotation, so symscan's answers can be checked against an
implementation simple enough to trust by reading it. Far too slow for real
tables.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import pandas

# Immvoke imports
from immvoke.Sources.Specificity.Annotate import WILDCARD, _allowedDifferences
from immvoke.Sources.Specificity.Genes import GeneGate


def _hamming(a, b):
    """
    Returns:
      list: (column, position) of every difference between two rows' column
      tuples, an X on either side never counting -- or None if any column's
      lengths differ.
    """
    if any(len(x) != len(y) for x, y in zip(a, b)):
        return None

    return [(col, pos) for col, (x, y) in enumerate(zip(a, b))
            for pos, (cx, cy) in enumerate(zip(x, y))
            if cx != cy and WILDCARD not in (cx, cy)]


def _levenshtein(a, b):
    """
    Returns:
      int: the plain edit distance between two strings.
    """
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1,
                               previous[j - 1] + (ca != cb)))
        previous = current

    return previous[-1]


def bruteForceAnnotation(queries, reference, compare_cols, max_distance,
                         distance_type='hamming', id_col='id', exact_cols=None,
                         threshold_is_fraction=False, gene_cols=None,
                         return_pairs=False):
    """
    The result table runSymscanAnnotation must produce, by brute force.

    Arguments: as Symscan.runSymscanAnnotation.

    Returns:
      pandas.DataFrame: as Symscan.runSymscanAnnotation.
    """
    hamming = distance_type == 'hamming'
    query_values = list(queries[compare_cols].fillna('').itertuples(index=False,
                                                                   name=None))
    ref_values = list(reference[compare_cols].fillna('').itertuples(index=False,
                                                                   name=None))
    exact = exact_cols or []
    gate = GeneGate(queries, reference, gene_cols) if gene_cols else None

    rows, pairs = [], []
    for q, values in enumerate(query_values):
        blank = '' in values
        length = sum(len(v) for v in values)
        limit = (_allowedDifferences(max_distance, length) if threshold_is_fraction
                 else int(max_distance))
        best = {}  # id -> (distance, reference position, mismatch positions)
        for r, candidate in enumerate(ref_values):
            if blank or '' in candidate:
                continue
            if any(queries[c].iloc[q] != reference[c].iloc[r] for c in exact):
                continue
            if gate is not None and not gate.agrees(q, r):
                continue
            if hamming:
                positions = _hamming(values, candidate)
                if positions is None:
                    continue
                distance = len(positions)
            else:
                positions, distance = None, _levenshtein(values[0], candidate[0])
            if distance > limit:
                continue
            hit = str(reference[id_col].iloc[r])
            if hit not in best or (distance, r) < best[hit][:2]:
                best[hit] = (distance, r, positions)

        order = sorted(best, key=lambda x: (best[x][0], x))
        row = {'query_row': q, 'n_hits_total': len(order),
               'min_mismatches': float(best[order[0]][0]) if order else float('nan'),
               'is_unique_hit': len(order) <= 1, 'hit_ids': ';'.join(order)}
        positions = {}
        if hamming:
            single = len(compare_cols) == 1
            positions = {x: ','.join(str(pos) if single else '%d:%d' % (col, pos)
                                     for col, pos in best[x][2])
                         for x in order}
            row['mismatch_positions'] = ';'.join(positions[x] for x in order)
        rows.append(row)
        for x in order:
            pair = {'query_row': q, 'hit_id': x, 'mismatches': best[x][0],
                    'identity': 1 - best[x][0] / length}
            if hamming:
                pair['mismatch_positions'] = positions[x]
            pairs.append(pair)

    result = pandas.DataFrame(rows)
    if not return_pairs:
        return result

    columns = ['query_row', 'hit_id', 'mismatches', 'identity']
    return result, pandas.DataFrame(
        pairs, columns=columns + (['mismatch_positions'] if hamming else []))
