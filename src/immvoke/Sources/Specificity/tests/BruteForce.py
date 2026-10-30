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
    result, pairs, _ = _bruteForce(queries, reference, compare_cols, max_distance,
                                   distance_type, id_col, exact_cols,
                                   threshold_is_fraction, gene_cols)

    return (result, pairs) if return_pairs else result


def bruteForceAllPairs(queries, reference, compare_cols, max_distance,
                       distance_type='hamming', id_col='id', exact_cols=None,
                       threshold_is_fraction=False, gene_cols=None):
    """
    The pairs Distances.allPairDistances must yield, concatenated, by brute
    force: every pair through the gates, whatever its distance.

    Arguments: as Symscan.runSymscanAnnotation.

    Returns:
      pandas.DataFrame: as Distances.allPairDistances's tables.
    """
    return _bruteForce(queries, reference, compare_cols, max_distance, distance_type,
                       id_col, exact_cols, threshold_is_fraction, gene_cols)[2]


def _bruteForce(queries, reference, compare_cols, max_distance, distance_type,
                id_col, exact_cols, threshold_is_fraction, gene_cols):
    """
    Returns:
      tuple: (the result table, its pairs, every gated pair).
    """
    hamming = distance_type == 'hamming'
    query_values = list(queries[compare_cols].fillna('').itertuples(index=False,
                                                                   name=None))
    ref_values = list(reference[compare_cols].fillna('').itertuples(index=False,
                                                                   name=None))
    exact = exact_cols or []
    gate = GeneGate(queries, reference, gene_cols) if gene_cols else None

    rows, pairs, all_pairs = [], [], []
    for q, values in enumerate(query_values):
        blank = '' in values
        length = sum(len(v) for v in values)
        limit = (_allowedDifferences(max_distance, length) if threshold_is_fraction
                 else int(max_distance))
        best = {}  # id -> (distance, reference position, mismatch positions)
        scored = {}  # the same, whatever the distance
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
            hit = str(reference[id_col].iloc[r])
            if hit not in scored or (distance, r) < scored[hit][:2]:
                scored[hit] = (distance, r, positions)
            if distance > limit:
                continue
            if hit not in best or (distance, r) < best[hit][:2]:
                best[hit] = (distance, r, positions)

        single = len(compare_cols) == 1

        def formatted(found, x):
            return ','.join(str(pos) if single else '%d:%d' % (col, pos)
                            for col, pos in found[x][2])

        order = sorted(best, key=lambda x: (best[x][0], x))
        row = {'query_row': q, 'n_hits_total': len(order),
               'min_mismatches': float(best[order[0]][0]) if order else float('nan'),
               'is_unique_hit': len(order) <= 1, 'hit_ids': ';'.join(order)}
        if hamming:
            row['mismatch_positions'] = ';'.join(formatted(best, x) for x in order)
        rows.append(row)
        for x in order:
            pair = {'query_row': q, 'reference_id': x, 'mismatches': best[x][0],
                    'identity': 1 - best[x][0] / length}
            if hamming:
                pair['mismatch_positions'] = formatted(best, x)
            pairs.append(pair)
        for x in sorted(scored, key=lambda x: (scored[x][0], x)):
            pair = {'query_row': q, 'reference_id': x, 'mismatches': scored[x][0],
                    'identity': 1 - scored[x][0] / length,
                    'within_threshold': scored[x][0] <= limit}
            if hamming:
                pair['mismatch_positions'] = formatted(scored, x)
            all_pairs.append(pair)

    extra = ['mismatch_positions'] if hamming else []
    columns = ['query_row', 'reference_id', 'mismatches', 'identity']

    return (pandas.DataFrame(rows), pandas.DataFrame(pairs, columns=columns + extra),
            pandas.DataFrame(all_pairs, columns=columns + ['within_threshold'] + extra))
