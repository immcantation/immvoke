"""
Row hashes, for specificity sources

Every specificity source tags each normalized record with
immvoke_row_hash: a short content hash of its upstream columns, so a row can
be recognized across downloads. Kept here rather than in each source, since
they all hash the same way.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import hashlib

import pandas


def rowHashes(frame):
    """
    Build a short content hash for each row of a frame.

    A row's key is its values as text joined by '|', the way
    '|'.join(str(x) for x in row) spells it, so hashes do not change with
    how they are computed. The values are turned to text a column at a time
    rather than with apply(..., axis=1), which builds a Series per row and
    took 45 of the 46 seconds of normalizing IEDB's 338,000 row tcr table.

    Arguments:
      frame (pandas.DataFrame): upstream columns, before provenance columns
        are added.

    Returns:
      pandas.Series: the first 12 hex characters of each row's SHA-256
      digest, on the frame's index.
    """
    if not len(frame) or not len(frame.columns):
        return pandas.Series([], dtype=str, index=frame.index)

    # str() on each value, as the old row-wise key did; pandas' own
    # astype(str) would leave missing values as NaN rather than 'nan'.
    columns = [map(str, frame.iloc[:, i].tolist()) for i in range(len(frame.columns))]

    return pandas.Series([hashlib.sha256('|'.join(row).encode('utf-8')).hexdigest()[:12]
                          for row in zip(*columns)], index=frame.index, dtype=str)
