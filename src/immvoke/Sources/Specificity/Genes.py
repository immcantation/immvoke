"""
V/J-gene call parsing and comparison for the annotation gene gate

An AIRR `v_call`/`j_call` cell is not always one clean name: IEDB alone mixes
allele level calls (`IGHV1-2*02`), gene level ones (`IGHV1-2`), several calls
in one cell (`IGHV1-2*02, IGHV1-2*04`) and blanks. Comparing those as plain
strings makes `IGHV1-2*02` and `IGHV1-2` disagree, so the gene gate parses
each cell into the genes it names and compares those instead: two calls
agree on the same gene, whatever the allele (IGHV1-2*02 == IGHV1-2*04 ==
IGHV1-2, but IGHV1-2 != IGHV1-69).

Two cells agree when they name a gene in common. A blank cell, or one naming
no recognizable gene, has nothing to compare and agrees with nothing -- not
even another blank.

The gene pattern is Change-O's (changeo/Gene.py, Immcantation), so a name
parses here the way Change-O's getGene would parse it. The allele is read
as an optional `*` suffix of the gene instead of with Change-O's allele
pattern, which also accepts `-` as the allele separator and so would read a
gene level `IGHV1-2` as gene `IGHV1` allele `2`.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
import re

import numpy

log = logging.getLogger(__name__)

#: Change-O's gene pattern, plus an optional `*allele` suffix.
_CALL_REGEX = re.compile(r'(?P<gene>(IG[HLK]|TR[ABGD])([VDJ][A-R0-9]+[-/\w]*))'
                         r'(\*(?P<allele>[.\w]+))?')


def parseCall(call):
    """
    Parse one gene call cell into the genes it names.

    Arguments:
      call (str): the cell, e.g. 'IGHV1-2*02, IGHV1-2*04'.

    Returns:
      frozenset: the gene names, alleles dropped (e.g. {'IGHV1-2'}); empty
      when the cell is blank or names no recognizable gene.
    """
    return frozenset(match.group('gene') for match in _CALL_REGEX.finditer(call or ''))


def callsAgree(a, b):
    """
    Whether two parsed calls (see parseCall) agree.

    Arguments:
      a (frozenset): one parsed call.
      b (frozenset): the other.

    Returns:
      bool: True if they name a gene in common. A blank or unrecognized call
      (empty) agrees with nothing.
    """
    return bool(a & b)


class GeneGate:
    """
    The gene agreement test between every query row and every reference
    row, over one or more gene call columns -- a pair passes only if it
    agrees (callsAgree) on each of them.
    """

    def __init__(self, queries, reference, gene_cols):
        """
        Arguments:
          queries (pandas.DataFrame): the query rows, in result order.
          reference (pandas.DataFrame): the reference rows, by position.
          gene_cols (list): gene call columns to gate on, e.g.
            ['v_call', 'j_call'].
        """
        self.query_calls = [self._parseColumn(queries, col, 'query')
                            for col in gene_cols]
        self.ref_calls = [self._parseColumn(reference, col, 'reference')
                          for col in gene_cols]

    @staticmethod
    def _parseColumn(frame, col, label):
        """
        Returns:
          list: parseCall of every value of `frame[col]`, in row order --
          each distinct value parsed once.
        """
        values = frame[col].fillna('').astype(str).tolist()
        parsed = {value: parseCall(value) for value in set(values)}
        unparsed = sum(1 for value in values if value.strip() and not parsed[value])
        if unparsed:
            log.warning("%d %s %s values name no recognizable gene; those rows "
                        'agree with no gene', unparsed, label, col)

        return [parsed[value] for value in values]

    def agrees(self, query_row, ref_row):
        """
        Arguments:
          query_row (int): a query position.
          ref_row (int): a reference position.

        Returns:
          bool: whether the pair agrees on every gene column.
        """
        return all(callsAgree(query_calls[query_row], ref_calls[ref_row])
                   for query_calls, ref_calls in zip(self.query_calls, self.ref_calls))

    def keep(self, query_rows, ref_rows):
        """
        Arguments:
          query_rows (numpy.ndarray): query positions, one per pair.
          ref_rows (numpy.ndarray): reference positions, one per pair.

        Returns:
          numpy.ndarray: one bool per pair, True where it agrees.
        """
        return numpy.fromiter((self.agrees(q, r) for q, r in zip(query_rows.tolist(),
                                                                 ref_rows.tolist())),
                              dtype=bool, count=len(query_rows))


def buildGeneGate(queries, reference, gene_cols):
    """
    Returns:
      GeneGate: for `gene_cols`, or None when `gene_cols` is empty or None
      -- no gene gate at all.
    """
    if not gene_cols:
        return None

    return GeneGate(queries, reference, gene_cols)
