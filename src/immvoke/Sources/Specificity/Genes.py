"""
V/J-gene call parsing and comparison for the annotation gene gate

An AIRR `v_call`/`j_call` cell is not always one clean name: IEDB alone mixes
allele level calls (`IGHV1-2*02`), gene level ones (`IGHV1-2`), several calls
in one cell (`IGHV1-2*02, IGHV1-2*04`) and blanks. Comparing those as plain
strings makes `IGHV1-2*02` and `IGHV1-2` disagree, so the gene gate parses
each cell into the names it holds and compares those instead, at a level the
caller chooses:

    allele  same gene and same allele (IGHV1-2*02 == IGHV1-2*02)
    gene    same gene, any allele     (IGHV1-2*02 == IGHV1-2*04)
    family  same family               (IGHV1-2*02 == IGHV1-69*01)

Two cells agree when any name in one agrees with any name in the other. A
blank cell, or one naming no recognizable gene, has nothing to compare and
agrees with nothing -- not even another blank. A gene level name (`IGHV1-2`) carries no allele, so under 'allele'
it agrees with nothing -- it is a gene level match, not an allele level one.

The gene and family patterns are Change-O's (changeo/Gene.py, Immcantation),
so a name parses here the way Change-O's getGene/getFamily would parse it.
The allele is read as an optional `*` suffix of the gene instead of with
Change-O's allele pattern, which also accepts `-` as the allele separator
and so would read a gene level `IGHV1-2` as gene `IGHV1` allele `2`.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging
import re

import numpy

# Immvoke imports
from immvoke.Exceptions import AnnotationError

log = logging.getLogger(__name__)

#: Levels a gene name can be compared at, strictest first.
GENE_MATCH_LEVELS = ('allele', 'gene', 'family')

#: Change-O's gene pattern, plus an optional `*allele` suffix.
_CALL_REGEX = re.compile(r'(?P<gene>(IG[HLK]|TR[ABGD])([VDJ][A-R0-9]+[-/\w]*))'
                         r'(\*(?P<allele>[.\w]+))?')
#: Change-O's family pattern.
_FAMILY_REGEX = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-R0-9]+))')


def _checkLevel(level):
    """
    Raises:
      AnnotationError: if `level` is not one of GENE_MATCH_LEVELS.
    """
    if level not in GENE_MATCH_LEVELS:
        raise AnnotationError('unknown gene_match %r; choose one of: %s'
                              % (level, ', '.join(GENE_MATCH_LEVELS)))


def parseCall(call, level='gene'):
    """
    Parse one gene call cell into the names it holds at `level`.

    Arguments:
      call (str): the cell, e.g. 'IGHV1-2*02, IGHV1-2*04'.
      level (str): one of GENE_MATCH_LEVELS.

    Returns:
      dict: name (the gene, or the family under 'family') to the set of its
      alleles -- None standing for "any allele" at 'gene'/'family' level; at
      'allele' level a gene level name keeps its key with no allele, so it
      can never agree. Empty when the cell is blank or names no recognizable
      gene.

    Raises:
      AnnotationError: if `level` is not recognized.
    """
    _checkLevel(level)
    names = {}
    for match in _CALL_REGEX.finditer(call or ''):
        name, allele = match.group('gene'), match.group('allele')
        if level == 'family':
            name = _FAMILY_REGEX.match(name).group(0)
        alleles = names.setdefault(name, set())
        if level != 'allele':
            alleles.add(None)
        elif allele is not None:
            alleles.add(allele)

    return names


def _loci(call):
    """
    Returns:
      set: the loci (IGH, IGK, ...) the gene names in one call cell belong to.
    """
    return {match.group(2) for match in _CALL_REGEX.finditer(call or '')}


def inferLocus(v_call, j_call):
    """
    The locus a chain's gene calls name, for a chain whose locus is blank.

    Arguments:
      v_call (str): the V-gene call cell.
      j_call (str): the J-gene call cell.

    Returns:
      str: the one locus the V call names (IGKV1-39 -> IGK), or, without a
      V call, the one the J call names; '' when the calls name none, more
      than one, or V and J disagree -- nothing is guessed.
    """
    v_loci, j_loci = _loci(v_call), _loci(j_call)
    if len(v_loci) == 1 and (not j_loci or j_loci == v_loci):
        return v_loci.pop()
    if not v_loci and len(j_loci) == 1:
        return j_loci.pop()

    return ''


def callsAgree(a, b):
    """
    Whether two parsed calls (see parseCall) agree.

    Arguments:
      a (dict): one parsed call.
      b (dict): the other, parsed at the same level.

    Returns:
      bool: True if they share a name with an allele in common -- any
      allele (None) at 'gene'/'family' level. A blank or unrecognized call
      (empty) agrees with nothing.
    """
    if not a or not b:
        return False

    for name in a.keys() & b.keys():
        alleles_a, alleles_b = a[name], b[name]
        if None in alleles_a or None in alleles_b or alleles_a & alleles_b:
            return True

    return False


class GeneGate:
    """
    The gene agreement test between every query row and every reference
    row, over one or more gene call columns -- a pair passes only if it
    agrees (callsAgree) on each of them.
    """

    def __init__(self, queries, reference, gene_cols, level='gene'):
        """
        Arguments:
          queries (pandas.DataFrame): the query rows, in result order.
          reference (pandas.DataFrame): the reference rows, by position.
          gene_cols (list): gene call columns to gate on, e.g.
            ['v_call', 'j_call'].
          level (str): one of GENE_MATCH_LEVELS.

        Raises:
          AnnotationError: if `level` is not recognized.
        """
        _checkLevel(level)
        self.query_calls = [self._parseColumn(queries, col, level, 'query')
                            for col in gene_cols]
        self.ref_calls = [self._parseColumn(reference, col, level, 'reference')
                          for col in gene_cols]

    @staticmethod
    def _parseColumn(frame, col, level, label):
        """
        Returns:
          list: parseCall of every value of `frame[col]`, in row order --
          each distinct value parsed once.
        """
        values = frame[col].fillna('').astype(str).tolist()
        parsed = {value: parseCall(value, level) for value in set(values)}
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


def buildGeneGate(queries, reference, gene_cols, level='gene'):
    """
    Returns:
      GeneGate: for `gene_cols` at `level`, or None when `gene_cols` is
      empty or None -- no gene gate at all.
    """
    if not gene_cols:
        _checkLevel(level)
        return None

    return GeneGate(queries, reference, gene_cols, level)
