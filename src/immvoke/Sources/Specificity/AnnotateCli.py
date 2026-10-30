"""
`immvoke annotate` -- the commandline front end for
Annotate.annotateSpecificity()

Kept apart from immvoke.Cli (which only imports addAnnotateParser/
handleAnnotate and wires them into the top level parser and dispatcher) so
that the annotation step -- argument definitions, help text, and the handler
itself -- lives alongside the rest of the specificity annotation code it
operates on. Every flag maps to one annotateSpecificity() parameter, and
handleAnnotate only forwards them.

Reads two AIRR rearrangement TSVs already on disk -- typically a user's own
receptor table as `--query` and a specificity database's downloaded table
(e.g. `immvoke specificity iedb download bcr`'s output) as `--reference` --
and writes the annotation (each query row's best hits) and a record of the
run into `--outdir`, plus every pair the filters admit and how close it
is with `--full-dist-matrix`. It does not itself download or convert anything; that
is handled by `immvoke specificity iedb download bcr`.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
from pathlib import Path

# Immvoke imports
from immvoke.Commandline import CommonHelpFormatter
from immvoke.Sources.Specificity.Annotate import (
    ANNOTATION_FILE,
    ANNOTATION_METADATA,
    DEFAULT_COUNT_FIELD,
    FULL_MATRIX_FILE,
    annotateSpecificity,
)


def addAnnotateParser(commands):
    """
    Add the `annotate` subcommand to the top level parser.

    Arguments:
      commands (argparse._SubParsersAction): the top level subparsers
        object, as built by `immvoke.Cli.getArgParser`.
    """
    parser = commands.add_parser(
        'annotate', help='annotate a receptor table against a reference table',
        description='Annotate each row of a query AIRR rearrangement table '
                    'with antigen specificity information from a fixed '
                    'reference AIRR rearrangement table, e.g. `immvoke '
                    'specificity iedb download bcr`\'s output.',
        formatter_class=CommonHelpFormatter)

    parser.add_argument('--query', type=Path, required=True,
                        help='AIRR rearrangement TSV to annotate, one row '
                             'per chain/sequence')
    parser.add_argument('--reference', type=Path, required=True,
                        help='AIRR rearrangement TSV to annotate against, '
                             'e.g. a specificity database\'s downloaded '
                             'table; its sequence_id is reported as each '
                             'hit\'s identifier (cell_id under --paired)')
    parser.add_argument('--sequence-region', dest='sequence_region',
                        choices=['CDR3H', 'VDJ', 'CDRs'], default='CDR3H',
                        help='which sequence to compare: CDR3H (the heavy '
                             'chain CDR3, cdr3_aa), VDJ '
                             '(sequence_aa, the full V(D)J amino acid '
                             'sequence) or CDRs (cdr1_aa, cdr2_aa and cdr3_aa '
                             'pooled -- --distance-metric hamming only)')
    parser.add_argument('--partition-by', dest='partition_by',
                        choices=['V', 'VJ', 'VJL'], default=None,
                        help='only compare row pairs with the same V-gene '
                             '(V), the same V- and J-gene (VJ), or the same '
                             'V-gene, J-gene and CDR3 length (VJL); unset, '
                             'every row pair is compared. Genes match by '
                             'name, any allele (IGHV1-2*02 = IGHV1-2*04); a '
                             'cell with several calls agrees if any one does, '
                             'a blank or unrecognized call never does')
    parser.add_argument('--paired', action='store_true',
                        help='require V-gene and J-gene agreement (same '
                             'gene, any allele) and exact CDR3 agreement on both '
                             'the heavy and the light chain '
                             'of the same reference antibody, and classify '
                             'each query cell_id as both/heavy_only/'
                             'light_only/conflicting/neither. Without it, '
                             'only heavy chain (locus IGH) rows are compared. '
                             '--sequence-region, --partition-by, '
                             '--identity-threshold and --distance-metric '
                             'have no effect here')
    parser.add_argument('--count-field', dest='count_field',
                        default=DEFAULT_COUNT_FIELD, metavar='COLUMN',
                        help='--paired only: when a query cell has several '
                             'heavy (or light) chains, keep the one with the '
                             'highest value in this AIRR count column (e.g. '
                             'umi_count); on a tie, a blank count, or no such '
                             'column, the first in the file')
    parser.add_argument('--identity-threshold', dest='identity_threshold',
                        type=float, default=None, metavar='FRACTION',
                        help='minimum sequence identity in (0.0, 1.0], e.g. '
                             '0.9 allows up to 10%% of the query\'s positions '
                             'to differ; unset, an exact match is required. '
                             'Under hamming (and an exact match) an X '
                             'residue matches any residue. Rows a threshold '
                             'makes too large for symscan (over 14 '
                             'differences, or a very long sequence) are '
                             'scored directly instead, more slowly')
    parser.add_argument('--distance-metric', dest='distance_metric',
                        choices=['hamming', 'levenshtein'], default='hamming',
                        help='hamming (substitutions only) or levenshtein '
                             '(substitutions, insertions and deletions -- '
                             'not with --sequence-region CDRs). Only matters '
                             'when --identity-threshold is below 1.0')
    parser.add_argument('--threads', type=int, default=1,
                        help='threads for symscan\'s thread pool (-1 for '
                             'every core; only the first symscan run in a '
                             'process sets it)')
    parser.add_argument('--outdir', type=Path, required=True,
                        help='folder to write into, created if needed, each '
                             'file named after the reference\'s database and '
                             'table (e.g. iedb_bcr_%s): <db>_<table>_%s (each '
                             'query row\'s best hits -- the highest identity, '
                             'all of them on a tie -- with their antigen, '
                             'epitope and antigen source species) and '
                             '<db>_<table>_%s (the options, inputs and '
                             'counts of this run). A reference not downloaded '
                             'by `immvoke specificity` is named by its file '
                             'name alone'
                             % (ANNOTATION_FILE, ANNOTATION_FILE, ANNOTATION_METADATA))
    parser.add_argument('--full-dist-matrix', dest='full_dist_matrix',
                        action='store_true',
                        help='also write <db>_<table>_%s: every pair the '
                             'filters admit and how close it is, one row per '
                             '(query sequence_id, hit_id, identity). '
                             '--identity-threshold does not apply to it, '
                             '--partition-by does; under --paired, each '
                             "antibody's hits and the chains they matched on"
                             % FULL_MATRIX_FILE)


def handleAnnotate(args):
    """
    Forward parsed arguments to annotateSpecificity().

    Arguments:
      args (argparse.Namespace): parsed arguments, as built by
        addAnnotateParser.

    Returns:
      int: process exit status.
    """
    annotateSpecificity(
        query=args.query, reference=args.reference,
        sequence_region=args.sequence_region, partition_by=args.partition_by,
        paired=args.paired, identity_threshold=args.identity_threshold,
        distance_metric=args.distance_metric, threads=args.threads,
        outdir=args.outdir, count_field=args.count_field,
        full_dist_matrix=args.full_dist_matrix)

    return 0
