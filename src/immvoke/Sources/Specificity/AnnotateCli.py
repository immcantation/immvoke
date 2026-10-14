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
(e.g. `immvoke specificity bcr`'s output) as `--reference` -- and writes
one antigen specificity annotation per query row. It does not itself
download or convert anything; that is handled by `immvoke specificity
bcr`.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
from pathlib import Path

# Immvoke imports
from immvoke.Commandline import CommonHelpFormatter
from immvoke.Sources.Specificity.Annotate import BACKENDS, annotateSpecificity
from immvoke.Sources.Specificity.Genes import GENE_MATCH_LEVELS


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
                    'specificity bcr`\'s output.',
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
                        choices=['CDR3', 'VDJ', 'CDRs'], default='CDR3',
                        help='which sequence to compare: CDR3 (cdr3_aa), VDJ '
                             '(sequence_aa, the full V(D)J amino acid '
                             'sequence) or CDRs (cdr1_aa, cdr2_aa and cdr3_aa '
                             'pooled -- --distance-metric hamming only)')
    parser.add_argument('--partition-by', dest='partition_by',
                        choices=['V', 'VJ', 'VJL'], default=None,
                        help='only compare row pairs with the same V-gene '
                             '(V), the same V- and J-gene (VJ), or the same '
                             'V-gene, J-gene and CDR3 length (VJL); unset, '
                             'every row pair is compared')
    parser.add_argument('--gene-match', dest='gene_match',
                        choices=GENE_MATCH_LEVELS, default='gene',
                        help='how --partition-by and --paired compare v_call/'
                             'j_call: allele (same gene and allele), gene '
                             '(same gene, any allele) or family. A cell with '
                             'several calls agrees if any one does; a blank '
                             'or unrecognized call is unknown and agrees '
                             'with any gene')
    parser.add_argument('--paired', action='store_true',
                        help='require V-gene and J-gene agreement (see '
                             '--gene-match) and exact CDR3 agreement on both '
                             'the heavy and the light chain '
                             'of the same reference antibody, and classify '
                             'each query cell_id as both/heavy_only/'
                             'light_only/conflicting/neither. Without it, '
                             'only heavy chain (locus IGH) rows are compared. '
                             '--sequence-region, --partition-by, '
                             '--identity-threshold and --distance-metric '
                             'have no effect here')
    parser.add_argument('--identity-threshold', dest='identity_threshold',
                        type=float, default=None, metavar='FRACTION',
                        help='minimum sequence identity in (0.0, 1.0], e.g. '
                             '0.9 allows up to 10%% of the query\'s positions '
                             'to differ; unset, an exact match is required')
    parser.add_argument('--distance-metric', dest='distance_metric',
                        choices=['hamming', 'levenshtein'], default='hamming',
                        help='hamming (substitutions only) or levenshtein '
                             '(substitutions, insertions and deletions -- '
                             'not with --sequence-region CDRs). Only matters '
                             'when --identity-threshold is below 1.0')
    parser.add_argument('--backend', choices=BACKENDS, default='symscan',
                        help='matching engine. symscan is a required '
                             'dependency and runs every mode; python runs '
                             'only when chosen here, never as a fallback, '
                             'and gives the same result. Rows symscan cannot '
                             'index (an '
                             '\'X\' wildcard, a distance above 14, or very '
                             'long sequences at a high tolerance) are '
                             'compared in Python automatically')
    parser.add_argument('--threads', type=int, default=1,
                        help='threads for symscan\'s thread pool (-1 for '
                             'every core; only the first symscan run in a '
                             'process sets it); the python backend is '
                             'single threaded')
    parser.add_argument('-o', '--out', type=Path, required=True,
                        help='write the annotation table here')


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
        distance_metric=args.distance_metric, gene_match=args.gene_match,
        threads=args.threads,
        output=args.out, backend=args.backend)

    return 0
