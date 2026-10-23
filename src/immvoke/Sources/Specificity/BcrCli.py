"""
`immvoke specificity bcr` -- the no-verb shortcut for fetching IEDB BCR

Sugar for `immvoke specificity iedb download bcr`: BCR is currently the
only table this project annotates against, so this drops both the source
name and the `download` verb -- naming the table is enough. Calling it
fetches IEDB's AIRR export and writes its BCR table as an AIRR rearrangement
TSV that keeps every column IEDB provides -- the epitope and assay context
included -- with every field the AIRR schema requires present, at
`<outdir>/iedb/bcr.tsv`.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
from argparse import Namespace
from pathlib import Path

# Immvoke imports
from immvoke.Commandline import CommonHelpFormatter

#: Written into when no --outdir is given; the table lands in its iedb/
#: folder (see immvoke.Cli.handleSpecificityDownload).
DEFAULT_OUTDIR = Path('databases')


def addBcrParser(dbs):
    """
    Add the `bcr` shortcut to the `specificity` subcommand.

    Arguments:
      dbs (argparse._SubParsersAction): the `specificity` subparsers object
        (its `db` level), as built by
        `immvoke.Cli._addSpecificityGroupParser`.
    """
    parser = dbs.add_parser(
        'bcr', help='fetch IEDB BCR receptor sequences as AIRR',
        description='Download IEDB\'s BCR (antibody) receptor table as an '
                    'AIRR rearrangement TSV, keeping IEDB\'s own columns '
                    '(epitope, assay, ...) alongside every AIRR required '
                    'field -- shorthand for `immvoke specificity iedb '
                    'download bcr`.',
        formatter_class=CommonHelpFormatter)
    parser.add_argument('--outdir', type=Path, default=DEFAULT_OUTDIR,
                        help='directory to write into; the table lands at '
                             'iedb/bcr.tsv inside it')


def handleBcr(args):
    """
    Fetch IEDB's `bcr` table and write it as AIRR, keeping every column.

    Delegates to `immvoke.Cli.handleSpecificityDownload` rather than
    duplicating it -- imported locally to avoid a circular import, since
    `immvoke.Cli` is what wires this module in.

    Arguments:
      args (argparse.Namespace): parsed arguments, as built by addBcrParser.

    Returns:
      int: process exit status.
    """
    from immvoke.Cli import handleSpecificityDownload

    return handleSpecificityDownload(Namespace(
        db='iedb', table='bcr', outdir=args.outdir, limit=None,
        dry_run=False, no_resume=False, strict_airr=False))
