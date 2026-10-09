"""
Specificity source registry

A second, separate registry from immvoke.Sources: specificity databases
(IEDB, CoV-AbDab, ASD, and others to follow) are SourceBase subclasses too,
but are kept apart from repertoire sources so the commandline can mount them
under their own `immvoke specificity <db>` tree instead of `immvoke
<source>`.
"""

# Info
__author__ = 'Pramod Shinde'

# Immvoke imports
from immvoke.Sources.Specificity.Asd import AsdSource
from immvoke.Sources.Specificity.CovAbDab import CovAbDabSource
from immvoke.Sources.Specificity.Iedb import IedbSource

#: Every specificity source immvoke knows about, by commandline name.
REGISTRY = {IedbSource.name: IedbSource, CovAbDabSource.name: CovAbDabSource,
           AsdSource.name: AsdSource}


def getSpecificitySource(name, client, schema=None):
    """
    Instantiate a specificity source by name.

    Arguments:
      name (str): the source name.
      client (HttpClient): the shared HTTP client.
      schema (SourceSchema): a preloaded snapshot, or None to load on demand.

    Returns:
      SourceBase: the source.

    Raises:
      KeyError: if the name is not registered.
    """
    if name not in REGISTRY:
        raise KeyError("unknown specificity source '%s'; known: %s"
                       % (name, ', '.join(sorted(REGISTRY))))

    return REGISTRY[name](client, schema=schema)
