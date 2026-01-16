"""
cifmolid - Assign molecule_id to mmCIF files based on covalent connectivity.

This package identifies connected components of chains in mmCIF files
using covalent bond information from _struct_conn.

Example:
    >>> from cifmolid import assign_molecule_id
    >>> mapping = assign_molecule_id("input.cif", "output.cif")
    >>> print(mapping)
    OrderedDict([('A', 0), ('B', 0), ('C', 1)])
"""

from cifmolid.cif import assign_molecule_id
from cifmolid.core import DEFAULT_COVALENT_TYPES, find_components
from cifmolid.output import to_csv, to_json, to_table, to_tsv

__version__ = "0.1.0"
__all__ = [
    "assign_molecule_id",
    "find_components",
    "DEFAULT_COVALENT_TYPES",
    "to_json",
    "to_csv",
    "to_tsv",
    "to_table",
    "__version__",
]
