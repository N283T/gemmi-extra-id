"""
cifmolid - Assign molecule_id to mmCIF files based on covalent connectivity.

This package identifies connected components of chains in mmCIF files
using covalent bond information from _struct_conn.

Example:
    >>> from cifmolid import assign_molecule_id
    >>> mapping = assign_molecule_id("input.cif", "output.cif")
    >>> print(mapping)
    OrderedDict([('A', 0), ('B', 0), ('C', 1)])

Extended example with all IDs:
    >>> from cifmolid import assign_extended_ids
    >>> result = assign_extended_ids("input.cif", "output.cif")
    >>> for chain, info in result.chain_info.items():
    ...     print(f"{chain}: entity={info.entity_id}, pn_unit={info.pn_unit_id}")
"""

from cifmolid.formatters import (
    to_csv,
    to_extended_csv,
    to_extended_json,
    to_extended_table,
    to_extended_tsv,
    to_json,
    to_table,
    to_tsv,
)
from cifmolid.graph import DEFAULT_COVALENT_TYPES, find_components, find_pn_units
from cifmolid.mmcif import AssignmentResult, ChainInfo, assign_extended_ids, assign_molecule_id

__version__ = "0.1.0"
__all__ = [
    # Core functions
    "assign_molecule_id",
    "assign_extended_ids",
    "find_components",
    "find_pn_units",
    # Data classes
    "ChainInfo",
    "AssignmentResult",
    # Constants
    "DEFAULT_COVALENT_TYPES",
    # Output formatters (basic)
    "to_json",
    "to_csv",
    "to_tsv",
    "to_table",
    # Output formatters (extended)
    "to_extended_json",
    "to_extended_csv",
    "to_extended_tsv",
    "to_extended_table",
    # Metadata
    "__version__",
]
