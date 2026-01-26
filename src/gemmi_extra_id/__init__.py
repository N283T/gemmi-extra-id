"""
gemmi_extra_id - Assign molecule_id to mmCIF files using Gemmi.

This package assigns molecule_id to chains in mmCIF files based on
covalent bond connectivity in _struct_conn.

Example:
    >>> from gemmi_extra_id import assign_molecule_id
    >>> mapping = assign_molecule_id("input.cif", "output.cif")
    >>> print(mapping)
    OrderedDict([('A', 0), ('B', 0), ('C', 1)])
"""

from gemmi_extra_id.formatters import to_json, write_output
from gemmi_extra_id.graph import DEFAULT_COVALENT_TYPES, find_components
from gemmi_extra_id.mmcif import (
    VALID_SWAP_TARGETS,
    assign_molecule_id,
    swap_auth_asym_id,
)

__version__ = "0.2.0"
__all__ = [
    # Core functions
    "assign_molecule_id",
    "swap_auth_asym_id",
    "find_components",
    # Constants
    "DEFAULT_COVALENT_TYPES",
    "VALID_SWAP_TARGETS",
    # Output formatters
    "to_json",
    "write_output",
    # Metadata
    "__version__",
]
