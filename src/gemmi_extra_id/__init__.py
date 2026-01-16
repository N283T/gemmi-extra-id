"""
gemmi_extra_id - Assign extra IDs to mmCIF files using Gemmi.

This package assigns AtomWorks-compatible IDs to mmCIF files:
- molecule_id: Connected component of covalent bond graph
- pn_unit_id: Same-type chains covalently linked
- chain_entity, pn_unit_entity, molecule_entity: Entity IDs

Example:
    >>> from gemmi_extra_id import assign_molecule_id
    >>> mapping = assign_molecule_id("input.cif", "output.cif")
    >>> print(mapping)
    OrderedDict([('A', 0), ('B', 0), ('C', 1)])

Extended example with all IDs:
    >>> from gemmi_extra_id import assign_extended_ids
    >>> result = assign_extended_ids("input.cif", "output.cif")
    >>> for chain, info in result.chain_info.items():
    ...     print(f"{chain}: entity={info.entity_id}, pn_unit={info.pn_unit_id}")
"""

from gemmi_extra_id.formatters import (
    to_csv,
    to_extended_csv,
    to_extended_json,
    to_extended_table,
    to_extended_tsv,
    to_json,
    to_table,
    to_tsv,
)
from gemmi_extra_id.graph import DEFAULT_COVALENT_TYPES, find_components, find_pn_units
from gemmi_extra_id.mmcif import (
    VALID_SWAP_TARGETS,
    AssignmentResult,
    ChainInfo,
    assign_extended_ids,
    assign_molecule_id,
    swap_auth_asym_id,
)

__version__ = "0.1.0"
__all__ = [
    # Core functions
    "assign_molecule_id",
    "assign_extended_ids",
    "swap_auth_asym_id",
    "find_components",
    "find_pn_units",
    # Data classes
    "ChainInfo",
    "AssignmentResult",
    # Constants
    "DEFAULT_COVALENT_TYPES",
    "VALID_SWAP_TARGETS",
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
