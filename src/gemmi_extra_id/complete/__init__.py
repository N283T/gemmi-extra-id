"""Complete mode: AtomWorks-compatible entity assignment.

This module implements AtomWorks-compatible entity ID assignment
using Weisfeiler-Lehman graph hashing with full residue-level
bond detection including inferred polymer bonds.

The complete mode produces entity IDs that exactly match AtomWorks output,
unlike the default loose mode which uses a faster approximation.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Set as AbstractSet

    from gemmi_extra_id.mmcif import AssignmentResult


def assign_extended_ids_complete(
    input_path: str | Path,
    output_path: str | Path | None = None,
    covalent_types: AbstractSet[str] | None = None,
) -> AssignmentResult:
    """Assign extended IDs using AtomWorks-compatible algorithm.

    This function implements the complete entity assignment algorithm
    that matches AtomWorks output exactly. It includes:
    - Inferred polymer bonds from coordinates
    - Weisfeiler-Lehman graph hashing at residue level
    - Inter-level bond hashing

    Args:
        input_path: Path to input mmCIF file.
        output_path: Path to output mmCIF file. If None, no file is written.
        covalent_types: Set of conn_type_id values to treat as covalent bonds.
            Defaults to {"covale", "disulf"}.

    Returns:
        AssignmentResult containing ChainInfo for each chain.

    Raises:
        NotImplementedError: Complete mode is not yet implemented.
        ImportError: If networkx is not installed.
    """
    # Check for networkx (required for complete mode)
    try:
        import networkx  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "Complete mode requires networkx. Install with: pip install gemmi-extra-id[complete]"
        ) from e

    raise NotImplementedError(
        "Complete mode is not yet implemented. "
        "Use default loose mode for now, or use --hash-entities for partial compatibility."
    )


__all__ = ["assign_extended_ids_complete"]
