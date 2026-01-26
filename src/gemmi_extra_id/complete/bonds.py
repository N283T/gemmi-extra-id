"""Polymer bond detection for complete mode.

This module provides functions to detect and infer polymer bonds
for AtomWorks-compatible entity assignment.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Final

if TYPE_CHECKING:
    from collections.abc import Sequence

# Import ResKey from cif_utils to avoid duplication
from gemmi_extra_id.complete.cif_utils import ResKey

# =============================================================================
# Chemical Component Type Constants (from AtomWorks)
# =============================================================================

# fmt: off
AA_LIKE_CHEM_TYPES: Final[frozenset[str]] = frozenset([
    chemtype.upper() for chemtype in (
        "D-beta-peptide, C-gamma linking", "D-gamma-peptide, C-delta linking",
        "D-peptide COOH carboxy terminus", "D-peptide NH3 amino terminus", "D-peptide linking",
        "L-beta-peptide, C-gamma linking", "L-gamma-peptide, C-delta linking",
        "L-peptide COOH carboxy terminus", "L-peptide NH3 amino terminus", "L-peptide linking",
        "peptide linking", "peptide-like",
    )
])
"""Set of amino acid-like chemical component types. All uppercase."""

RNA_LIKE_CHEM_TYPES: Final[frozenset[str]] = frozenset([
    chemtype.upper() for chemtype in (
        "L-RNA linking", "RNA OH 3 prime terminus", "RNA OH 5 prime terminus", "RNA linking"
    )
])
"""Set of RNA-like chemical component types. All uppercase."""

DNA_LIKE_CHEM_TYPES: Final[frozenset[str]] = frozenset([
    chemtype.upper() for chemtype in (
        "DNA OH 3 prime terminus", "DNA OH 5 prime terminus", "DNA linking", "L-DNA linking"
    )
])
"""Set of DNA-like chemical component types. All uppercase."""

NA_LIKE_CHEM_TYPES: Final[frozenset[str]] = RNA_LIKE_CHEM_TYPES | DNA_LIKE_CHEM_TYPES
"""DNA or RNA-like chemical component types."""

CHEM_TYPE_POLYMERIZATION_ATOMS: Final[dict[str, tuple[str, str]]] = {
    chemtype.upper(): atoms for chemtype, atoms in {
        # peptide bonds
        "peptide linking": ("C", "N"),
        "L-peptide linking": ("C", "N"),
        "D-peptide linking": ("C", "N"),
        "L-beta-peptide, C-gamma linking": ("CG", "N"),
        "D-beta-peptide, C-gamma linking": ("CG", "N"),
        "L-gamma-peptide, C-delta linking": ("CD", "N"),
        "D-gamma-peptide, C-delta linking": ("CD", "N"),
        # phosphodiester bonds
        "DNA linking": ("O3'", "P"),
        "L-DNA linking": ("O3'", "P"),
        "RNA linking": ("O3'", "P"),
        "L-RNA linking": ("O3'", "P"),
    }.items()
}
"""A mapping of chemical component types to the atoms that link when part of a polymer."""
# fmt: on


def detect_sequence_gaps(
    residues: Sequence[ResKey],
) -> list[int]:
    """Detect gaps in a sequence by checking residue ID differences.

    A gap is detected when the difference between consecutive residue
    sequence numbers is greater than 1. This matches AtomWorks behavior.

    Args:
        residues: List of (seq_id, ins_code) tuples in sequence order.

    Returns:
        List of indices where gaps occur (i.e., no bond between residues[i] and residues[i+1]).

    Example:
        >>> residues = [("1", "."), ("2", "."), ("5", "."), ("6", ".")]
        >>> detect_sequence_gaps(residues)
        [1]  # Gap between index 1 (seq 2) and index 2 (seq 5)
    """
    gaps: list[int] = []

    for i in range(len(residues) - 1):
        seq_a, _ = residues[i]
        seq_b, _ = residues[i + 1]

        # Try to parse as integers
        try:
            seq_a_int = int(seq_a)
            seq_b_int = int(seq_b)
            if seq_b_int - seq_a_int > 1:
                gaps.append(i)
        except ValueError:
            # Non-numeric sequence IDs - assume no gap
            # (this is a conservative approach)
            pass

    return gaps


def get_inferred_polymer_bonds(
    chain_residues: list[ResKey],
    residue_chem_types: dict[ResKey, str] | None = None,
) -> list[tuple[ResKey, ResKey]]:
    """Infer polymer bonds from consecutive residues.

    This function infers bonds between consecutive residues in a polymer chain,
    respecting sequence gaps. It matches AtomWorks behavior.

    Args:
        chain_residues: List of (seq_id, ins_code) tuples for a single chain,
            in sequence order.
        residue_chem_types: Optional mapping from ResKey to chem_comp_type.
            Used to determine polymerization atoms. If None, assumes standard
            peptide/nucleic acid bonds.

    Returns:
        List of (ResKey, ResKey) tuples representing inferred bonds.

    Example:
        >>> residues = [("1", "."), ("2", "."), ("3", ".")]
        >>> get_inferred_polymer_bonds(residues)
        [(("1", "."), ("2", ".")), (("2", "."), ("3", "."))]
    """
    if len(chain_residues) < 2:
        return []

    # Detect gaps in sequence
    gaps = set(detect_sequence_gaps(chain_residues))

    # Infer bonds between consecutive residues (except at gaps)
    bonds: list[tuple[ResKey, ResKey]] = []

    for i in range(len(chain_residues) - 1):
        if i in gaps:
            # Gap detected - no bond
            continue

        res_a = chain_residues[i]
        res_b = chain_residues[i + 1]
        bonds.append((res_a, res_b))

    return bonds


def get_polymerization_atoms(
    chem_type: str,
) -> tuple[str, str] | None:
    """Get the atom names that form polymer bonds for a given chem_comp_type.

    Args:
        chem_type: Chemical component type (e.g., "L-PEPTIDE LINKING").

    Returns:
        Tuple of (leaving_atom, joining_atom) names, or None if not a polymer type.
        For peptides: ("C", "N") meaning C of residue i bonds to N of residue i+1.
        For nucleic acids: ("O3'", "P") meaning O3' of residue i bonds to P of residue i+1.
    """
    return CHEM_TYPE_POLYMERIZATION_ATOMS.get(chem_type.upper())


def is_aa_like(chem_type: str) -> bool:
    """Check if a chem_comp_type is amino acid-like."""
    return chem_type.upper() in AA_LIKE_CHEM_TYPES


def is_na_like(chem_type: str) -> bool:
    """Check if a chem_comp_type is nucleic acid-like (DNA or RNA)."""
    return chem_type.upper() in NA_LIKE_CHEM_TYPES


def is_rna_like(chem_type: str) -> bool:
    """Check if a chem_comp_type is RNA-like."""
    return chem_type.upper() in RNA_LIKE_CHEM_TYPES


def is_dna_like(chem_type: str) -> bool:
    """Check if a chem_comp_type is DNA-like."""
    return chem_type.upper() in DNA_LIKE_CHEM_TYPES


__all__ = [
    "AA_LIKE_CHEM_TYPES",
    "RNA_LIKE_CHEM_TYPES",
    "DNA_LIKE_CHEM_TYPES",
    "NA_LIKE_CHEM_TYPES",
    "CHEM_TYPE_POLYMERIZATION_ATOMS",
    "detect_sequence_gaps",
    "get_inferred_polymer_bonds",
    "get_polymerization_atoms",
    "is_aa_like",
    "is_na_like",
    "is_rna_like",
    "is_dna_like",
]
