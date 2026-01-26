"""CCD (Chemical Component Dictionary) reader for bond information.

This module provides functions to read bond information from CCD files.
The bond information is used for WL graph hashing to compute entity IDs.
"""

from __future__ import annotations

import functools
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi

if TYPE_CHECKING:
    from collections.abc import Sequence

logger = logging.getLogger(__name__)

# Bond type mapping: CCD value_order -> integer
BOND_ORDER_MAP: dict[str, int] = {
    "SING": 1,
    "DOUB": 2,
    "TRIP": 3,
    "AROM": 5,  # Aromatic
    "DELO": 1,  # Delocalized (treat as single)
    "QUAD": 4,  # Quadruple (rare)
}

# Default CCD mirror path (can be set via environment variable)
DEFAULT_CCD_PATH = os.environ.get("CCD_MIRROR_PATH", "")


@dataclass(frozen=True)
class CCDBond:
    """A bond in a CCD component."""

    atom1: str
    atom2: str
    order: int  # 1=single, 2=double, 3=triple, 5=aromatic


@dataclass(frozen=True)
class CCDComponent:
    """A chemical component from CCD with atoms and bonds."""

    code: str
    atoms: tuple[str, ...]  # atom names
    elements: tuple[str, ...]  # element symbols
    bonds: tuple[CCDBond, ...]


def _get_ccd_path(code: str, ccd_path: str) -> Path:
    """Get path to CCD file for a component code."""
    return Path(ccd_path) / code[0] / code / f"{code}.cif"


@functools.lru_cache(maxsize=1000)
def load_ccd_component(code: str, ccd_path: str = DEFAULT_CCD_PATH) -> CCDComponent | None:
    """Load a CCD component from file.

    Args:
        code: The CCD component code (e.g., "ALA", "GLY").
        ccd_path: Path to the CCD mirror directory.

    Returns:
        CCDComponent with atoms and bonds, or None if not found.
    """
    if not ccd_path:
        logger.debug(f"CCD path not set, skipping {code}")
        return None

    path = _get_ccd_path(code, ccd_path)
    if not path.exists():
        logger.debug(f"CCD file not found: {path}")
        return None

    try:
        doc = gemmi.cif.read(str(path))
        block = doc[0]
    except Exception as e:
        logger.warning(f"Failed to read CCD file {path}: {e}")
        return None

    # Extract atoms
    atom_loop = block.find_loop("_chem_comp_atom.atom_id")
    element_loop = block.find_loop("_chem_comp_atom.type_symbol")

    if not atom_loop or not element_loop:
        logger.debug(f"No atom data in CCD {code}")
        return None

    atoms = tuple(str(a) for a in atom_loop)
    elements = tuple(str(e) for e in element_loop)

    # Extract bonds
    bond_atom1_loop = block.find_loop("_chem_comp_bond.atom_id_1")
    bond_atom2_loop = block.find_loop("_chem_comp_bond.atom_id_2")
    bond_order_loop = block.find_loop("_chem_comp_bond.value_order")

    bonds: list[CCDBond] = []
    if bond_atom1_loop and bond_atom2_loop and bond_order_loop:
        for a1, a2, order in zip(bond_atom1_loop, bond_atom2_loop, bond_order_loop, strict=True):
            bond_order = BOND_ORDER_MAP.get(str(order), 1)
            bonds.append(CCDBond(str(a1), str(a2), bond_order))

    return CCDComponent(
        code=code,
        atoms=atoms,
        elements=elements,
        bonds=tuple(bonds),
    )


def get_residue_bonds(
    res_name: str,
    atom_names: Sequence[str],
    ccd_path: str = DEFAULT_CCD_PATH,
) -> list[tuple[int, int, int]]:
    """Get bonds for a residue based on CCD information.

    Args:
        res_name: The residue name (e.g., "ALA").
        atom_names: The atom names in the residue.
        ccd_path: Path to the CCD mirror directory.

    Returns:
        List of (atom1_idx, atom2_idx, bond_order) tuples.
        Indices are relative to the atom_names list.
    """
    component = load_ccd_component(res_name, ccd_path)
    if component is None:
        return []

    # Create atom name to index mapping
    name_to_idx: dict[str, int] = {}
    for i, name in enumerate(atom_names):
        if name not in name_to_idx:
            name_to_idx[name] = i

    # Map CCD bonds to local indices
    bonds: list[tuple[int, int, int]] = []
    for ccd_bond in component.bonds:
        idx1 = name_to_idx.get(ccd_bond.atom1)
        idx2 = name_to_idx.get(ccd_bond.atom2)

        if idx1 is not None and idx2 is not None:
            bonds.append((idx1, idx2, ccd_bond.order))

    return bonds


def clear_cache() -> None:
    """Clear the CCD component cache."""
    load_ccd_component.cache_clear()
