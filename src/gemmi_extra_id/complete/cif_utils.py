"""CIF data reading utilities for complete mode.

This module provides functions to read CIF data required for
AtomWorks-compatible entity assignment.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING

import gemmi

if TYPE_CHECKING:
    from collections.abc import Set as AbstractSet


# Mapping of chem_comp types to their polymerization atoms (atom1, atom2)
# where atom1 in residue N connects to atom2 in residue N+1.
# These bonds are standard polymer backbone bonds and should not be included
# in the inter-level bond hash. From AtomWorks CHEM_TYPE_POLYMERIZATION_ATOMS.
CHEM_TYPE_POLYMERIZATION_ATOMS: dict[str, tuple[str, str]] = {
    # Peptide bonds (C of residue N connects to N of residue N+1)
    "PEPTIDE LINKING": ("C", "N"),
    "L-PEPTIDE LINKING": ("C", "N"),
    "D-PEPTIDE LINKING": ("C", "N"),
    "L-BETA-PEPTIDE, C-GAMMA LINKING": ("CG", "N"),
    "D-BETA-PEPTIDE, C-GAMMA LINKING": ("CG", "N"),
    "L-GAMMA-PEPTIDE, C-DELTA LINKING": ("CD", "N"),
    "D-GAMMA-PEPTIDE, C-DELTA LINKING": ("CD", "N"),
    # Phosphodiester bonds (O3' of residue N connects to P of residue N+1)
    "DNA LINKING": ("O3'", "P"),
    "L-DNA LINKING": ("O3'", "P"),
    "RNA LINKING": ("O3'", "P"),
    "L-RNA LINKING": ("O3'", "P"),
}

# Type alias for residue key (seq_id, ins_code)
ResKey = tuple[str, str]


@dataclass(frozen=True)
class ChainMetadata:
    """Chain metadata from CIF for complete mode entity assignment.

    This is distinct from mmcif.ChainMetadata which contains assigned IDs.
    """

    chain_id: str  # label_asym_id
    entity_id: str
    entity_type: str  # polymer, non-polymer, water, etc.
    is_polymer: bool
    polymer_type: str | None = None  # polypeptide(L), polydeoxyribonucleotide, etc.


@dataclass
class AtomBond:
    """Atom-level bond from struct_conn."""

    chain1: str
    seq1: str
    ins1: str
    comp1: str
    atom1: str
    chain2: str
    seq2: str
    ins2: str
    comp2: str
    atom2: str
    conn_type: str


# CIF tags
_ENTITY_ID = "_entity.id"
_ENTITY_TYPE = "_entity.type"

_ENTITY_POLY_ENTITY_ID = "_entity_poly.entity_id"
_ENTITY_POLY_TYPE = "_entity_poly.type"
_ENTITY_POLY_STRAND_ID = "_entity_poly.pdbx_strand_id"

_ENTITY_POLY_SEQ_ENTITY_ID = "_entity_poly_seq.entity_id"
_ENTITY_POLY_SEQ_NUM = "_entity_poly_seq.num"
_ENTITY_POLY_SEQ_MON_ID = "_entity_poly_seq.mon_id"

_STRUCT_ASYM_ID = "_struct_asym.id"
_STRUCT_ASYM_ENTITY_ID = "_struct_asym.entity_id"

_STRUCT_CONN_TYPE = "_struct_conn.conn_type_id"
_STRUCT_CONN_P1_CHAIN = "_struct_conn.ptnr1_label_asym_id"
_STRUCT_CONN_P1_SEQ = "_struct_conn.ptnr1_label_seq_id"
_STRUCT_CONN_P1_INS = "_struct_conn.pdbx_ptnr1_PDB_ins_code"
_STRUCT_CONN_P1_COMP = "_struct_conn.ptnr1_label_comp_id"
_STRUCT_CONN_P1_ATOM = "_struct_conn.ptnr1_label_atom_id"
_STRUCT_CONN_P2_CHAIN = "_struct_conn.ptnr2_label_asym_id"
_STRUCT_CONN_P2_SEQ = "_struct_conn.ptnr2_label_seq_id"
_STRUCT_CONN_P2_INS = "_struct_conn.pdbx_ptnr2_PDB_ins_code"
_STRUCT_CONN_P2_COMP = "_struct_conn.ptnr2_label_comp_id"
_STRUCT_CONN_P2_ATOM = "_struct_conn.ptnr2_label_atom_id"


def _normalize_cif_value(value: str | None) -> str:
    """Normalize CIF missing values to '.'."""
    if value is None or value in ("?", ".", ""):
        return "."
    return value


def get_canonical_sequences(
    block: gemmi.cif.Block,
) -> dict[str, list[tuple[ResKey, str]]]:
    """Get canonical sequences from _entity_poly_seq.

    Reads the _entity_poly_seq category to get the complete canonical
    sequence for each polymer entity. This includes residues that may
    be missing from the _atom_site records.

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping entity_id to list of (ResKey, mon_id) tuples.
        ResKey is (seq_num, ".") since _entity_poly_seq doesn't have ins_code.
        The list is ordered by sequence number.

    Example:
        >>> seqs = get_canonical_sequences(block)
        >>> seqs["1"]  # Entity 1's sequence
        [(("1", "."), "MET"), (("2", "."), "ALA"), ...]
    """
    result: dict[str, list[tuple[ResKey, str]]] = {}

    # Try loop format
    entity_id_col = block.find_loop(_ENTITY_POLY_SEQ_ENTITY_ID)
    if entity_id_col:
        loop = entity_id_col.get_loop()
        if loop:
            tags = list(loop.tags)
            entity_idx = tags.index(_ENTITY_POLY_SEQ_ENTITY_ID)
            num_idx = tags.index(_ENTITY_POLY_SEQ_NUM) if _ENTITY_POLY_SEQ_NUM in tags else None
            mon_idx = (
                tags.index(_ENTITY_POLY_SEQ_MON_ID) if _ENTITY_POLY_SEQ_MON_ID in tags else None
            )

            if num_idx is not None and mon_idx is not None:
                for row_idx in range(loop.length()):
                    entity_id = _normalize_cif_value(loop[row_idx, entity_idx])
                    seq_num = _normalize_cif_value(loop[row_idx, num_idx])
                    mon_id = _normalize_cif_value(loop[row_idx, mon_idx])

                    if entity_id == "." or seq_num == "." or mon_id == ".":
                        continue

                    if entity_id not in result:
                        result[entity_id] = []

                    # ResKey uses (seq_num, ins_code); _entity_poly_seq has no ins_code
                    res_key: ResKey = (seq_num, ".")
                    result[entity_id].append((res_key, mon_id))

    return result


def get_entity_info(
    block: gemmi.cif.Block,
) -> dict[str, str]:
    """Get entity_id -> entity_type mapping from _entity.

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping entity_id to entity_type (e.g., "polymer", "non-polymer", "water").
    """
    result: dict[str, str] = {}

    # Try loop format
    entity_id_col = block.find_loop(_ENTITY_ID)
    if entity_id_col:
        loop = entity_id_col.get_loop()
        if loop:
            tags = list(loop.tags)
            id_idx = tags.index(_ENTITY_ID)
            type_idx = tags.index(_ENTITY_TYPE) if _ENTITY_TYPE in tags else None

            if type_idx is not None:
                for row_idx in range(loop.length()):
                    entity_id = _normalize_cif_value(loop[row_idx, id_idx])
                    entity_type = _normalize_cif_value(loop[row_idx, type_idx])

                    if entity_id != ".":
                        result[entity_id] = entity_type if entity_type != "." else "unknown"
    else:
        # Try single-value format
        entity_id = block.find_value(_ENTITY_ID)
        entity_type = block.find_value(_ENTITY_TYPE)
        if entity_id and entity_id not in ("?", "."):
            etype = entity_type if entity_type and entity_type not in ("?", ".") else "unknown"
            result[entity_id] = etype

    return result


def get_polymer_info(
    block: gemmi.cif.Block,
) -> dict[str, tuple[str, list[str]]]:
    """Get polymer information from _entity_poly.

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping entity_id to (polymer_type, list of chain_ids).
        polymer_type is e.g., "polypeptide(L)", "polydeoxyribonucleotide".
    """
    result: dict[str, tuple[str, list[str]]] = {}

    # Try loop format
    entity_id_col = block.find_loop(_ENTITY_POLY_ENTITY_ID)
    if entity_id_col:
        loop = entity_id_col.get_loop()
        if loop:
            tags = list(loop.tags)
            entity_idx = tags.index(_ENTITY_POLY_ENTITY_ID)
            type_idx = tags.index(_ENTITY_POLY_TYPE) if _ENTITY_POLY_TYPE in tags else None
            strand_idx = (
                tags.index(_ENTITY_POLY_STRAND_ID) if _ENTITY_POLY_STRAND_ID in tags else None
            )

            for row_idx in range(loop.length()):
                entity_id = _normalize_cif_value(loop[row_idx, entity_idx])
                poly_type = (
                    _normalize_cif_value(loop[row_idx, type_idx]) if type_idx is not None else "."
                )
                strand_id = (
                    _normalize_cif_value(loop[row_idx, strand_idx])
                    if strand_idx is not None
                    else "."
                )

                if entity_id == ".":
                    continue

                # Parse strand_id (comma-separated chain IDs)
                chains = [c.strip() for c in strand_id.split(",") if c.strip() and c.strip() != "."]

                result[entity_id] = (poly_type if poly_type != "." else "unknown", chains)
    else:
        # Try single-value format
        entity_id = block.find_value(_ENTITY_POLY_ENTITY_ID)
        poly_type = block.find_value(_ENTITY_POLY_TYPE)
        strand_id = block.find_value(_ENTITY_POLY_STRAND_ID)

        if entity_id and entity_id not in ("?", "."):
            ptype = poly_type if poly_type and poly_type not in ("?", ".") else "unknown"
            chains = []
            if strand_id and strand_id not in ("?", "."):
                chains = [c.strip() for c in strand_id.split(",") if c.strip()]
            result[entity_id] = (ptype, chains)

    return result


def get_struct_asym_mapping(
    block: gemmi.cif.Block,
) -> dict[str, str]:
    """Get chain_id -> entity_id mapping from _struct_asym.

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping chain_id (label_asym_id) to entity_id.
    """
    result: dict[str, str] = {}

    # Try loop format
    asym_id_col = block.find_loop(_STRUCT_ASYM_ID)
    if asym_id_col:
        loop = asym_id_col.get_loop()
        if loop:
            tags = list(loop.tags)
            id_idx = tags.index(_STRUCT_ASYM_ID)
            entity_idx = (
                tags.index(_STRUCT_ASYM_ENTITY_ID) if _STRUCT_ASYM_ENTITY_ID in tags else None
            )

            if entity_idx is not None:
                for row_idx in range(loop.length()):
                    asym_id = _normalize_cif_value(loop[row_idx, id_idx])
                    entity_id = _normalize_cif_value(loop[row_idx, entity_idx])

                    if asym_id != "." and entity_id != ".":
                        result[asym_id] = entity_id
    else:
        # Try single-value format
        asym_id = block.find_value(_STRUCT_ASYM_ID)
        entity_id = block.find_value(_STRUCT_ASYM_ENTITY_ID)
        if asym_id and asym_id not in ("?", ".") and entity_id and entity_id not in ("?", "."):
            result[asym_id] = entity_id

    return result


def get_chain_info_complete(
    block: gemmi.cif.Block,
) -> dict[str, ChainMetadata]:
    """Get complete chain information from CIF.

    Combines information from _entity, _entity_poly, and _struct_asym
    to build comprehensive chain information.

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping chain_id (label_asym_id) to ChainMetadata.
    """
    # Get entity type info
    entity_types = get_entity_info(block)

    # Get polymer info
    polymer_info = get_polymer_info(block)

    # Get chain -> entity mapping
    asym_mapping = get_struct_asym_mapping(block)

    # Build ChainMetadata for each chain
    result: dict[str, ChainMetadata] = {}

    for chain_id, entity_id in asym_mapping.items():
        entity_type = entity_types.get(entity_id, "unknown")
        is_polymer = entity_type == "polymer"

        polymer_type: str | None = None
        if entity_id in polymer_info:
            polymer_type = polymer_info[entity_id][0]

        result[chain_id] = ChainMetadata(
            chain_id=chain_id,
            entity_id=entity_id,
            entity_type=entity_type,
            is_polymer=is_polymer,
            polymer_type=polymer_type,
        )

    return result


def get_struct_conn_bonds(
    block: gemmi.cif.Block,
    bond_types: AbstractSet[str] | None = None,
) -> list[AtomBond]:
    """Get atom-level bonds from _struct_conn.

    Extracts bond information at the atomic level, including
    atom names and residue information. This is more detailed
    than the residue-level bonds used in loose mode.

    Args:
        block: mmCIF data block.
        bond_types: Set of conn_type_id values to include.
            Defaults to {"covale"} (covalent bonds only).

    Returns:
        List of AtomBond objects representing each bond.
    """
    if bond_types is None:
        bond_types = frozenset({"covale"})
    else:
        bond_types = frozenset(t.lower() for t in bond_types)

    result: list[AtomBond] = []

    # Find the struct_conn loop
    conn_type_col = block.find_loop(_STRUCT_CONN_TYPE)
    if not conn_type_col:
        return result

    loop = conn_type_col.get_loop()
    if not loop:
        return result

    tags = list(loop.tags)

    # Get column indices
    def get_idx(tag: str) -> int | None:
        return tags.index(tag) if tag in tags else None

    type_idx = get_idx(_STRUCT_CONN_TYPE)
    p1_chain_idx = get_idx(_STRUCT_CONN_P1_CHAIN)
    p1_seq_idx = get_idx(_STRUCT_CONN_P1_SEQ)
    p1_ins_idx = get_idx(_STRUCT_CONN_P1_INS)
    p1_comp_idx = get_idx(_STRUCT_CONN_P1_COMP)
    p1_atom_idx = get_idx(_STRUCT_CONN_P1_ATOM)
    p2_chain_idx = get_idx(_STRUCT_CONN_P2_CHAIN)
    p2_seq_idx = get_idx(_STRUCT_CONN_P2_SEQ)
    p2_ins_idx = get_idx(_STRUCT_CONN_P2_INS)
    p2_comp_idx = get_idx(_STRUCT_CONN_P2_COMP)
    p2_atom_idx = get_idx(_STRUCT_CONN_P2_ATOM)

    # Must have at least type and chain IDs
    if type_idx is None or p1_chain_idx is None or p2_chain_idx is None:
        return result

    for row_idx in range(loop.length()):
        conn_type = _normalize_cif_value(loop[row_idx, type_idx]).lower()

        if conn_type not in bond_types:
            continue

        # Get partner 1 info
        chain1 = _normalize_cif_value(loop[row_idx, p1_chain_idx])
        seq1 = _normalize_cif_value(loop[row_idx, p1_seq_idx]) if p1_seq_idx is not None else "."
        ins1 = _normalize_cif_value(loop[row_idx, p1_ins_idx]) if p1_ins_idx is not None else "."
        comp1 = _normalize_cif_value(loop[row_idx, p1_comp_idx]) if p1_comp_idx is not None else "."
        atom1 = _normalize_cif_value(loop[row_idx, p1_atom_idx]) if p1_atom_idx is not None else "."

        # Get partner 2 info
        chain2 = _normalize_cif_value(loop[row_idx, p2_chain_idx])
        seq2 = _normalize_cif_value(loop[row_idx, p2_seq_idx]) if p2_seq_idx is not None else "."
        ins2 = _normalize_cif_value(loop[row_idx, p2_ins_idx]) if p2_ins_idx is not None else "."
        comp2 = _normalize_cif_value(loop[row_idx, p2_comp_idx]) if p2_comp_idx is not None else "."
        atom2 = _normalize_cif_value(loop[row_idx, p2_atom_idx]) if p2_atom_idx is not None else "."

        # Skip if missing chain info
        if chain1 == "." or chain2 == ".":
            continue

        result.append(
            AtomBond(
                chain1=chain1,
                seq1=seq1,
                ins1=ins1,
                comp1=comp1,
                atom1=atom1,
                chain2=chain2,
                seq2=seq2,
                ins2=ins2,
                comp2=comp2,
                atom2=atom2,
                conn_type=conn_type,
            )
        )

    return result


def get_intra_chain_bonds(
    block: gemmi.cif.Block,
    bond_types: AbstractSet[str] | None = None,
) -> dict[str, list[tuple[ResKey, ResKey]]]:
    """Get intra-chain bonds from _struct_conn.

    Extracts bonds where both partners are in the same chain.
    This is used to add disulfide bonds and other intra-chain
    covalent bonds to the chain_entity graph.

    Args:
        block: mmCIF data block.
        bond_types: Set of conn_type_id values to include.
            Defaults to {"disulf", "covale"} (disulfide and covalent bonds).

    Returns:
        Dict mapping chain_id to list of (ResKey, ResKey) bonds within that chain.
        Each bond is represented as a tuple of two ResKeys. Bond direction is not
        significant (bonds are treated as undirected).

    Example:
        >>> intra_bonds = get_intra_chain_bonds(block)
        >>> intra_bonds["A"]  # Intra-chain bonds in chain A
        [(("22", "."), ("87", "."))]  # e.g., disulfide between CYS 22 and CYS 87
    """
    if bond_types is None:
        bond_types = frozenset({"disulf", "covale"})
    else:
        bond_types = frozenset(t.lower() for t in bond_types)

    # Get all atom-level bonds matching the bond types
    atom_bonds = get_struct_conn_bonds(block, bond_types)

    # Filter to intra-chain bonds only and convert to ResKey pairs
    result: dict[str, list[tuple[ResKey, ResKey]]] = {}
    seen: dict[str, set[tuple[ResKey, ResKey]]] = {}

    for bond in atom_bonds:
        # Only include bonds within the same chain
        if bond.chain1 != bond.chain2:
            continue

        chain_id = bond.chain1
        res_key1: ResKey = (bond.seq1, bond.ins1)
        res_key2: ResKey = (bond.seq2, bond.ins2)

        # Skip if same residue (self-loop)
        if res_key1 == res_key2:
            continue

        # Normalize bond order for deduplication
        bond_pair = (min(res_key1, res_key2), max(res_key1, res_key2))

        if chain_id not in seen:
            seen[chain_id] = set()
            result[chain_id] = []

        if bond_pair not in seen[chain_id]:
            seen[chain_id].add(bond_pair)
            result[chain_id].append((res_key1, res_key2))

    return result


# Type alias for atom-level bond tuple (res_id, res_name, atom_name) x 2
AtomBondTuple = tuple[str, str, str, str, str, str]


@lru_cache(maxsize=1024)
def _get_chem_type_from_ccd(comp_id: str, ccd_path: str) -> str | None:
    """Get chemical component type from CCD.

    Args:
        comp_id: Component ID (e.g., "ALA", "DAL").
        ccd_path: Path to CCD mirror directory.

    Returns:
        Chemical component type (e.g., "L-PEPTIDE LINKING") or None if not found.
    """
    if not ccd_path:
        return None

    # Try different path structures for CCD
    # Structure 1: {ccd_path}/{first_letter_upper}/{comp_id}/{comp_id}.cif
    cif_file = os.path.join(ccd_path, comp_id[0].upper(), comp_id, f"{comp_id}.cif")

    if not os.path.exists(cif_file):
        # Structure 2: {ccd_path}/{first_letter_lower}/{comp_id}.cif
        cif_file = os.path.join(ccd_path, comp_id[0].lower(), f"{comp_id}.cif")

    if not os.path.exists(cif_file):
        return None

    try:
        doc = gemmi.cif.read(cif_file)
        block = doc[0]
        chem_type = block.find_value("_chem_comp.type")
        if chem_type and chem_type not in ("?", "."):
            # Normalize: remove quotes and convert to uppercase
            return chem_type.strip('"').strip("'").upper()
    except Exception:
        pass

    return None


def _is_polymer_backbone_bond(
    comp1: str,
    atom1: str,
    comp2: str,
    atom2: str,
    ccd_path: str,
) -> bool:
    """Check if a bond is a standard polymer backbone bond.

    Polymer backbone bonds (peptide bonds, phosphodiester bonds) are inferred
    from the polymer structure and should not be included in the inter-level
    bond hash.

    Args:
        comp1: Component ID of first residue.
        atom1: Atom name in first residue.
        comp2: Component ID of second residue.
        atom2: Atom name in second residue.
        ccd_path: Path to CCD mirror directory.

    Returns:
        True if this is a polymer backbone bond, False otherwise.
    """
    if not ccd_path:
        return False

    chem_type1 = _get_chem_type_from_ccd(comp1, ccd_path)
    chem_type2 = _get_chem_type_from_ccd(comp2, ccd_path)

    if not chem_type1 or not chem_type2:
        return False

    # Get polymerization atoms for each residue
    poly_atoms1 = CHEM_TYPE_POLYMERIZATION_ATOMS.get(chem_type1)
    poly_atoms2 = CHEM_TYPE_POLYMERIZATION_ATOMS.get(chem_type2)

    if not poly_atoms1 or not poly_atoms2:
        return False

    # Check if this is a polymer backbone bond:
    # - atom1 is the "exit" atom of comp1 (poly_atoms1[0]) and
    #   atom2 is the "entry" atom of comp2 (poly_atoms1[1])
    # - OR the reverse direction
    exit_atom1, entry_atom1 = poly_atoms1
    exit_atom2, entry_atom2 = poly_atoms2

    # Check forward direction: comp1 exit -> comp2 entry
    if atom1 == exit_atom1 and atom2 == entry_atom2:
        return True

    # Check reverse direction: comp2 exit -> comp1 entry
    if atom2 == exit_atom2 and atom1 == entry_atom1:
        return True

    return False


def get_intra_chain_atom_bonds(
    block: gemmi.cif.Block,
    bond_types: AbstractSet[str] | None = None,
    ccd_path: str = "",
    exclude_polymer_backbone: bool = True,
) -> dict[str, list[AtomBondTuple]]:
    """Get intra-chain atom-level bonds from _struct_conn.

    Returns bonds with full atomic information (res_id, res_name, atom_name)
    for use in inter-level bond hash computation. This matches the format
    used by AtomWorks's generate_inter_level_bond_hash().

    By default, excludes polymer backbone bonds (peptide bonds, phosphodiester
    bonds) as these are inferred from the polymer structure by AtomWorks.

    Args:
        block: mmCIF data block.
        bond_types: Set of conn_type_id values to include.
            Defaults to {"covale"} (covalent bonds only, like AtomWorks).
        ccd_path: Path to CCD mirror directory for chem_type lookup.
            Required if exclude_polymer_backbone is True.
        exclude_polymer_backbone: If True, exclude standard polymer backbone
            bonds (peptide bonds, phosphodiester bonds) from the result.
            Defaults to True for AtomWorks compatibility.

    Returns:
        Dict mapping chain_id to list of bond tuples. Each tuple contains:
        (res_id1, res_name1, atom_name1, res_id2, res_name2, atom_name2)

    Example:
        >>> atom_bonds = get_intra_chain_atom_bonds(block, ccd_path="/path/to/ccd")
        >>> atom_bonds["A"]
        [("82", "GLU", "OE2", "165", "TYR", "CZ"), ...]
    """
    if bond_types is None:
        bond_types = frozenset({"covale"})
    else:
        bond_types = frozenset(t.lower() for t in bond_types)

    # Get all atom-level bonds matching the bond types
    atom_bonds = get_struct_conn_bonds(block, bond_types)

    # Filter to intra-chain bonds and convert to AtomBondTuple
    result: dict[str, list[AtomBondTuple]] = {}
    seen: dict[str, set[AtomBondTuple]] = {}

    for bond in atom_bonds:
        # Only include bonds within the same chain
        if bond.chain1 != bond.chain2:
            continue

        # Skip polymer backbone bonds if requested
        if exclude_polymer_backbone and ccd_path:
            if _is_polymer_backbone_bond(bond.comp1, bond.atom1, bond.comp2, bond.atom2, ccd_path):
                continue

        chain_id = bond.chain1

        # Create atom-level bond tuple
        bond_tuple: AtomBondTuple = (
            bond.seq1,
            bond.comp1,
            bond.atom1,
            bond.seq2,
            bond.comp2,
            bond.atom2,
        )

        # Normalize bond order for deduplication (smaller half first)
        half1 = (bond.seq1, bond.comp1, bond.atom1)
        half2 = (bond.seq2, bond.comp2, bond.atom2)
        if half1 > half2:
            bond_tuple = (
                bond.seq2,
                bond.comp2,
                bond.atom2,
                bond.seq1,
                bond.comp1,
                bond.atom1,
            )

        if chain_id not in seen:
            seen[chain_id] = set()
            result[chain_id] = []

        if bond_tuple not in seen[chain_id]:
            seen[chain_id].add(bond_tuple)
            result[chain_id].append(bond_tuple)

    return result


# 8-tuple: (pn_entity1, res_id1, res_name1, atom_name1, pn_entity2, res_id2, res_name2, atom_name2)
InterChainAtomBondTuple = tuple[str, str, str, str, str, str, str, str]


def get_inter_chain_atom_bonds(
    block: gemmi.cif.Block,
    bond_types: AbstractSet[str] | None = None,
    ccd_path: str = "",
    exclude_polymer_backbone: bool = True,
) -> list[tuple[str, str, AtomBondTuple]]:
    """Get inter-chain atom-level bonds from _struct_conn.

    Returns bonds with full atomic information including chain IDs for use in
    molecule_entity inter-level bond hash computation.

    Args:
        block: mmCIF data block.
        bond_types: Set of conn_type_id values to include.
            Defaults to {"covale"} (covalent bonds only, like AtomWorks).
        ccd_path: Path to CCD mirror directory for chem_type lookup.
        exclude_polymer_backbone: If True, exclude standard polymer backbone
            bonds. Defaults to True for AtomWorks compatibility.

    Returns:
        List of tuples: (chain1, chain2, (res_id1, res_name1, atom_name1,
                                           res_id2, res_name2, atom_name2))

    Example:
        >>> inter_bonds = get_inter_chain_atom_bonds(block)
        >>> inter_bonds
        [("E", "I", ("8", "LYS", "NZ", "201", "8HB", "C1")), ...]
    """
    if bond_types is None:
        bond_types = frozenset({"covale"})
    else:
        bond_types = frozenset(t.lower() for t in bond_types)

    # Get all atom-level bonds matching the bond types
    atom_bonds = get_struct_conn_bonds(block, bond_types)

    # Filter to inter-chain bonds and convert to tuple format
    result: list[tuple[str, str, AtomBondTuple]] = []
    seen: set[tuple[str, str, AtomBondTuple]] = set()

    for bond in atom_bonds:
        # Only include bonds between different chains
        if bond.chain1 == bond.chain2:
            continue

        # Skip polymer backbone bonds if requested
        if exclude_polymer_backbone and ccd_path:
            if _is_polymer_backbone_bond(bond.comp1, bond.atom1, bond.comp2, bond.atom2, ccd_path):
                continue

        # Create atom-level bond tuple (6-tuple)
        bond_tuple: AtomBondTuple = (
            bond.seq1,
            bond.comp1,
            bond.atom1,
            bond.seq2,
            bond.comp2,
            bond.atom2,
        )

        # Normalize bond order for deduplication (smaller chain first)
        if bond.chain1 < bond.chain2:
            entry = (bond.chain1, bond.chain2, bond_tuple)
        else:
            # Swap chain order and bond tuple halves
            entry = (
                bond.chain2,
                bond.chain1,
                (bond.seq2, bond.comp2, bond.atom2, bond.seq1, bond.comp1, bond.atom1),
            )

        if entry not in seen:
            seen.add(entry)
            result.append(entry)

    return result


__all__ = [
    "ResKey",
    "ChainMetadata",
    "AtomBond",
    "AtomBondTuple",
    "InterChainAtomBondTuple",
    "get_canonical_sequences",
    "get_entity_info",
    "get_polymer_info",
    "get_struct_asym_mapping",
    "get_chain_info_complete",
    "get_struct_conn_bonds",
    "get_intra_chain_bonds",
    "get_intra_chain_atom_bonds",
    "get_inter_chain_atom_bonds",
]
