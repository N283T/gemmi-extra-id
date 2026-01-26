"""CIF data reading utilities for complete mode.

This module provides functions to read CIF data required for
AtomWorks-compatible entity assignment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import gemmi

if TYPE_CHECKING:
    from collections.abc import Set as AbstractSet

# Type alias for residue key (seq_id, ins_code)
ResKey = tuple[str, str]


@dataclass
class ChainInfo:
    """Enhanced chain information for complete mode."""

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
) -> dict[str, ChainInfo]:
    """Get complete chain information from CIF.

    Combines information from _entity, _entity_poly, and _struct_asym
    to build comprehensive chain information.

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping chain_id (label_asym_id) to ChainInfo.
    """
    # Get entity type info
    entity_types = get_entity_info(block)

    # Get polymer info
    polymer_info = get_polymer_info(block)

    # Get chain -> entity mapping
    asym_mapping = get_struct_asym_mapping(block)

    # Build ChainInfo for each chain
    result: dict[str, ChainInfo] = {}

    for chain_id, entity_id in asym_mapping.items():
        entity_type = entity_types.get(entity_id, "unknown")
        is_polymer = entity_type == "polymer"

        polymer_type: str | None = None
        if entity_id in polymer_info:
            polymer_type = polymer_info[entity_id][0]

        result[chain_id] = ChainInfo(
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


__all__ = [
    "ResKey",
    "ChainInfo",
    "AtomBond",
    "get_canonical_sequences",
    "get_entity_info",
    "get_polymer_info",
    "get_struct_asym_mapping",
    "get_chain_info_complete",
    "get_struct_conn_bonds",
]
