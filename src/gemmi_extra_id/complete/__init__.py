"""Complete mode: AtomWorks-compatible entity assignment.

This module implements AtomWorks-compatible entity ID assignment
using Weisfeiler-Lehman graph hashing with full residue-level
bond detection including inferred polymer bonds.

The complete mode produces entity IDs that exactly match AtomWorks output,
unlike the default loose mode which uses a faster approximation.
"""

from __future__ import annotations

import os
from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi
import numpy as np

from gemmi_extra_id.complete.bonds import get_inferred_polymer_bonds
from gemmi_extra_id.complete.cif_utils import (
    ResKey,
    get_canonical_sequences,
    get_chain_info_complete,
    get_inter_chain_atom_bonds,
    get_intra_chain_bonds,
    get_struct_conn_bonds,
)
from gemmi_extra_id.complete.entities import (
    EntityAssigner,
    find_molecules,
    find_pn_units,
)
from gemmi_extra_id.mmcif import AssignmentResult, ChainInfo

if TYPE_CHECKING:
    from collections.abc import Set as AbstractSet


def _get_atoms_from_atom_site(
    block: gemmi.cif.Block,
    chain_id: str,
) -> tuple[list[str], list[str], list[str], list[int], list[str]]:
    """Get atomic data from _atom_site for a chain.

    Returns:
        Tuple of (atom_names, elements, res_names, res_ids, chem_types).
        chem_types is currently empty (not available in atom_site).
    """
    atom_names: list[str] = []
    elements: list[str] = []
    res_names: list[str] = []
    res_ids: list[int] = []
    chem_types: list[str] = []

    label_col = block.find_loop("_atom_site.label_asym_id")
    if not label_col:
        return atom_names, elements, res_names, res_ids, chem_types

    loop = label_col.get_loop()
    if not loop:
        return atom_names, elements, res_names, res_ids, chem_types

    tags = list(loop.tags)

    def get_idx(tag: str) -> int | None:
        return tags.index(tag) if tag in tags else None

    asym_idx = get_idx("_atom_site.label_asym_id")
    atom_idx = get_idx("_atom_site.label_atom_id")
    elem_idx = get_idx("_atom_site.type_symbol")
    comp_idx = get_idx("_atom_site.label_comp_id")
    seq_idx = get_idx("_atom_site.label_seq_id")

    if asym_idx is None or atom_idx is None:
        return atom_names, elements, res_names, res_ids, chem_types

    # Track current residue number
    current_res_num = 0
    prev_seq_id: str | None = None

    for row_idx in range(loop.length()):
        if loop[row_idx, asym_idx] != chain_id:
            continue

        atom_name = loop[row_idx, atom_idx]
        element = loop[row_idx, elem_idx] if elem_idx is not None else ""
        res_name = loop[row_idx, comp_idx] if comp_idx is not None else ""

        # Handle sequence ID for residue numbering
        seq_id = loop[row_idx, seq_idx] if seq_idx is not None else "."
        if seq_id in (".", "?", ""):
            # For non-polymer (water, ligands), increment on each residue change
            # We'll track by position - each new atom with same res_name is same residue
            # Actually for simplicity, we increment res_num when res_name changes
            # or when seq_id is different
            pass

        # Assign residue ID (use seq_id if numeric, otherwise track internally)
        try:
            res_id = int(seq_id)
        except ValueError:
            # Non-numeric seq_id: use incremental numbering
            if seq_id != prev_seq_id:
                current_res_num += 1
            res_id = current_res_num

        prev_seq_id = seq_id

        atom_names.append(atom_name)
        elements.append(element)
        res_names.append(res_name)
        res_ids.append(res_id)
        chem_types.append("")  # Chem type not in atom_site

    return atom_names, elements, res_names, res_ids, chem_types


def _get_residue_sequence_from_atom_site(
    block: gemmi.cif.Block,
    chain_id: str,
    entity_type: str = "",
) -> list[tuple[int, str]]:
    """Get residue sequence from _atom_site for a chain.

    For water chains (entity_type="water"), uses auth_seq_id to distinguish
    individual water molecules, since all waters have label_seq_id=".".

    Returns:
        List of (res_id, res_name) tuples for each residue.
    """
    result: list[tuple[int, str]] = []
    seen: set[tuple[int, str]] = set()

    label_col = block.find_loop("_atom_site.label_asym_id")
    if not label_col:
        return result

    loop = label_col.get_loop()
    if not loop:
        return result

    tags = list(loop.tags)

    def get_idx(tag: str) -> int | None:
        return tags.index(tag) if tag in tags else None

    asym_idx = get_idx("_atom_site.label_asym_id")
    comp_idx = get_idx("_atom_site.label_comp_id")
    seq_idx = get_idx("_atom_site.label_seq_id")
    auth_seq_idx = get_idx("_atom_site.auth_seq_id")

    if asym_idx is None or comp_idx is None:
        return result

    # For water chains, use auth_seq_id to distinguish molecules
    is_water = entity_type.lower() == "water"
    use_auth_seq = is_water and auth_seq_idx is not None

    # Track current residue number for non-polymer
    current_res_num = 0
    prev_key: tuple[str, str] | None = None

    for row_idx in range(loop.length()):
        if loop[row_idx, asym_idx] != chain_id:
            continue

        res_name = loop[row_idx, comp_idx]

        # Use auth_seq_id for water, label_seq_id otherwise
        if use_auth_seq:
            seq_id = loop[row_idx, auth_seq_idx]
        else:
            seq_id = loop[row_idx, seq_idx] if seq_idx is not None else "."

        # Assign residue ID
        try:
            res_id = int(seq_id)
        except ValueError:
            # Non-numeric seq_id: use incremental numbering
            key = (seq_id, res_name)
            if key != prev_key:
                current_res_num += 1
            res_id = current_res_num
            prev_key = key

        entry = (res_id, res_name)
        if entry not in seen:
            seen.add(entry)
            result.append(entry)

    return result


def _get_residues_from_atom_site(
    block: gemmi.cif.Block,
    chain_id: str,
    entity_type: str = "",
) -> list[tuple[ResKey, str]]:
    """Get residues from _atom_site for a chain.

    Fallback when canonical sequence is not available.
    Returns list of (ResKey, mon_id) tuples.

    For water chains (entity_type="water"), uses auth_seq_id to distinguish
    individual water molecules, since all waters have label_seq_id=".".
    """
    result: list[tuple[ResKey, str]] = []
    seen: set[ResKey] = set()

    # Find the atom_site loop
    label_asym_col = block.find_loop("_atom_site.label_asym_id")
    if not label_asym_col:
        return result

    loop = label_asym_col.get_loop()
    if not loop:
        return result

    tags = list(loop.tags)

    def get_idx(tag: str) -> int | None:
        return tags.index(tag) if tag in tags else None

    asym_idx = get_idx("_atom_site.label_asym_id")
    seq_idx = get_idx("_atom_site.label_seq_id")
    auth_seq_idx = get_idx("_atom_site.auth_seq_id")
    ins_idx = get_idx("_atom_site.pdbx_PDB_ins_code")
    comp_idx = get_idx("_atom_site.label_comp_id")

    if asym_idx is None or seq_idx is None or comp_idx is None:
        return result

    # For water chains, use auth_seq_id to distinguish molecules
    is_water = entity_type.lower() == "water"
    use_auth_seq = is_water and auth_seq_idx is not None

    for row_idx in range(loop.length()):
        asym_id = loop[row_idx, asym_idx]
        if asym_id != chain_id:
            continue

        # For water, use auth_seq_id; otherwise use label_seq_id
        seq_id = loop[row_idx, auth_seq_idx] if use_auth_seq else loop[row_idx, seq_idx]

        ins_code = loop[row_idx, ins_idx] if ins_idx is not None else "."
        comp_id = loop[row_idx, comp_idx]

        # Normalize
        if seq_id in ("?", ".", ""):
            seq_id = "."
        if ins_code in ("?", ".", "", None):
            ins_code = "."

        res_key: ResKey = (seq_id, ins_code)
        if res_key not in seen:
            seen.add(res_key)
            result.append((res_key, comp_id))

    return result


def _get_inter_chain_bonds(
    block: gemmi.cif.Block,
    covalent_types: set[str],
) -> list[tuple[str, str]]:
    """Get inter-chain bonds (chain_id pairs) from struct_conn."""
    atom_bonds = get_struct_conn_bonds(block, covalent_types)
    result: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for bond in atom_bonds:
        if bond.chain1 != bond.chain2:
            # Normalize order for consistency
            pair: tuple[str, str] = (min(bond.chain1, bond.chain2), max(bond.chain1, bond.chain2))
            if pair not in seen:
                seen.add(pair)
                result.append(pair)

    return result


def assign_extended_ids_complete(
    input_path: str | Path,
    output_path: str | Path | None = None,
    covalent_types: AbstractSet[str] | None = None,
    ccd_path: str = "",
) -> AssignmentResult:
    """Assign extended IDs using AtomWorks-compatible algorithm.

    This function implements the complete entity assignment algorithm
    that matches AtomWorks output exactly. It includes:
    - Inferred polymer bonds from coordinates
    - Weisfeiler-Lehman graph hashing at residue level (or atomic level with CCD)
    - Inter-level bond hashing

    Args:
        input_path: Path to input mmCIF file.
        output_path: Path to output mmCIF file. If None, no file is written.
        covalent_types: Set of conn_type_id values to treat as covalent bonds.
            Defaults to {"covale"}. Note: disulf bonds are excluded by default
            to match AtomWorks behavior where disulfide bonds don't affect
            molecule connectivity or entity hashing.
        ccd_path: Path to CCD mirror directory. If provided, uses atomic-level
            hashing with CCD bond information for more accurate entity matching.
            Can also be set via CCD_MIRROR_PATH environment variable.

    Returns:
        AssignmentResult containing ChainInfo for each chain.

    Raises:
        ImportError: If networkx is not installed.
        ValueError: If the CIF file cannot be read.
    """
    # Check for networkx (required for complete mode)
    try:
        import networkx  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "Complete mode requires networkx. Install with: pip install gemmi-extra-id[complete]"
        ) from e

    # Default covalent types - only "covale" to match AtomWorks behavior
    # AtomWorks excludes disulfide bonds from molecule connectivity and entity hashing
    if covalent_types is None:
        covalent_types = {"covale"}
    covalent_types_set = {t.lower() for t in covalent_types}

    # Read CIF file
    input_path = Path(input_path)
    doc = gemmi.cif.read(str(input_path))
    if not doc:
        raise ValueError(f"Could not read CIF file: {input_path}")

    block = doc[0]

    # Get chain metadata
    chain_metadata = get_chain_info_complete(block)

    # Get canonical sequences
    canonical_sequences = get_canonical_sequences(block)

    # Get inter-chain bonds
    inter_chain_bonds = _get_inter_chain_bonds(block, covalent_types_set)

    # Get intra-chain bonds (disulfide, etc.) for chain_entity calculation
    intra_chain_struct_conn = get_intra_chain_bonds(block, covalent_types_set)

    # Determine CCD path (parameter takes precedence over environment variable)
    # Need this early for polymer backbone bond filtering
    effective_ccd_path = ccd_path or os.environ.get("CCD_MIRROR_PATH", "")
    use_atomic_hash = bool(effective_ccd_path)

    # Note: intra-chain atom-level bonds are not used for chain_entity hash
    # (AtomWorks only uses them for pn_unit/molecule levels)

    # Get inter-chain atom-level bonds for molecule_entity inter-level bond hash
    inter_chain_atom_bonds = get_inter_chain_atom_bonds(
        block, covalent_types_set, ccd_path=effective_ccd_path, exclude_polymer_backbone=True
    )

    # Build entity_id -> canonical sequence mapping
    entity_sequences: dict[str, list[tuple[ResKey, str]]] = canonical_sequences

    # Initialize entity assigner
    assigner = EntityAssigner()

    # Build auth_asym_id mapping once (O(m) total, not O(n*m))
    auth_asym_mapping: dict[str, str] = {}
    label_col = block.find_loop("_atom_site.label_asym_id")
    if label_col:
        loop = label_col.get_loop()
        if loop:
            tags = list(loop.tags)
            if "_atom_site.auth_asym_id" in tags:
                label_idx = tags.index("_atom_site.label_asym_id")
                auth_idx = tags.index("_atom_site.auth_asym_id")
                for row_idx in range(loop.length()):
                    label = loop[row_idx, label_idx]
                    if label not in auth_asym_mapping:
                        auth_asym_mapping[label] = loop[row_idx, auth_idx]

    # Process each chain to compute chain_entity
    # Sort chain IDs using numpy's sort order to match AtomWorks behavior
    # np.unique returns sorted values lexicographically ('AA' comes after 'A' but before 'B')
    chain_ids = [str(x) for x in np.sort(list(chain_metadata.keys()))]
    chain_entities: dict[str, int] = {}

    for chain_id in chain_ids:
        meta = chain_metadata[chain_id]
        entity_id = meta.entity_id

        if use_atomic_hash:
            # Use atomic-level hashing with CCD template atoms (normalized)
            # This ensures chains with same sequence have identical hashes

            # For polymer chains, use canonical sequence from _entity_poly_seq
            # For non-polymer, use atom_site
            if entity_id in entity_sequences and entity_sequences[entity_id]:
                # Use canonical sequence (numbered from 1)
                canon_seq = entity_sequences[entity_id]
                residue_sequence = [(i + 1, mon_id) for i, (_, mon_id) in enumerate(canon_seq)]
            else:
                # Fallback to atom_site for non-polymer chains
                residue_sequence = _get_residue_sequence_from_atom_site(
                    block, chain_id, meta.entity_type
                )

            # Note: Inter-level bonds (from struct_conn) are NOT included in chain_entity hash.
            # AtomWorks only uses them for pn_unit_entity and molecule_entity levels.
            # Chain entity is determined purely by residue sequence and atom composition.

            if not residue_sequence:
                # Empty chain
                chain_entities[chain_id] = assigner.assign_chain_entity_normalized(
                    [], None, effective_ccd_path, inter_level_bonds=None
                )
            else:
                chain_entities[chain_id] = assigner.assign_chain_entity_normalized(
                    residue_sequence=residue_sequence,
                    chem_types=None,  # Could extract from entity_poly if needed
                    ccd_path=effective_ccd_path,
                    inter_level_bonds=None,
                )
        else:
            # Use residue-level hashing (original behavior)
            # Get residues for this chain
            if entity_id in entity_sequences and entity_sequences[entity_id]:
                residues_with_names = entity_sequences[entity_id]
            else:
                # Fallback to atom_site (pass entity_type for water handling)
                residues_with_names = _get_residues_from_atom_site(
                    block, chain_id, meta.entity_type
                )

            if not residues_with_names:
                # Empty chain
                chain_entities[chain_id] = assigner.assign_chain_entity([], {}, [])
                continue

            # Build residue list and names mapping
            residues = [r[0] for r in residues_with_names]
            residue_names = {r[0]: r[1] for r in residues_with_names}

            # Infer polymer bonds for this chain
            polymer_bonds = get_inferred_polymer_bonds(residues)

            # Get intra-chain struct_conn bonds (disulfide, etc.) for this chain
            struct_conn_bonds = intra_chain_struct_conn.get(chain_id, [])

            # Combine polymer bonds and struct_conn bonds
            intra_bonds = polymer_bonds + struct_conn_bonds

            # Compute chain entity
            chain_entities[chain_id] = assigner.assign_chain_entity(
                residues, residue_names, intra_bonds
            )

    # Find PN units
    is_polymer = {cid: chain_metadata[cid].is_polymer for cid in chain_ids}
    pn_unit_mapping = find_pn_units(chain_ids, is_polymer, inter_chain_bonds)

    # Get unique PN unit IDs in order
    seen_pn_units: set[str] = set()
    pn_unit_ids: list[str] = []
    for chain_id in chain_ids:
        pn_id = pn_unit_mapping[chain_id]
        if pn_id not in seen_pn_units:
            seen_pn_units.add(pn_id)
            pn_unit_ids.append(pn_id)

    # Compute pn_unit_entity for each PN unit
    pn_unit_entities: dict[str, int] = {}
    for pn_unit_id in pn_unit_ids:
        pn_unit_entities[pn_unit_id] = assigner.assign_pn_unit_entity(
            pn_unit_id, chain_entities, inter_chain_bonds
        )

    # Find molecules (group connected PN units)
    # Inter-PN-unit bonds are inter-chain bonds where chains are in different PN units
    inter_pn_unit_bonds: list[tuple[str, str]] = []
    seen_pn_bonds: set[tuple[str, str]] = set()
    for chain_a, chain_b in inter_chain_bonds:
        pn_a = pn_unit_mapping.get(chain_a)
        pn_b = pn_unit_mapping.get(chain_b)
        if pn_a and pn_b and pn_a != pn_b:
            pair: tuple[str, str] = (min(pn_a, pn_b), max(pn_a, pn_b))
            if pair not in seen_pn_bonds:
                seen_pn_bonds.add(pair)
                inter_pn_unit_bonds.append(pair)

    molecule_mapping = find_molecules(pn_unit_ids, inter_pn_unit_bonds)

    # Get unique molecule IDs and their PN units
    molecule_pn_units: dict[int, list[str]] = {}
    for pn_unit_id in pn_unit_ids:
        mol_id = molecule_mapping[pn_unit_id]
        if mol_id not in molecule_pn_units:
            molecule_pn_units[mol_id] = []
        molecule_pn_units[mol_id].append(pn_unit_id)

    # Compute molecule_entity for each molecule
    molecule_entities: dict[int, int] = {}
    for mol_id, pn_units in molecule_pn_units.items():
        # Get inter-pn-unit bonds within this molecule
        mol_pn_set = set(pn_units)
        mol_inter_bonds = [
            (a, b) for a, b in inter_pn_unit_bonds if a in mol_pn_set and b in mol_pn_set
        ]

        # Build inter-level bonds for molecule entity (8-tuple format)
        # (pn_unit_entity1, res_id1, res_name1, atom_name1,
        #  pn_unit_entity2, res_id2, res_name2, atom_name2)
        mol_inter_level_bonds: list[tuple[str, str, str, str, str, str, str, str]] = []
        for chain1, chain2, bond_tuple in inter_chain_atom_bonds:
            pn1 = pn_unit_mapping.get(chain1)
            pn2 = pn_unit_mapping.get(chain2)
            # Only include bonds between different PN units within this molecule
            if pn1 and pn2 and pn1 != pn2 and pn1 in mol_pn_set and pn2 in mol_pn_set:
                # bond_tuple is (res_id1, res_name1, atom_name1, res_id2, res_name2, atom_name2)
                res_id1, res_name1, atom_name1, res_id2, res_name2, atom_name2 = bond_tuple
                # Create 8-tuple with pn_unit_entity values
                inter_bond: tuple[str, str, str, str, str, str, str, str] = (
                    str(pn_unit_entities[pn1]),
                    str(res_id1),
                    res_name1,
                    atom_name1,
                    str(pn_unit_entities[pn2]),
                    str(res_id2),
                    res_name2,
                    atom_name2,
                )
                mol_inter_level_bonds.append(inter_bond)

        molecule_entities[mol_id] = assigner.assign_molecule_entity(
            pn_units, pn_unit_entities, mol_inter_bonds, mol_inter_level_bonds or None
        )

    # Build ChainInfo for each chain
    chain_info: OrderedDict[str, ChainInfo] = OrderedDict()

    for chain_id in chain_ids:
        meta = chain_metadata[chain_id]
        pn_unit_id = pn_unit_mapping[chain_id]
        mol_id = molecule_mapping[pn_unit_id]

        # Get auth_asym_id from pre-built mapping (fallback to chain_id)
        auth_asym_id = auth_asym_mapping.get(chain_id, chain_id)

        chain_info[chain_id] = ChainInfo(
            label_asym_id=chain_id,
            auth_asym_id=auth_asym_id,
            entity_id=meta.entity_id,
            entity_type=meta.entity_type,
            molecule_id=mol_id,
            pn_unit_id=pn_unit_id,
            chain_entity=str(chain_entities[chain_id]),
            pn_unit_entity=str(pn_unit_entities[pn_unit_id]),
            molecule_entity=str(molecule_entities[mol_id]),
        )

    result = AssignmentResult(chain_info=chain_info)

    # Write output if requested
    if output_path is not None:
        _write_complete_output(block, chain_info, output_path)

    return result


def _write_complete_output(
    block: gemmi.cif.Block,
    chain_info: OrderedDict[str, ChainInfo],
    output_path: str | Path,
) -> None:
    """Write CIF output with assigned IDs."""
    # Add molecule_id to atom_site
    label_col = block.find_loop("_atom_site.label_asym_id")
    if not label_col:
        return

    loop = label_col.get_loop()
    if not loop:
        return

    tags = list(loop.tags)
    label_idx = tags.index("_atom_site.label_asym_id")

    # Check if molecule_id column exists, add if not
    mol_id_tag = "_atom_site.molecule_id"
    if mol_id_tag not in tags:
        loop.add_tag(mol_id_tag)
        tags = list(loop.tags)

    mol_idx = tags.index(mol_id_tag)

    # Set molecule_id for each atom
    for row_idx in range(loop.length()):
        chain_id = loop[row_idx, label_idx]
        if chain_id in chain_info:
            loop[row_idx, mol_idx] = str(chain_info[chain_id].molecule_id)
        else:
            loop[row_idx, mol_idx] = "."

    # Write output
    output_path = Path(output_path)
    doc = gemmi.cif.Document()
    doc.add_copied_block(block)
    doc.write_file(str(output_path))


__all__ = ["assign_extended_ids_complete"]
