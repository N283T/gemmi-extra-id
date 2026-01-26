"""Complete mode: AtomWorks-compatible entity assignment.

This module implements AtomWorks-compatible entity ID assignment
using Weisfeiler-Lehman graph hashing with full residue-level
bond detection including inferred polymer bonds.

The complete mode produces entity IDs that exactly match AtomWorks output,
unlike the default loose mode which uses a faster approximation.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi

from gemmi_extra_id.complete.bonds import get_inferred_polymer_bonds
from gemmi_extra_id.complete.cif_utils import (
    ResKey,
    get_canonical_sequences,
    get_chain_info_complete,
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


def _get_residues_from_atom_site(
    block: gemmi.cif.Block,
    chain_id: str,
) -> list[tuple[ResKey, str]]:
    """Get residues from _atom_site for a chain.

    Fallback when canonical sequence is not available.
    Returns list of (ResKey, mon_id) tuples.
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
    ins_idx = get_idx("_atom_site.pdbx_PDB_ins_code")
    comp_idx = get_idx("_atom_site.label_comp_id")

    if asym_idx is None or seq_idx is None or comp_idx is None:
        return result

    for row_idx in range(loop.length()):
        asym_id = loop[row_idx, asym_idx]
        if asym_id != chain_id:
            continue

        seq_id = loop[row_idx, seq_idx]
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

    for bond in atom_bonds:
        if bond.chain1 != bond.chain2:
            # Normalize order for consistency
            pair = tuple(sorted([bond.chain1, bond.chain2]))
            if pair not in result:
                result.append(pair)

    return result


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

    # Default covalent types
    if covalent_types is None:
        covalent_types = {"covale", "disulf"}
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

    # Build entity_id -> canonical sequence mapping
    entity_sequences: dict[str, list[tuple[ResKey, str]]] = canonical_sequences

    # Initialize entity assigner
    assigner = EntityAssigner()

    # Process each chain to compute chain_entity
    chain_ids = list(chain_metadata.keys())
    chain_entities: dict[str, int] = {}
    chain_residues: dict[str, list[tuple[ResKey, str]]] = {}

    for chain_id in chain_ids:
        meta = chain_metadata[chain_id]
        entity_id = meta.entity_id

        # Get residues for this chain
        if entity_id in entity_sequences and entity_sequences[entity_id]:
            residues_with_names = entity_sequences[entity_id]
        else:
            # Fallback to atom_site
            residues_with_names = _get_residues_from_atom_site(block, chain_id)

        chain_residues[chain_id] = residues_with_names

        if not residues_with_names:
            # Empty chain
            chain_entities[chain_id] = assigner.assign_chain_entity([], {}, [])
            continue

        # Build residue list and names mapping
        residues = [r[0] for r in residues_with_names]
        residue_names = {r[0]: r[1] for r in residues_with_names}

        # Infer polymer bonds for this chain
        intra_bonds = get_inferred_polymer_bonds(residues)

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
    for chain_a, chain_b in inter_chain_bonds:
        pn_a = pn_unit_mapping.get(chain_a)
        pn_b = pn_unit_mapping.get(chain_b)
        if pn_a and pn_b and pn_a != pn_b:
            pair = tuple(sorted([pn_a, pn_b]))
            if pair not in inter_pn_unit_bonds:
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

        molecule_entities[mol_id] = assigner.assign_molecule_entity(
            pn_units, pn_unit_entities, mol_inter_bonds
        )

    # Build ChainInfo for each chain
    chain_info: OrderedDict[str, ChainInfo] = OrderedDict()

    for chain_id in chain_ids:
        meta = chain_metadata[chain_id]
        pn_unit_id = pn_unit_mapping[chain_id]
        mol_id = molecule_mapping[pn_unit_id]

        # Get auth_asym_id from atom_site (fallback to chain_id)
        auth_asym_id = chain_id  # Default
        label_col = block.find_loop("_atom_site.label_asym_id")
        if label_col:
            loop = label_col.get_loop()
            if loop:
                tags = list(loop.tags)
                if "_atom_site.auth_asym_id" in tags:
                    label_idx = tags.index("_atom_site.label_asym_id")
                    auth_idx = tags.index("_atom_site.auth_asym_id")
                    for row_idx in range(loop.length()):
                        if loop[row_idx, label_idx] == chain_id:
                            auth_asym_id = loop[row_idx, auth_idx]
                            break

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
