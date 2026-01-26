"""Atomic-level WL hash computation using CCD bond information.

This module provides functions to compute WL (Weisfeiler-Lehman) hashes
at the atomic level, using bond information from the CCD.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import networkx as nx

from .ccd import get_residue_bonds, load_ccd_component

if TYPE_CHECKING:
    from collections.abc import Sequence


# Leaving atoms that should be removed for non-terminal residues
# OXT is the terminal carboxyl oxygen
LEAVING_ATOMS = frozenset({"OXT", "HXT"})

logger = logging.getLogger(__name__)

# Polymer bond atoms for amino acids (C of residue i to N of residue i+1)
AA_POLYMER_BOND = ("C", "N")
# Polymer bond atoms for nucleic acids (O3' of residue i to P of residue i+1)
NA_POLYMER_BOND = ("O3'", "P")

# Chemical types that use AA polymer bonds
AA_CHEM_TYPES = frozenset(
    {
        "L-PEPTIDE LINKING",
        "D-PEPTIDE LINKING",
        "PEPTIDE LINKING",
        "L-PEPTIDE NH3 AMINO TERMINUS",
        "L-PEPTIDE COOH CARBOXY TERMINUS",
        "D-PEPTIDE NH3 AMINO TERMINUS",
        "D-PEPTIDE COOH CARBOXY TERMINUS",
    }
)

# Chemical types that use NA polymer bonds
NA_CHEM_TYPES = frozenset(
    {
        "RNA LINKING",
        "DNA LINKING",
        "RNA OH 3 PRIME TERMINUS",
        "RNA OH 5 PRIME TERMINUS",
        "DNA OH 3 PRIME TERMINUS",
        "DNA OH 5 PRIME TERMINUS",
    }
)


def _sum_string_arrays(*arrays: Sequence[str]) -> list[str]:
    """Concatenate string arrays element-wise."""
    return ["|".join(strs) for strs in zip(*arrays, strict=True)]


def build_atom_graph(
    atom_names: Sequence[str],
    elements: Sequence[str],
    res_names: Sequence[str],
    res_ids: Sequence[int],
    chem_types: Sequence[str] | None = None,
    ccd_path: str = "",
    add_polymer_bonds: bool = True,
) -> nx.Graph:
    """Build an atomic-level bond graph for a chain.

    Args:
        atom_names: Atom names (e.g., ["N", "CA", "C", "O", ...]).
        elements: Element symbols (e.g., ["N", "C", "C", "O", ...]).
        res_names: Residue names per atom (e.g., ["ALA", "ALA", ...]).
        res_ids: Residue IDs per atom.
        chem_types: Chemical types per residue (optional, for polymer bond detection).
        ccd_path: Path to CCD mirror directory.
        add_polymer_bonds: Whether to add polymer bonds between residues.

    Returns:
        NetworkX graph with nodes for atoms and edges for bonds.
        Node attributes: "node_data" = element|atom_name
        Edge attributes: "bond_type" = bond order (1, 2, 3, 5)
    """
    n_atoms = len(atom_names)

    # Create graph with all nodes
    graph = nx.Graph()
    graph.add_nodes_from(range(n_atoms))

    # Set node attributes (element|atom_name)
    node_data = _sum_string_arrays(elements, atom_names)
    nx.set_node_attributes(graph, {i: node_data[i] for i in range(n_atoms)}, "node_data")

    # Group atoms by residue
    res_starts: dict[int, int] = {}  # res_id -> first atom index
    res_ends: dict[int, int] = {}  # res_id -> last atom index + 1

    for i, rid in enumerate(res_ids):
        if rid not in res_starts:
            res_starts[rid] = i
        res_ends[rid] = i + 1

    # Add intra-residue bonds from CCD
    for rid in res_starts:
        start = res_starts[rid]
        end = res_ends[rid]

        res_atom_names = list(atom_names[start:end])
        res_name = res_names[start]

        bonds = get_residue_bonds(res_name, res_atom_names, ccd_path)
        for local_a1, local_a2, order in bonds:
            global_a1 = start + local_a1
            global_a2 = start + local_a2
            graph.add_edge(global_a1, global_a2, bond_type=order)

    # Add polymer bonds between consecutive residues
    if add_polymer_bonds:
        sorted_rids = sorted(res_starts.keys())
        for i in range(len(sorted_rids) - 1):
            rid1 = sorted_rids[i]
            rid2 = sorted_rids[i + 1]

            # Skip if not consecutive
            if rid2 - rid1 > 1:
                continue

            # Determine polymer bond atoms based on chemical type
            bond_atoms = None
            if chem_types is not None:
                ctype = chem_types[res_starts[rid1]]
                if ctype in AA_CHEM_TYPES:
                    bond_atoms = AA_POLYMER_BOND
                elif ctype in NA_CHEM_TYPES:
                    bond_atoms = NA_POLYMER_BOND
            else:
                # Default to AA polymer bond
                bond_atoms = AA_POLYMER_BOND

            if bond_atoms is None:
                continue

            # Find atom indices
            start1, end1 = res_starts[rid1], res_ends[rid1]
            start2, end2 = res_starts[rid2], res_ends[rid2]

            atom1_idx = None
            atom2_idx = None

            for j in range(start1, end1):
                if atom_names[j] == bond_atoms[0]:
                    atom1_idx = j
                    break

            for j in range(start2, end2):
                if atom_names[j] == bond_atoms[1]:
                    atom2_idx = j
                    break

            if atom1_idx is not None and atom2_idx is not None:
                graph.add_edge(atom1_idx, atom2_idx, bond_type=1)  # Single bond

    return graph


def hash_atom_graph(
    graph: nx.Graph,
    iterations: int = 3,
    digest_size: int = 16,
) -> str:
    """Compute WL hash for an atomic bond graph.

    Args:
        graph: NetworkX graph with node_data and bond_type attributes.
        iterations: Number of WL iterations.
        digest_size: Digest size for WL hash.

    Returns:
        Hash string including WL hash, node count, and node label distribution.
    """
    # Compute WL hash
    wl_hash = nx.algorithms.graph_hashing.weisfeiler_lehman_graph_hash(
        graph,
        node_attr="node_data",
        edge_attr="bond_type",
        iterations=iterations,
        digest_size=digest_size,
    )

    # Add node count
    wl_hash += f"_{len(graph.nodes)}"

    # Add node label distribution
    node_attrs = nx.get_node_attributes(graph, "node_data")
    label_counts: dict[str, int] = {}
    for label in node_attrs.values():
        label_counts[label] = label_counts.get(label, 0) + 1

    # Sort labels for deterministic output
    sorted_labels = sorted(label_counts.items())
    label_str = ",".join(f"{label}:{count}" for label, count in sorted_labels)
    wl_hash += f"_{label_str}"

    return wl_hash


def build_normalized_atom_graph(
    residue_sequence: Sequence[tuple[int, str]],
    chem_types: Sequence[str] | None = None,
    ccd_path: str = "",
    remove_hydrogens: bool = True,
    add_polymer_bonds: bool = True,
) -> nx.Graph:
    """Build an atomic-level bond graph using CCD template atoms.

    Unlike build_atom_graph which uses actual atoms from the structure,
    this function uses CCD template atoms. This ensures that chains with
    the same residue sequence have identical atom graphs, regardless of
    which atoms are actually present in the structure.

    Args:
        residue_sequence: List of (res_id, res_name) tuples for each residue.
        chem_types: Chemical types per residue (optional, for polymer bond detection).
        ccd_path: Path to CCD mirror directory.
        remove_hydrogens: Whether to remove hydrogen atoms from templates.
        add_polymer_bonds: Whether to add polymer bonds between residues.

    Returns:
        NetworkX graph with nodes for template atoms and edges for bonds.
    """
    if not residue_sequence:
        return nx.Graph()

    # Build atom lists from CCD templates
    all_atom_names: list[str] = []
    all_elements: list[str] = []
    all_res_ids: list[int] = []
    all_res_names: list[str] = []

    # Track residue boundaries
    res_starts: dict[int, int] = {}
    res_ends: dict[int, int] = {}

    n_residues = len(residue_sequence)
    atom_offset = 0

    for res_idx, (res_id, res_name) in enumerate(residue_sequence):
        comp = load_ccd_component(res_name, ccd_path)

        if comp is None:
            # No CCD data - skip this residue or use placeholder
            logger.debug(f"No CCD template for {res_name}, skipping")
            continue

        # Determine if this is a terminal residue (for leaving atom handling)
        is_last = res_idx == n_residues - 1

        # Filter atoms
        res_atoms: list[str] = []
        res_elements: list[str] = []

        for atom_name, element in zip(comp.atoms, comp.elements, strict=True):
            # Skip hydrogens if requested
            if remove_hydrogens and element == "H":
                continue
            # Skip leaving atoms for non-terminal residues
            if not is_last and atom_name in LEAVING_ATOMS:
                continue
            res_atoms.append(atom_name)
            res_elements.append(element)

        if not res_atoms:
            continue

        res_starts[res_id] = atom_offset
        res_ends[res_id] = atom_offset + len(res_atoms)

        all_atom_names.extend(res_atoms)
        all_elements.extend(res_elements)
        all_res_ids.extend([res_id] * len(res_atoms))
        all_res_names.extend([res_name] * len(res_atoms))

        atom_offset += len(res_atoms)

    if not all_atom_names:
        return nx.Graph()

    # Create graph with all nodes
    n_atoms = len(all_atom_names)
    graph = nx.Graph()
    graph.add_nodes_from(range(n_atoms))

    # Set node attributes (element|atom_name)
    node_data = _sum_string_arrays(all_elements, all_atom_names)
    nx.set_node_attributes(graph, {i: node_data[i] for i in range(n_atoms)}, "node_data")

    # Add intra-residue bonds from CCD
    for res_id in res_starts:
        start = res_starts[res_id]
        end = res_ends[res_id]

        res_atom_names = all_atom_names[start:end]
        res_name = all_res_names[start]

        bonds = get_residue_bonds(res_name, res_atom_names, ccd_path)
        for local_a1, local_a2, order in bonds:
            global_a1 = start + local_a1
            global_a2 = start + local_a2
            graph.add_edge(global_a1, global_a2, bond_type=order)

    # Add polymer bonds between consecutive residues
    if add_polymer_bonds:
        sorted_rids = sorted(res_starts.keys())
        res_id_to_chem_type: dict[int, str] = {}
        if chem_types:
            for i, (res_id, _) in enumerate(residue_sequence):
                if i < len(chem_types):
                    res_id_to_chem_type[res_id] = chem_types[i]

        for i in range(len(sorted_rids) - 1):
            rid1 = sorted_rids[i]
            rid2 = sorted_rids[i + 1]

            # Skip if not consecutive
            if rid2 - rid1 > 1:
                continue

            # Determine polymer bond atoms based on chemical type
            bond_atoms = None
            ctype = res_id_to_chem_type.get(rid1, "")
            if ctype in AA_CHEM_TYPES:
                bond_atoms = AA_POLYMER_BOND
            elif ctype in NA_CHEM_TYPES:
                bond_atoms = NA_POLYMER_BOND
            else:
                # Default to AA polymer bond
                bond_atoms = AA_POLYMER_BOND

            if bond_atoms is None:
                continue

            # Find atom indices
            start1, end1 = res_starts[rid1], res_ends[rid1]
            start2, end2 = res_starts[rid2], res_ends[rid2]

            atom1_idx = None
            atom2_idx = None

            for j in range(start1, end1):
                if all_atom_names[j] == bond_atoms[0]:
                    atom1_idx = j
                    break

            for j in range(start2, end2):
                if all_atom_names[j] == bond_atoms[1]:
                    atom2_idx = j
                    break

            if atom1_idx is not None and atom2_idx is not None:
                graph.add_edge(atom1_idx, atom2_idx, bond_type=1)

    return graph


def compute_chain_atom_hash(
    atom_names: Sequence[str],
    elements: Sequence[str],
    res_names: Sequence[str],
    res_ids: Sequence[int],
    chem_types: Sequence[str] | None = None,
    ccd_path: str = "",
) -> str:
    """Compute atomic-level WL hash for a chain using actual atoms.

    This is a convenience function that builds the graph and computes the hash.

    Args:
        atom_names: Atom names for all atoms in the chain.
        elements: Element symbols for all atoms.
        res_names: Residue names per atom.
        res_ids: Residue IDs per atom.
        chem_types: Chemical types per residue (optional).
        ccd_path: Path to CCD mirror directory.

    Returns:
        WL hash string.
    """
    graph = build_atom_graph(
        atom_names=atom_names,
        elements=elements,
        res_names=res_names,
        res_ids=res_ids,
        chem_types=chem_types,
        ccd_path=ccd_path,
    )
    return hash_atom_graph(graph)


def compute_chain_atom_hash_normalized(
    residue_sequence: Sequence[tuple[int, str]],
    chem_types: Sequence[str] | None = None,
    ccd_path: str = "",
    remove_hydrogens: bool = True,
) -> str:
    """Compute atomic-level WL hash for a chain using CCD template atoms.

    This function uses CCD templates to normalize atom lists, ensuring that
    chains with the same residue sequence have identical hashes regardless
    of which atoms are actually present in the structure.

    Args:
        residue_sequence: List of (res_id, res_name) tuples for each residue.
        chem_types: Chemical types per residue (optional).
        ccd_path: Path to CCD mirror directory.
        remove_hydrogens: Whether to remove hydrogen atoms from templates.

    Returns:
        WL hash string.
    """
    graph = build_normalized_atom_graph(
        residue_sequence=residue_sequence,
        chem_types=chem_types,
        ccd_path=ccd_path,
        remove_hydrogens=remove_hydrogens,
    )
    return hash_atom_graph(graph)
