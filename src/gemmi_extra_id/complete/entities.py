"""Entity computation for complete mode.

This module provides functions to compute entity IDs at different levels
(chain, pn_unit, molecule) using AtomWorks-compatible algorithms.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from gemmi_extra_id.complete.cif_utils import ResKey
from gemmi_extra_id.complete.hash import (
    HashToIntMapper,
    generate_inter_level_bond_hash,
    hash_graph,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


def _check_networkx() -> None:
    """Check if networkx is available."""
    if not HAS_NETWORKX:
        raise ImportError(
            "Complete mode requires networkx. Install with: pip install gemmi-extra-id[complete]"
        )


def compute_chain_entity_hash(
    residues: Sequence[ResKey],
    residue_names: Mapping[ResKey, str],
    intra_bonds: Sequence[tuple[ResKey, ResKey]],
) -> str:
    """Compute chain entity hash using WL graph hashing.

    Builds a graph where:
    - Nodes = residues (indexed 0 to n-1)
    - Edges = intra-chain bonds (polymer bonds)
    - Node attribute = residue name (mon_id)

    Args:
        residues: Ordered list of residue keys in the chain.
        residue_names: Mapping from ResKey to residue name (mon_id).
        intra_bonds: List of (ResKey, ResKey) tuples representing bonds.

    Returns:
        WL hash string for the chain entity.

    Example:
        >>> residues = [("1", "."), ("2", "."), ("3", ".")]
        >>> names = {("1", "."): "ALA", ("2", "."): "GLY", ("3", "."): "ALA"}
        >>> bonds = [(("1", "."), ("2", ".")), (("2", "."), ("3", "."))]
        >>> compute_chain_entity_hash(residues, names, bonds)
        '<wl_hash>_3_ALA:2,GLY:1'
    """
    _check_networkx()

    if not residues:
        return ""

    # Create residue index mapping
    res_to_idx: dict[ResKey, int] = {r: i for i, r in enumerate(residues)}

    # Build graph
    graph = nx.Graph()
    graph.add_nodes_from(range(len(residues)))

    # Add edges from bonds
    for res_a, res_b in intra_bonds:
        if res_a in res_to_idx and res_b in res_to_idx:
            graph.add_edge(res_to_idx[res_a], res_to_idx[res_b])

    # Set node attributes (residue names)
    node_attrs = {i: residue_names.get(r, "UNK") for i, r in enumerate(residues)}
    nx.set_node_attributes(graph, node_attrs, "node")

    return hash_graph(graph, node_attr="node")


def find_pn_units(
    chain_ids: Sequence[str],
    is_polymer: Mapping[str, bool],
    inter_chain_bonds: Sequence[tuple[str, str]],
) -> dict[str, str]:
    """Find PN (Polymer/Non-polymer) units.

    Groups non-polymer chains that are covalently bonded together.
    Each polymer chain is its own PN unit.

    Args:
        chain_ids: Ordered list of chain IDs.
        is_polymer: Mapping from chain_id to is_polymer flag.
        inter_chain_bonds: List of (chain_id, chain_id) tuples for inter-chain bonds.

    Returns:
        Mapping from chain_id to pn_unit_id.
        pn_unit_id is a comma-separated string of chain_ids in the unit.

    Example:
        >>> chain_ids = ["A", "B", "C", "D"]
        >>> is_polymer = {"A": True, "B": True, "C": False, "D": False}
        >>> bonds = [("C", "D")]  # C and D are bonded
        >>> find_pn_units(chain_ids, is_polymer, bonds)
        {"A": "A", "B": "B", "C": "C,D", "D": "C,D"}
    """
    _check_networkx()

    # Build graph of non-polymer chains connected by bonds
    non_polymer_graph = nx.Graph()
    non_polymer_chains = [c for c in chain_ids if not is_polymer.get(c, False)]
    non_polymer_graph.add_nodes_from(non_polymer_chains)

    # Add edges for bonds between non-polymer chains
    for chain_a, chain_b in inter_chain_bonds:
        if chain_a in non_polymer_chains and chain_b in non_polymer_chains:
            non_polymer_graph.add_edge(chain_a, chain_b)

    # Find connected components for non-polymer chains
    result: dict[str, str] = {}

    # Polymer chains are their own PN unit
    for chain_id in chain_ids:
        if is_polymer.get(chain_id, False):
            result[chain_id] = chain_id

    # Non-polymer chains grouped by connected components
    for component in nx.connected_components(non_polymer_graph):
        sorted_chains = sorted(component)
        pn_unit_id = ",".join(sorted_chains)
        for chain_id in sorted_chains:
            result[chain_id] = pn_unit_id

    return result


def compute_pn_unit_entity_hash(
    pn_unit_id: str,
    chain_entities: Mapping[str, int],
    inter_chain_bonds: Sequence[tuple[str, str]],
    inter_level_bonds: Sequence[tuple[str, str, str, str, str, str, str, str]]
    | None = None,
) -> str:
    """Compute PN unit entity hash using WL graph hashing.

    Builds a graph where:
    - Nodes = chains in the PN unit
    - Edges = inter-chain bonds within the PN unit
    - Node attribute = chain_entity

    Args:
        pn_unit_id: The PN unit ID (comma-separated chain IDs).
        chain_entities: Mapping from chain_id to chain_entity.
        inter_chain_bonds: List of (chain_id, chain_id) tuples for all inter-chain bonds.
        inter_level_bonds: Optional list of atomic-level bond info for inter-level hash.

    Returns:
        WL hash string for the PN unit entity.
    """
    _check_networkx()

    chain_ids = pn_unit_id.split(",")
    if not chain_ids:
        return ""

    chain_set = set(chain_ids)
    chain_to_idx = {c: i for i, c in enumerate(chain_ids)}

    # Build graph
    graph = nx.Graph()
    graph.add_nodes_from(range(len(chain_ids)))

    # Add edges for inter-chain bonds within this PN unit
    for chain_a, chain_b in inter_chain_bonds:
        if chain_a in chain_set and chain_b in chain_set and chain_a != chain_b:
            graph.add_edge(chain_to_idx[chain_a], chain_to_idx[chain_b])

    # Set node attributes (chain entities)
    node_attrs = {i: str(chain_entities.get(c, 0)) for i, c in enumerate(chain_ids)}
    nx.set_node_attributes(graph, node_attrs, "node")

    wl_hash = hash_graph(graph, node_attr="node")

    # Add inter-level bond hash if provided
    if inter_level_bonds:
        inter_hash = generate_inter_level_bond_hash(inter_level_bonds)
        if inter_hash:
            wl_hash += "_" + inter_hash

    return wl_hash


def find_molecules(
    pn_unit_ids: Sequence[str],
    inter_pn_unit_bonds: Sequence[tuple[str, str]],
) -> dict[str, int]:
    """Find molecules by grouping connected PN units.

    Args:
        pn_unit_ids: Ordered list of unique PN unit IDs.
        inter_pn_unit_bonds: List of (pn_unit_id, pn_unit_id) tuples for bonds.

    Returns:
        Mapping from pn_unit_id to molecule_id (0-indexed).

    Example:
        >>> pn_unit_ids = ["A", "B", "C,D"]
        >>> bonds = [("A", "C,D")]  # A and C,D are bonded
        >>> find_molecules(pn_unit_ids, bonds)
        {"A": 0, "B": 1, "C,D": 0}
    """
    _check_networkx()

    # Build graph of PN units
    graph = nx.Graph()
    graph.add_nodes_from(pn_unit_ids)
    graph.add_edges_from(inter_pn_unit_bonds)

    # Find connected components
    result: dict[str, int] = {}
    for mol_id, component in enumerate(nx.connected_components(graph)):
        for pn_unit_id in component:
            result[pn_unit_id] = mol_id

    return result


def compute_molecule_entity_hash(
    molecule_id: int,
    pn_unit_ids: Sequence[str],
    pn_unit_entities: Mapping[str, int],
    inter_pn_unit_bonds: Sequence[tuple[str, str]],
    inter_level_bonds: Sequence[tuple[str, str, str, str, str, str, str, str]]
    | None = None,
) -> str:
    """Compute molecule entity hash using WL graph hashing.

    Builds a graph where:
    - Nodes = PN units in the molecule
    - Edges = inter-PN-unit bonds
    - Node attribute = pn_unit_entity

    Args:
        molecule_id: The molecule ID (for filtering, not used in hash).
        pn_unit_ids: List of PN unit IDs in this molecule.
        pn_unit_entities: Mapping from pn_unit_id to pn_unit_entity.
        inter_pn_unit_bonds: List of bonds between PN units in this molecule.
        inter_level_bonds: Optional list of atomic-level bond info for inter-level hash.

    Returns:
        WL hash string for the molecule entity.
    """
    _check_networkx()

    if not pn_unit_ids:
        return ""

    pn_set = set(pn_unit_ids)
    pn_to_idx = {p: i for i, p in enumerate(pn_unit_ids)}

    # Build graph
    graph = nx.Graph()
    graph.add_nodes_from(range(len(pn_unit_ids)))

    # Add edges for bonds between PN units in this molecule
    for pn_a, pn_b in inter_pn_unit_bonds:
        if pn_a in pn_set and pn_b in pn_set and pn_a != pn_b:
            graph.add_edge(pn_to_idx[pn_a], pn_to_idx[pn_b])

    # Set node attributes (PN unit entities)
    node_attrs = {i: str(pn_unit_entities.get(p, 0)) for i, p in enumerate(pn_unit_ids)}
    nx.set_node_attributes(graph, node_attrs, "node")

    wl_hash = hash_graph(graph, node_attr="node")

    # Add inter-level bond hash if provided
    if inter_level_bonds:
        inter_hash = generate_inter_level_bond_hash(inter_level_bonds)
        if inter_hash:
            wl_hash += "_" + inter_hash

    return wl_hash


class EntityAssigner:
    """Assigns entity IDs at all levels using AtomWorks-compatible algorithm.

    This class orchestrates entity assignment across chain, pn_unit, and molecule levels.
    """

    def __init__(self) -> None:
        """Initialize the entity assigner."""
        self.chain_entity_mapper = HashToIntMapper()
        self.pn_unit_entity_mapper = HashToIntMapper()
        self.molecule_entity_mapper = HashToIntMapper()

    def assign_chain_entity(
        self,
        residues: Sequence[ResKey],
        residue_names: Mapping[ResKey, str],
        intra_bonds: Sequence[tuple[ResKey, ResKey]],
    ) -> int:
        """Assign chain entity ID.

        Args:
            residues: Ordered list of residue keys in the chain.
            residue_names: Mapping from ResKey to residue name.
            intra_bonds: List of intra-chain bonds.

        Returns:
            Chain entity ID (0-indexed integer).
        """
        hash_value = compute_chain_entity_hash(residues, residue_names, intra_bonds)
        return self.chain_entity_mapper.get_or_create(hash_value)

    def assign_pn_unit_entity(
        self,
        pn_unit_id: str,
        chain_entities: Mapping[str, int],
        inter_chain_bonds: Sequence[tuple[str, str]],
        inter_level_bonds: Sequence[tuple[str, str, str, str, str, str, str, str]]
        | None = None,
    ) -> int:
        """Assign PN unit entity ID.

        Args:
            pn_unit_id: The PN unit ID.
            chain_entities: Mapping from chain_id to chain_entity.
            inter_chain_bonds: List of inter-chain bonds.
            inter_level_bonds: Optional atomic-level bond info.

        Returns:
            PN unit entity ID (0-indexed integer).
        """
        hash_value = compute_pn_unit_entity_hash(
            pn_unit_id, chain_entities, inter_chain_bonds, inter_level_bonds
        )
        return self.pn_unit_entity_mapper.get_or_create(hash_value)

    def assign_molecule_entity(
        self,
        molecule_id: int,
        pn_unit_ids: Sequence[str],
        pn_unit_entities: Mapping[str, int],
        inter_pn_unit_bonds: Sequence[tuple[str, str]],
        inter_level_bonds: Sequence[tuple[str, str, str, str, str, str, str, str]]
        | None = None,
    ) -> int:
        """Assign molecule entity ID.

        Args:
            molecule_id: The molecule ID.
            pn_unit_ids: List of PN unit IDs in this molecule.
            pn_unit_entities: Mapping from pn_unit_id to pn_unit_entity.
            inter_pn_unit_bonds: List of bonds between PN units.
            inter_level_bonds: Optional atomic-level bond info.

        Returns:
            Molecule entity ID (0-indexed integer).
        """
        hash_value = compute_molecule_entity_hash(
            molecule_id, pn_unit_ids, pn_unit_entities, inter_pn_unit_bonds, inter_level_bonds
        )
        return self.molecule_entity_mapper.get_or_create(hash_value)

    def reset(self) -> None:
        """Reset all mappers to initial state."""
        self.chain_entity_mapper.reset()
        self.pn_unit_entity_mapper.reset()
        self.molecule_entity_mapper.reset()


__all__ = [
    "compute_chain_entity_hash",
    "compute_pn_unit_entity_hash",
    "compute_molecule_entity_hash",
    "find_pn_units",
    "find_molecules",
    "EntityAssigner",
]
