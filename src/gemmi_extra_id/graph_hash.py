"""Graph hashing for entity identification.

This module provides Weisfeiler-Lehman graph hashing for computing entity IDs
in a manner compatible with AtomWorks. The hashing is hierarchical:

1. chain_entity: Hash of residue graph (nodes=residues, attr=res_name)
2. pn_unit_entity: Hash of chain graph (nodes=chains, attr=chain_entity)
3. molecule_entity: Hash of pn_unit graph (nodes=pn_units, attr=pn_unit_entity)

This module requires networkx. Install with: pip install gemmi-extra-id[hash]
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Mapping
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence


def _check_networkx_available() -> None:
    """Check if networkx is available, raise ImportError if not."""
    try:
        import networkx  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "networkx is required for graph hashing. Install with: pip install gemmi-extra-id[hash]"
        ) from e


def hash_graph(
    nodes: Sequence[Hashable],
    edges: Iterable[tuple[Hashable, Hashable]],
    node_attrs: Mapping[Hashable, str],
    iterations: int = 3,
) -> str:
    """
    Compute Weisfeiler-Lehman graph hash.

    Args:
        nodes: Sequence of node identifiers.
        edges: Iterable of (node_a, node_b) pairs representing edges.
        node_attrs: Mapping from node to its attribute value (used for hashing).
        iterations: Number of WL iterations (default 3).

    Returns:
        Graph hash string combining WL hash, node count, and attribute distribution.
    """
    _check_networkx_available()
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # Set node attributes for WL hashing
    for node in G.nodes():
        G.nodes[node]["label"] = node_attrs.get(node, "")

    # Compute WL hash
    wl_hash = nx.weisfeiler_lehman_graph_hash(G, node_attr="label", iterations=iterations)

    # Add node count and attribute distribution for AtomWorks compatibility
    n_nodes = len(nodes)
    attr_counts = Counter(node_attrs.get(n, "") for n in nodes)
    attr_str = "_".join(f"{k}:{v}" for k, v in sorted(attr_counts.items()))

    return f"{wl_hash}_{n_nodes}_{attr_str}"


class KeyToIntMapper:
    """Maps arbitrary hashable keys to 0-indexed integers.

    Compatible with AtomWorks common.py KeyToIntMapper.
    """

    def __init__(self) -> None:
        self._mapping: dict[Hashable, int] = {}
        self._next_id: int = 0

    def get_id(self, key: Hashable) -> int:
        """Get or create an integer ID for the given key."""
        if key not in self._mapping:
            self._mapping[key] = self._next_id
            self._next_id += 1
        return self._mapping[key]

    def __len__(self) -> int:
        return len(self._mapping)


def compute_chain_entity(
    residues: Sequence[tuple[str, str]],
    residue_names: Mapping[tuple[str, str], str],
    intra_chain_bonds: Iterable[tuple[tuple[str, str], tuple[str, str]]],
) -> str:
    """
    Compute chain entity hash from residue graph.

    Args:
        residues: Sequence of (seq_id, ins_code) tuples identifying residues.
        residue_names: Mapping from residue ID to residue name (comp_id).
        intra_chain_bonds: Edges between residues within the chain.

    Returns:
        Chain entity hash string.
    """
    return hash_graph(
        nodes=list(residues),
        edges=intra_chain_bonds,
        node_attrs={r: residue_names.get(r, "") for r in residues},
    )


def compute_pn_unit_entity(
    chain_ids: Sequence[str],
    chain_entities: Mapping[str, str | int],
    inter_chain_bonds: Iterable[tuple[str, str]],
) -> str:
    """
    Compute pn_unit entity hash from chain graph.

    Args:
        chain_ids: Sequence of chain IDs in the pn_unit.
        chain_entities: Mapping from chain ID to its chain_entity.
        inter_chain_bonds: Edges between chains in the pn_unit.

    Returns:
        PN unit entity hash string.
    """
    return hash_graph(
        nodes=list(chain_ids),
        edges=inter_chain_bonds,
        node_attrs={c: str(chain_entities.get(c, "")) for c in chain_ids},
    )


def compute_molecule_entity(
    pn_unit_ids: Sequence[str],
    pn_unit_entities: Mapping[str, str | int],
    inter_pn_unit_bonds: Iterable[tuple[str, str]],
) -> str:
    """
    Compute molecule entity hash from pn_unit graph.

    Args:
        pn_unit_ids: Sequence of pn_unit IDs in the molecule.
        pn_unit_entities: Mapping from pn_unit ID to its pn_unit_entity.
        inter_pn_unit_bonds: Edges between pn_units in the molecule.

    Returns:
        Molecule entity hash string.
    """
    return hash_graph(
        nodes=list(pn_unit_ids),
        edges=inter_pn_unit_bonds,
        node_attrs={p: str(pn_unit_entities.get(p, "")) for p in pn_unit_ids},
    )


def compute_hash_entities(
    chain_order: Sequence[str],
    molecule_mapping: Mapping[str, int],
    pn_unit_mapping: Mapping[str, str],
    chain_residues: Mapping[str, Sequence[tuple[str, str]]],
    residue_names: Mapping[str, Mapping[tuple[str, str], str]],
    intra_chain_bonds: Mapping[str, Iterable[tuple[tuple[str, str], tuple[str, str]]]],
    inter_chain_bonds: Iterable[tuple[str, str]],
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    """
    Compute all hash-based entity IDs.

    This is the main entry point for computing hash entities. It computes
    chain_entity, pn_unit_entity, and molecule_entity for all chains using
    Weisfeiler-Lehman graph hashing.

    Args:
        chain_order: Sequence of chain IDs in order.
        molecule_mapping: Mapping from chain ID to molecule ID.
        pn_unit_mapping: Mapping from chain ID to pn_unit_id.
        chain_residues: Mapping from chain ID to list of (seq_id, ins_code) tuples.
        residue_names: Mapping from chain ID to {(seq_id, ins_code): res_name}.
        intra_chain_bonds: Mapping from chain ID to intra-chain residue bonds.
        inter_chain_bonds: All inter-chain covalent bonds.

    Returns:
        Tuple of (chain_entity_map, pn_unit_entity_map, molecule_entity_map),
        where each maps chain_id to its integer entity ID.
    """
    _check_networkx_available()

    # Step 1: Compute chain_entity for each chain
    chain_entity_mapper = KeyToIntMapper()
    chain_entity_hash: dict[str, str] = {}

    for chain in chain_order:
        residues = chain_residues.get(chain, [])
        names = residue_names.get(chain, {})
        bonds = intra_chain_bonds.get(chain, [])
        chain_entity_hash[chain] = compute_chain_entity(residues, names, bonds)

    # Map hashes to integers
    chain_entity_map: dict[str, int] = {}
    for chain in chain_order:
        chain_entity_map[chain] = chain_entity_mapper.get_id(chain_entity_hash[chain])

    # Step 2: Compute pn_unit_entity for each pn_unit
    # Group chains by pn_unit
    pn_unit_to_chains: dict[str, list[str]] = {}
    for chain in chain_order:
        pn_unit_id = pn_unit_mapping.get(chain, chain)
        if pn_unit_id not in pn_unit_to_chains:
            pn_unit_to_chains[pn_unit_id] = []
        pn_unit_to_chains[pn_unit_id].append(chain)

    # Filter inter-chain bonds to those within each pn_unit
    inter_chain_edges = list(inter_chain_bonds)

    pn_unit_entity_mapper = KeyToIntMapper()
    pn_unit_entity_hash: dict[str, str] = {}

    for pn_unit_id, chains in pn_unit_to_chains.items():
        chain_set = set(chains)
        pn_unit_bonds = [(a, b) for a, b in inter_chain_edges if a in chain_set and b in chain_set]
        pn_unit_entity_hash[pn_unit_id] = compute_pn_unit_entity(
            chains, chain_entity_map, pn_unit_bonds
        )

    # Map pn_unit hashes to chain entities
    pn_unit_entity_map: dict[str, int] = {}
    for chain in chain_order:
        pn_unit_id = pn_unit_mapping.get(chain, chain)
        pn_unit_entity_map[chain] = pn_unit_entity_mapper.get_id(pn_unit_entity_hash[pn_unit_id])

    # Step 3: Compute molecule_entity for each molecule
    # Group pn_units by molecule
    molecule_to_pn_units: dict[int, list[str]] = {}
    pn_unit_to_molecule: dict[str, int] = {}
    for chain in chain_order:
        mol_id = molecule_mapping.get(chain, 0)
        pn_unit_id = pn_unit_mapping.get(chain, chain)
        if mol_id not in molecule_to_pn_units:
            molecule_to_pn_units[mol_id] = []
        if pn_unit_id not in pn_unit_to_molecule:
            molecule_to_pn_units[mol_id].append(pn_unit_id)
            pn_unit_to_molecule[pn_unit_id] = mol_id

    # Determine inter-pn_unit bonds (bonds that cross pn_unit boundaries within a molecule)
    chain_to_pn_unit = {c: pn_unit_mapping.get(c, c) for c in chain_order}

    molecule_entity_mapper = KeyToIntMapper()
    molecule_entity_hash: dict[int, str] = {}

    for mol_id, pn_units in molecule_to_pn_units.items():
        # Find bonds between different pn_units in this molecule
        mol_chains = [c for c in chain_order if molecule_mapping.get(c, 0) == mol_id]
        mol_chain_set = set(mol_chains)

        inter_pn_unit_bonds: list[tuple[str, str]] = []
        for a, b in inter_chain_edges:
            if a in mol_chain_set and b in mol_chain_set:
                pn_a = chain_to_pn_unit[a]
                pn_b = chain_to_pn_unit[b]
                if pn_a != pn_b and (pn_a, pn_b) not in inter_pn_unit_bonds:
                    inter_pn_unit_bonds.append((pn_a, pn_b))

        # Get pn_unit_entity values for pn_units in this molecule
        pn_unit_entities_for_mol = {
            pn: pn_unit_entity_mapper.get_id(pn_unit_entity_hash[pn]) for pn in pn_units
        }

        molecule_entity_hash[mol_id] = compute_molecule_entity(
            pn_units, pn_unit_entities_for_mol, inter_pn_unit_bonds
        )

    # Map molecule hashes to chain entities
    molecule_entity_map: dict[str, int] = {}
    for chain in chain_order:
        mol_id = molecule_mapping.get(chain, 0)
        molecule_entity_map[chain] = molecule_entity_mapper.get_id(molecule_entity_hash[mol_id])

    return chain_entity_map, pn_unit_entity_map, molecule_entity_map
