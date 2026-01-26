"""Weisfeiler-Lehman graph hashing for complete mode.

This module provides AtomWorks-compatible graph hashing functions
for entity assignment.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

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


def hash_graph(
    graph: nx.Graph,
    node_attr: str | None = None,
    edge_attr: str | None = None,
    iterations: int = 3,
    digest_size: int = 16,
) -> str:
    """Compute AtomWorks-compatible Weisfeiler-Lehman graph hash.

    This function computes a hash for a given graph using the WL algorithm
    and additionally adds node/edge attribute counts to handle edge cases
    where WL fails (e.g., disconnected graphs).

    This matches AtomWorks' hash_graph() function exactly.

    Args:
        graph: NetworkX graph to hash.
        node_attr: Node attribute name to include in hashing. If None, ignored.
        edge_attr: Edge attribute name to include in hashing. If None, ignored.
        iterations: Number of WL iterations. Default is 3.
        digest_size: Size of hash digest. Default is 16.

    Returns:
        Hash string in format: "<wl_hash>_<node_count>_<attr>:<count>,..."

    Example:
        >>> import networkx as nx
        >>> G = nx.path_graph(3)
        >>> nx.set_node_attributes(G, {0: "A", 1: "B", 2: "A"}, "node")
        >>> hash_graph(G, node_attr="node")
        '<wl_hash>_3_A:2,B:1'
    """
    _check_networkx()

    # Compute base WL hash
    wl_hash = nx.algorithms.graph_hashing.weisfeiler_lehman_graph_hash(
        graph,
        node_attr=node_attr,
        edge_attr=edge_attr,
        iterations=iterations,
        digest_size=digest_size,
    )

    if node_attr is not None:
        # Add number of nodes
        wl_hash += f"_{len(graph.nodes)}"

        # Add node attribute counts (sorted for determinism)
        node_attr_dict = nx.get_node_attributes(graph, node_attr)
        attr_counts = Counter(node_attr_dict.values())
        sorted_counts = sorted(attr_counts.items(), key=lambda x: (str(x[0]), x[1]))
        wl_hash += "_" + ",".join(f"{attr}:{count}" for attr, count in sorted_counts)

    if edge_attr is not None:
        # Add number of edges
        wl_hash += f"_{len(graph.edges)}"

    return wl_hash


def generate_inter_level_bond_hash(
    bonds: Sequence[tuple[str, ...]] | Sequence[tuple[str, str, str, str, str, str]],
) -> str:
    """Generate hash for inter-level bonds.

    This function creates a hash from atomic-level bond information
    to distinguish entities that have the same graph structure but
    different bond connectivity patterns.

    Supports two tuple formats:
    - 6-tuple: (res_id1, res_name1, atom_name1, res_id2, res_name2, atom_name2)
      Used for chain entity level (no lower_level_entity)
    - 8-tuple: (entity1, res_id1, res_name1, atom_name1, entity2, res_id2, res_name2, atom_name2)
      Used for pn_unit/molecule levels (with lower_level_entity)

    Args:
        bonds: Sequence of bond tuples with atomic information.

    Returns:
        Hash string representing the bond pattern.

    Example:
        >>> # 6-tuple format (chain entity level)
        >>> bonds = [("10", "CYS", "SG", "20", "CYS", "SG")]
        >>> generate_inter_level_bond_hash(bonds)
        '<hash_string>'

        >>> # 8-tuple format (pn_unit/molecule level)
        >>> bonds = [("1", "10", "CYS", "SG", "2", "20", "CYS", "SG")]
        >>> generate_inter_level_bond_hash(bonds)
        '<hash_string>'
    """
    if not bonds:
        return ""

    # Determine tuple format from first bond
    bond_len = len(bonds[0])
    if bond_len == 6:
        half_len = 3
    elif bond_len == 8:
        half_len = 4
    else:
        raise ValueError(f"Expected 6 or 8 element bond tuple, got {bond_len}")

    # Sort bond tuples for determinism
    # Each bond is normalized so smaller tuple comes first
    normalized_bonds = []
    for bond in bonds:
        # Split into two halves
        half1 = bond[:half_len]
        half2 = bond[half_len:]
        # Sort so smaller tuple comes first
        if half1 <= half2:
            normalized_bonds.append((half1, half2))
        else:
            normalized_bonds.append((half2, half1))

    # Sort all bonds
    sorted_bonds = sorted(normalized_bonds)

    # Create hash from sorted bond tuples
    bond_str = str(sorted_bonds)
    return hashlib.md5(bond_str.encode(), usedforsecurity=False).hexdigest()


class HashToIntMapper:
    """Map hash values to consecutive integers.

    This class provides AtomWorks-compatible hash-to-integer mapping
    for entity IDs. Each unique hash gets assigned a consecutive integer
    starting from 0.
    """

    def __init__(self) -> None:
        """Initialize the mapper."""
        self._hash_to_int: dict[str, int] = {}
        self._next_id: int = 0

    def get_or_create(self, hash_value: str) -> int:
        """Get or create an integer ID for a hash value.

        Args:
            hash_value: Hash string to map.

        Returns:
            Integer ID (0-indexed, assigned in order of first appearance).
        """
        if hash_value not in self._hash_to_int:
            self._hash_to_int[hash_value] = self._next_id
            self._next_id += 1
        return self._hash_to_int[hash_value]

    def get(self, hash_value: str) -> int | None:
        """Get the integer ID for a hash value, or None if not mapped.

        Args:
            hash_value: Hash string to look up.

        Returns:
            Integer ID if mapped, None otherwise.
        """
        return self._hash_to_int.get(hash_value)

    @property
    def count(self) -> int:
        """Number of unique hashes mapped."""
        return self._next_id

    def reset(self) -> None:
        """Reset the mapper to empty state."""
        self._hash_to_int.clear()
        self._next_id = 0


__all__ = [
    "hash_graph",
    "generate_inter_level_bond_hash",
    "HashToIntMapper",
    "HAS_NETWORKX",
]
