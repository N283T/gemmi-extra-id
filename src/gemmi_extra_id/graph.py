"""Graph algorithms for molecular connectivity analysis.

This module provides algorithms for finding connected components in molecular
structures based on covalent bond connectivity. It is used to assign molecule_id
to chains in mmCIF files.

Example:
    >>> from gemmi_extra_id.graph import find_components
    >>> nodes = ["A", "B", "C"]
    >>> edges = [("A", "B")]
    >>> find_components(nodes, edges)
    OrderedDict([('A', 0), ('B', 0), ('C', 1)])
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterable

# AtomWorks-compatible: only "covale" by default
DEFAULT_COVALENT_TYPES: frozenset[str] = frozenset({"covale"})


def find_components(
    nodes: Iterable[str],
    edges: Iterable[tuple[str, str]],
) -> OrderedDict[str, int]:
    """
    Find connected components using iterative DFS.

    Args:
        nodes: Iterable of node identifiers (preserves order for component numbering).
        edges: Iterable of (node_a, node_b) pairs representing connections.

    Returns:
        OrderedDict mapping each node to its component ID (0-indexed).
        Component IDs are assigned in the order nodes are first encountered.
    """
    node_list = list(nodes)
    adjacency: dict[str, set[str]] = {node: set() for node in node_list}

    for a, b in edges:
        if a not in adjacency:
            adjacency[a] = set()
            node_list.append(a)
        if b not in adjacency:
            adjacency[b] = set()
            node_list.append(b)
        adjacency[a].add(b)
        adjacency[b].add(a)

    mapping: OrderedDict[str, int] = OrderedDict()
    component_id = 0

    for node in node_list:
        if node in mapping:
            continue
        stack = [node]
        while stack:
            current = stack.pop()
            if current in mapping:
                continue
            mapping[current] = component_id
            stack.extend(adjacency.get(current, set()))
        component_id += 1

    return mapping
