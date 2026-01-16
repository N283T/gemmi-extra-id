"""Graph algorithms for molecular connectivity analysis.

This module provides algorithms for finding connected components in molecular
structures based on covalent bond connectivity. It is used to assign molecule_id
and pn_unit_id to chains in mmCIF files.

The main functions are:
- find_components: Find connected components using iterative DFS
- find_pn_units: Find same-type connected components (pn_units)

Example:
    >>> from gemmi_extra_id.graph import find_components
    >>> nodes = ["A", "B", "C"]
    >>> edges = [("A", "B")]
    >>> find_components(nodes, edges)
    OrderedDict([('A', 0), ('B', 0), ('C', 1)])
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from collections.abc import Iterable
from typing import TypeAlias

# Type aliases for clarity
NodeId: TypeAlias = str
ComponentId: TypeAlias = int
Edge: TypeAlias = tuple[NodeId, NodeId]
EntityType: TypeAlias = str
PnUnitId: TypeAlias = str  # Comma-separated chain IDs

DEFAULT_COVALENT_TYPES: frozenset[str] = frozenset({"covale", "disulf"})


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


def find_pn_units(
    nodes: Iterable[str],
    edges: Iterable[tuple[str, str]],
    entity_types: dict[str, str],
) -> dict[str, str]:
    """
    Find pn_units (same-type connected components).

    A pn_unit is a group of covalently linked chains that share the same entity type
    (e.g., polymer, non-polymer, branched, water).

    Args:
        nodes: Iterable of node identifiers (chain IDs).
        edges: Iterable of (node_a, node_b) pairs representing covalent connections.
        entity_types: Dict mapping each node to its entity type.

    Returns:
        Dict mapping each node to its pn_unit_id (comma-separated sorted list of
        chain IDs in the same pn_unit).
    """
    node_list = list(nodes)
    edge_list = list(edges)

    # Group nodes by entity_type
    type_to_nodes: dict[str, list[str]] = defaultdict(list)
    for node in node_list:
        etype = entity_types.get(node, "unknown")
        type_to_nodes[etype].append(node)

    # Build pn_unit_id mapping
    pn_unit_mapping: dict[str, str] = {}

    for _etype, typed_nodes in type_to_nodes.items():
        typed_node_set = set(typed_nodes)

        # Filter edges to only include edges between nodes of this type
        typed_edges = [(a, b) for a, b in edge_list if a in typed_node_set and b in typed_node_set]

        # Find connected components within this type
        components = find_components(typed_nodes, typed_edges)

        # Group nodes by component
        component_to_nodes: dict[int, list[str]] = defaultdict(list)
        for node, comp_id in components.items():
            component_to_nodes[comp_id].append(node)

        # Assign pn_unit_id as comma-separated sorted chain list
        for comp_nodes in component_to_nodes.values():
            pn_unit_id = ",".join(sorted(comp_nodes))
            for node in comp_nodes:
                pn_unit_mapping[node] = pn_unit_id

    return pn_unit_mapping
