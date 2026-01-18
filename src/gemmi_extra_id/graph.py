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

# AtomWorks-compatible: only "covale" by default (disulf is not considered intermolecular)
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


def find_pn_units(
    nodes: Iterable[str],
    edges: Iterable[tuple[str, str]],
    entity_types: dict[str, str],
    is_polymer: dict[str, bool] | None = None,
) -> dict[str, str]:
    """
    Find pn_units (same-type connected components for non-polymers).

    AtomWorks-compatible behavior:
    - Polymer chains always have pn_unit_id = chain_id (never grouped)
    - Non-polymer chains are grouped by entity_type and connected components

    Args:
        nodes: Iterable of node identifiers (chain IDs).
        edges: Iterable of (node_a, node_b) pairs representing covalent connections.
        entity_types: Dict mapping each node to its entity type.
        is_polymer: Dict mapping each node to whether it's a polymer chain.
            If None, uses legacy behavior (all chains can be grouped by type).

    Returns:
        Dict mapping each node to its pn_unit_id (comma-separated sorted list of
        chain IDs in the same pn_unit).
    """
    node_list = list(nodes)
    edge_list = list(edges)

    pn_unit_mapping: dict[str, str] = {}

    # If is_polymer is provided, use AtomWorks-compatible behavior
    if is_polymer is not None:
        # Polymer chains: each is its own pn_unit
        for node in node_list:
            if is_polymer.get(node, False):
                pn_unit_mapping[node] = node

        # Non-polymer chains: find connected components (NO entity_type grouping)
        # AtomWorks groups non-polymers by connectivity only, regardless of entity_type
        non_polymer_nodes = [n for n in node_list if not is_polymer.get(n, False)]
        non_polymer_set = set(non_polymer_nodes)

        # Filter edges to non-polymer chains only
        non_polymer_edges = [
            (a, b) for a, b in edge_list if a in non_polymer_set and b in non_polymer_set
        ]

        # Find connected components among all non-polymers directly
        components = find_components(non_polymer_nodes, non_polymer_edges)

        component_to_nodes: dict[int, list[str]] = defaultdict(list)
        for node, comp_id in components.items():
            component_to_nodes[comp_id].append(node)

        for comp_nodes in component_to_nodes.values():
            pn_unit_id = ",".join(sorted(comp_nodes))
            for node in comp_nodes:
                pn_unit_mapping[node] = pn_unit_id

        return pn_unit_mapping

    # Legacy behavior: group all chains by entity_type (for backward compatibility)
    type_to_nodes = defaultdict(list)
    for node in node_list:
        etype = entity_types.get(node, "unknown")
        type_to_nodes[etype].append(node)

    for _etype, typed_nodes in type_to_nodes.items():
        typed_node_set = set(typed_nodes)
        typed_edges = [(a, b) for a, b in edge_list if a in typed_node_set and b in typed_node_set]
        components = find_components(typed_nodes, typed_edges)

        component_to_nodes: dict[int, list[str]] = defaultdict(list)
        for node, comp_id in components.items():
            component_to_nodes[comp_id].append(node)

        for comp_nodes in component_to_nodes.values():
            pn_unit_id = ",".join(sorted(comp_nodes))
            for node in comp_nodes:
                pn_unit_mapping[node] = pn_unit_id

    return pn_unit_mapping
