"""Tests for cifmolid.core module."""

from collections import OrderedDict

from cifmolid.core import find_components


class TestFindComponents:
    """Tests for find_components function."""

    def test_single_node_no_edges(self) -> None:
        """Single node with no edges forms one component."""
        result = find_components(["A"], [])
        assert result == OrderedDict([("A", 0)])

    def test_multiple_isolated_nodes(self) -> None:
        """Multiple isolated nodes each form their own component."""
        result = find_components(["A", "B", "C"], [])
        assert result == OrderedDict([("A", 0), ("B", 1), ("C", 2)])

    def test_two_connected_nodes(self) -> None:
        """Two connected nodes form one component."""
        result = find_components(["A", "B"], [("A", "B")])
        assert result == OrderedDict([("A", 0), ("B", 0)])

    def test_chain_connectivity(self) -> None:
        """Chain A-B-C forms one component."""
        result = find_components(["A", "B", "C"], [("A", "B"), ("B", "C")])
        assert result == OrderedDict([("A", 0), ("B", 0), ("C", 0)])

    def test_two_separate_components(self) -> None:
        """Two separate groups form two components."""
        result = find_components(["A", "B", "C", "D"], [("A", "B"), ("C", "D")])
        assert result == OrderedDict([("A", 0), ("B", 0), ("C", 1), ("D", 1)])

    def test_preserves_node_order(self) -> None:
        """Component IDs are assigned in node encounter order."""
        result = find_components(["X", "Y", "Z"], [("Z", "Y")])
        # X is first, so component 0
        # Y is second, connected to Z, so component 1
        # Z is third, connected to Y, so component 1
        assert result == OrderedDict([("X", 0), ("Y", 1), ("Z", 1)])

    def test_edge_with_unknown_node(self) -> None:
        """Edges can reference nodes not in initial list."""
        result = find_components(["A"], [("A", "B")])
        assert result == OrderedDict([("A", 0), ("B", 0)])

    def test_empty_input(self) -> None:
        """Empty nodes and edges returns empty mapping."""
        result = find_components([], [])
        assert result == OrderedDict()

    def test_cycle(self) -> None:
        """Cyclic graph forms one component."""
        result = find_components(["A", "B", "C"], [("A", "B"), ("B", "C"), ("C", "A")])
        assert set(result.keys()) == {"A", "B", "C"}
        assert result["A"] == result["B"] == result["C"] == 0

    def test_complex_graph(self) -> None:
        """Complex graph with multiple components."""
        nodes = ["A", "B", "C", "D", "E", "F"]
        edges = [("A", "B"), ("B", "C"), ("D", "E")]  # A-B-C and D-E, F isolated
        result = find_components(nodes, edges)

        assert result["A"] == result["B"] == result["C"]
        assert result["D"] == result["E"]
        assert result["F"] != result["A"]
        assert result["F"] != result["D"]
