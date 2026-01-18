"""Tests for gemmi_extra_id.graph module."""

from collections import OrderedDict

from gemmi_extra_id.graph import find_components, find_pn_units


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


class TestFindPnUnits:
    """Tests for find_pn_units function."""

    def test_single_type_no_edges(self) -> None:
        """Each chain forms its own pn_unit when isolated."""
        nodes = ["A", "B", "C"]
        edges: list[tuple[str, str]] = []
        entity_types = {"A": "polymer", "B": "polymer", "C": "polymer"}
        result = find_pn_units(nodes, edges, entity_types)
        assert result == {"A": "A", "B": "B", "C": "C"}

    def test_single_type_with_edges_legacy(self) -> None:
        """Legacy behavior: Chains of same type connected by edges share pn_unit."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B")]
        entity_types = {"A": "polymer", "B": "polymer", "C": "polymer"}
        # Without is_polymer param, uses legacy behavior
        result = find_pn_units(nodes, edges, entity_types)
        assert result["A"] == result["B"] == "A,B"
        assert result["C"] == "C"

    def test_mixed_types_with_edges(self) -> None:
        """Edges between different types do not affect pn_units."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B")]  # A is polymer, B is non-polymer
        entity_types = {"A": "polymer", "B": "non-polymer", "C": "polymer"}
        result = find_pn_units(nodes, edges, entity_types)
        # A and B are different types, so they are separate pn_units
        assert result["A"] == "A"
        assert result["B"] == "B"
        assert result["C"] == "C"

    def test_multiple_components_same_type_legacy(self) -> None:
        """Legacy behavior: Multiple connected components within same type."""
        nodes = ["A", "B", "C", "D"]
        edges = [("A", "B"), ("C", "D")]
        entity_types = {"A": "polymer", "B": "polymer", "C": "polymer", "D": "polymer"}
        # Without is_polymer param, uses legacy behavior
        result = find_pn_units(nodes, edges, entity_types)
        assert result["A"] == result["B"] == "A,B"
        assert result["C"] == result["D"] == "C,D"

    def test_mixed_types_complex(self) -> None:
        """Complex case with multiple entity types and edges."""
        nodes = ["A", "B", "C", "D", "E"]
        edges = [("A", "C"), ("B", "C"), ("D", "E")]  # A,B=polymer linked to C=branched
        entity_types = {
            "A": "polymer",
            "B": "polymer",
            "C": "branched",
            "D": "non-polymer",
            "E": "non-polymer",
        }
        result = find_pn_units(nodes, edges, entity_types)
        # Within polymer type: A and B have no direct edge
        assert result["A"] == "A"
        assert result["B"] == "B"
        # branched: C is alone
        assert result["C"] == "C"
        # non-polymer: D and E are connected
        assert result["D"] == result["E"] == "D,E"

    def test_empty_input(self) -> None:
        """Empty nodes and edges returns empty mapping."""
        result = find_pn_units([], [], {})
        assert result == {}


class TestFindPnUnitsAtomWorksCompatible:
    """Tests for find_pn_units with AtomWorks-compatible behavior (is_polymer param)."""

    def test_polymers_never_grouped(self) -> None:
        """AtomWorks: Polymers are never grouped, even when connected."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B")]  # A-B connected by covale
        entity_types = {"A": "polymer", "B": "polymer", "C": "polymer"}
        is_polymer = {"A": True, "B": True, "C": True}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        # Each polymer is its own pn_unit
        assert result["A"] == "A"
        assert result["B"] == "B"
        assert result["C"] == "C"

    def test_non_polymers_can_be_grouped(self) -> None:
        """AtomWorks: Non-polymers can be grouped when connected."""
        nodes = ["A", "B", "C", "D"]
        edges = [("C", "D")]  # C-D connected (non-polymers)
        entity_types = {
            "A": "polymer",
            "B": "polymer",
            "C": "non-polymer",
            "D": "non-polymer",
        }
        is_polymer = {"A": True, "B": True, "C": False, "D": False}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        # Polymers are separate
        assert result["A"] == "A"
        assert result["B"] == "B"
        # Non-polymers are grouped
        assert result["C"] == result["D"] == "C,D"

    def test_polymer_nonpolymer_edge_ignored(self) -> None:
        """AtomWorks: Edge between polymer and non-polymer doesn't affect grouping."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B"), ("B", "C")]  # A(polymer)-B(non-polymer)-C(non-polymer)
        entity_types = {"A": "polymer", "B": "non-polymer", "C": "non-polymer"}
        is_polymer = {"A": True, "B": False, "C": False}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        # Polymer A is separate
        assert result["A"] == "A"
        # Non-polymers B-C are grouped
        assert result["B"] == result["C"] == "B,C"

    def test_488d_case_polymer_covale(self) -> None:
        """488d.cif case: chains B-D are polymers connected by covale."""
        nodes = ["A", "B", "C", "D"]
        edges = [("B", "D")]  # B-D polymer covale
        entity_types = {"A": "polymer", "B": "polymer", "C": "polymer", "D": "polymer"}
        is_polymer = {"A": True, "B": True, "C": True, "D": True}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        # Each polymer is separate (AtomWorks behavior)
        assert result["A"] == "A"
        assert result["B"] == "B"
        assert result["C"] == "C"
        assert result["D"] == "D"

    def test_mixed_types_non_polymer_grouping(self) -> None:
        """Non-polymers of different entity types ARE grouped when connected."""
        nodes = ["A", "B", "C", "D"]
        edges = [("B", "C")]  # B(branched)-C(non-polymer) connected
        entity_types = {
            "A": "polymer",
            "B": "branched",
            "C": "non-polymer",
            "D": "non-polymer",
        }
        is_polymer = {"A": True, "B": False, "C": False, "D": False}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        # Polymer is separate
        assert result["A"] == "A"
        # branched and non-polymer are grouped by connectivity (not entity_type)
        assert result["B"] == result["C"] == "B,C"
        # D is isolated
        assert result["D"] == "D"

    def test_non_polymers_connected_grouped(self) -> None:
        """Non-polymers connected by edges are grouped regardless of entity_type."""
        nodes = ["A", "B", "C", "D", "E"]
        edges = [("B", "C"), ("D", "E")]  # B-C, D-E both non-polymer
        entity_types = {
            "A": "polymer",
            "B": "non-polymer",
            "C": "non-polymer",
            "D": "non-polymer",
            "E": "non-polymer",
        }
        is_polymer = {"A": True, "B": False, "C": False, "D": False, "E": False}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        assert result["A"] == "A"
        assert result["B"] == result["C"] == "B,C"
        assert result["D"] == result["E"] == "D,E"

    def test_146d_case_branched_nonpolymer_grouped(self) -> None:
        """146d.cif case: branched (C,D) and non-polymer (H) connected via covale."""
        # In 146d, chains C and D are branched, H is non-polymer
        # H-C and H-D are connected via covale bonds
        # AtomWorks groups them all together: pn_unit_id = "C,D,H"
        nodes = ["A", "B", "C", "D", "H"]
        edges = [("C", "H"), ("D", "H")]  # H connects to both C and D
        entity_types = {
            "A": "polymer",
            "B": "polymer",
            "C": "branched",
            "D": "branched",
            "H": "non-polymer",
        }
        is_polymer = {"A": True, "B": True, "C": False, "D": False, "H": False}
        result = find_pn_units(nodes, edges, entity_types, is_polymer)
        # Polymers are separate
        assert result["A"] == "A"
        assert result["B"] == "B"
        # All non-polymers connected form one group
        assert result["C"] == result["D"] == result["H"] == "C,D,H"
