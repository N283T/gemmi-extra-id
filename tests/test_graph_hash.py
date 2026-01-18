"""Tests for gemmi_extra_id.graph_hash module."""

import pytest

# Skip all tests if networkx is not available
networkx = pytest.importorskip("networkx")

# ruff: noqa: E402
from gemmi_extra_id.graph_hash import (
    KeyToIntMapper,
    compute_chain_entity,
    compute_hash_entities,
    compute_molecule_entity,
    compute_pn_unit_entity,
    hash_graph,
)


class TestHashGraph:
    """Tests for hash_graph function."""

    def test_empty_graph(self) -> None:
        """Empty graph produces a hash."""
        h = hash_graph(nodes=[], edges=[], node_attrs={})
        assert isinstance(h, str)
        assert len(h) > 0

    def test_single_node(self) -> None:
        """Single node graph produces a hash."""
        h = hash_graph(nodes=["A"], edges=[], node_attrs={"A": "ALA"})
        assert isinstance(h, str)
        assert "1" in h  # node count included

    def test_two_nodes_same_attr(self) -> None:
        """Two nodes with same attribute produce a hash."""
        h = hash_graph(nodes=["A", "B"], edges=[], node_attrs={"A": "ALA", "B": "ALA"})
        assert isinstance(h, str)
        assert "2" in h  # node count included

    def test_two_nodes_different_attr(self) -> None:
        """Two nodes with different attributes produce a hash."""
        h = hash_graph(nodes=["A", "B"], edges=[], node_attrs={"A": "ALA", "B": "GLY"})
        assert isinstance(h, str)

    def test_connected_vs_disconnected(self) -> None:
        """Connected and disconnected graphs produce different hashes."""
        h1 = hash_graph(nodes=["A", "B"], edges=[("A", "B")], node_attrs={"A": "ALA", "B": "ALA"})
        h2 = hash_graph(nodes=["A", "B"], edges=[], node_attrs={"A": "ALA", "B": "ALA"})
        assert h1 != h2

    def test_same_structure_same_hash(self) -> None:
        """Isomorphic graphs with same attributes produce same hash."""
        h1 = hash_graph(
            nodes=["A", "B", "C"],
            edges=[("A", "B"), ("B", "C")],
            node_attrs={"A": "ALA", "B": "GLY", "C": "ALA"},
        )
        h2 = hash_graph(
            nodes=["X", "Y", "Z"],
            edges=[("X", "Y"), ("Y", "Z")],
            node_attrs={"X": "ALA", "Y": "GLY", "Z": "ALA"},
        )
        assert h1 == h2


class TestKeyToIntMapper:
    """Tests for KeyToIntMapper class."""

    def test_first_key_gets_zero(self) -> None:
        """First key gets ID 0."""
        mapper = KeyToIntMapper()
        assert mapper.get_id("first") == 0

    def test_same_key_same_id(self) -> None:
        """Same key always returns same ID."""
        mapper = KeyToIntMapper()
        id1 = mapper.get_id("key")
        id2 = mapper.get_id("key")
        assert id1 == id2 == 0

    def test_different_keys_different_ids(self) -> None:
        """Different keys get different IDs."""
        mapper = KeyToIntMapper()
        id1 = mapper.get_id("key1")
        id2 = mapper.get_id("key2")
        assert id1 != id2
        assert id1 == 0
        assert id2 == 1

    def test_sequential_ids(self) -> None:
        """IDs are assigned sequentially."""
        mapper = KeyToIntMapper()
        for i in range(5):
            assert mapper.get_id(f"key{i}") == i

    def test_len(self) -> None:
        """len() returns number of unique keys."""
        mapper = KeyToIntMapper()
        mapper.get_id("a")
        mapper.get_id("b")
        mapper.get_id("a")  # duplicate
        assert len(mapper) == 2


class TestComputeChainEntity:
    """Tests for compute_chain_entity function."""

    def test_single_residue(self) -> None:
        """Single residue chain produces a hash."""
        residues = [("1", ".")]
        names = {("1", "."): "ALA"}
        bonds: list[tuple[tuple[str, str], tuple[str, str]]] = []
        h = compute_chain_entity(residues, names, bonds)
        assert isinstance(h, str)

    def test_polymer_chain(self) -> None:
        """Polymer chain with multiple residues."""
        residues = [("1", "."), ("2", "."), ("3", ".")]
        names = {("1", "."): "ALA", ("2", "."): "GLY", ("3", "."): "VAL"}
        bonds = [(("1", "."), ("2", ".")), (("2", "."), ("3", "."))]
        h = compute_chain_entity(residues, names, bonds)
        assert isinstance(h, str)

    def test_same_sequence_same_hash(self) -> None:
        """Same sequence produces same hash."""
        residues = [("1", "."), ("2", ".")]
        names = {("1", "."): "ALA", ("2", "."): "GLY"}
        bonds = [(("1", "."), ("2", "."))]
        h1 = compute_chain_entity(residues, names, bonds)
        h2 = compute_chain_entity(residues, names, bonds)
        assert h1 == h2


class TestComputePnUnitEntity:
    """Tests for compute_pn_unit_entity function."""

    def test_single_chain(self) -> None:
        """Single chain pn_unit."""
        h = compute_pn_unit_entity(["A"], {"A": 0}, [])
        assert isinstance(h, str)

    def test_multiple_chains(self) -> None:
        """Multiple chains in pn_unit."""
        h = compute_pn_unit_entity(["A", "B"], {"A": 0, "B": 1}, [("A", "B")])
        assert isinstance(h, str)


class TestComputeMoleculeEntity:
    """Tests for compute_molecule_entity function."""

    def test_single_pn_unit(self) -> None:
        """Single pn_unit molecule."""
        h = compute_molecule_entity(["A"], {"A": 0}, [])
        assert isinstance(h, str)

    def test_multiple_pn_units(self) -> None:
        """Multiple pn_units in molecule."""
        h = compute_molecule_entity(["A", "B"], {"A": 0, "B": 1}, [("A", "B")])
        assert isinstance(h, str)


class TestComputeHashEntities:
    """Tests for compute_hash_entities function."""

    def test_simple_structure(self) -> None:
        """Simple structure with one chain."""
        chain_order = ["A"]
        molecule_mapping = {"A": 0}
        pn_unit_mapping = {"A": "A"}
        chain_residues = {"A": [("1", "."), ("2", ".")]}
        residue_names = {"A": {("1", "."): "ALA", ("2", "."): "GLY"}}
        intra_chain_bonds: dict[str, list[tuple[tuple[str, str], tuple[str, str]]]] = {
            "A": [(("1", "."), ("2", "."))]
        }
        inter_chain_bonds: list[tuple[str, str]] = []

        chain_entity, pn_unit_entity, molecule_entity = compute_hash_entities(
            chain_order=chain_order,
            molecule_mapping=molecule_mapping,
            pn_unit_mapping=pn_unit_mapping,
            chain_residues=chain_residues,
            residue_names=residue_names,
            intra_chain_bonds=intra_chain_bonds,
            inter_chain_bonds=inter_chain_bonds,
        )

        assert "A" in chain_entity
        assert "A" in pn_unit_entity
        assert "A" in molecule_entity
        assert chain_entity["A"] == 0  # First unique chain entity
        assert pn_unit_entity["A"] == 0  # First unique pn_unit entity
        assert molecule_entity["A"] == 0  # First unique molecule entity

    def test_two_identical_chains(self) -> None:
        """Two identical chains get same entity IDs."""
        chain_order = ["A", "B"]
        molecule_mapping = {"A": 0, "B": 1}  # Different molecules
        pn_unit_mapping = {"A": "A", "B": "B"}
        chain_residues = {
            "A": [("1", "."), ("2", ".")],
            "B": [("1", "."), ("2", ".")],
        }
        residue_names = {
            "A": {("1", "."): "ALA", ("2", "."): "GLY"},
            "B": {("1", "."): "ALA", ("2", "."): "GLY"},
        }
        intra_chain_bonds: dict[str, list[tuple[tuple[str, str], tuple[str, str]]]] = {
            "A": [(("1", "."), ("2", "."))],
            "B": [(("1", "."), ("2", "."))],
        }
        inter_chain_bonds: list[tuple[str, str]] = []

        chain_entity, pn_unit_entity, molecule_entity = compute_hash_entities(
            chain_order=chain_order,
            molecule_mapping=molecule_mapping,
            pn_unit_mapping=pn_unit_mapping,
            chain_residues=chain_residues,
            residue_names=residue_names,
            intra_chain_bonds=intra_chain_bonds,
            inter_chain_bonds=inter_chain_bonds,
        )

        # Identical chains should have same entity IDs
        assert chain_entity["A"] == chain_entity["B"]
        assert pn_unit_entity["A"] == pn_unit_entity["B"]
        assert molecule_entity["A"] == molecule_entity["B"]

    def test_two_different_chains(self) -> None:
        """Two different chains get different entity IDs."""
        chain_order = ["A", "B"]
        molecule_mapping = {"A": 0, "B": 1}
        pn_unit_mapping = {"A": "A", "B": "B"}
        chain_residues = {
            "A": [("1", "."), ("2", ".")],
            "B": [("1", "."), ("2", "."), ("3", ".")],  # Different length
        }
        residue_names = {
            "A": {("1", "."): "ALA", ("2", "."): "GLY"},
            "B": {("1", "."): "ALA", ("2", "."): "GLY", ("3", "."): "VAL"},
        }
        intra_chain_bonds: dict[str, list[tuple[tuple[str, str], tuple[str, str]]]] = {
            "A": [(("1", "."), ("2", "."))],
            "B": [(("1", "."), ("2", ".")), (("2", "."), ("3", "."))],
        }
        inter_chain_bonds: list[tuple[str, str]] = []

        chain_entity, pn_unit_entity, molecule_entity = compute_hash_entities(
            chain_order=chain_order,
            molecule_mapping=molecule_mapping,
            pn_unit_mapping=pn_unit_mapping,
            chain_residues=chain_residues,
            residue_names=residue_names,
            intra_chain_bonds=intra_chain_bonds,
            inter_chain_bonds=inter_chain_bonds,
        )

        # Different chains should have different entity IDs
        assert chain_entity["A"] != chain_entity["B"]
        assert pn_unit_entity["A"] != pn_unit_entity["B"]
        assert molecule_entity["A"] != molecule_entity["B"]

    def test_connected_molecule(self) -> None:
        """Two chains in same molecule via inter-chain bond."""
        chain_order = ["A", "B"]
        molecule_mapping = {"A": 0, "B": 0}  # Same molecule
        pn_unit_mapping = {"A": "A", "B": "B"}  # Different pn_units
        chain_residues = {
            "A": [("1", ".")],
            "B": [("1", ".")],
        }
        residue_names = {
            "A": {("1", "."): "ALA"},
            "B": {("1", "."): "GLY"},
        }
        intra_chain_bonds: dict[str, list[tuple[tuple[str, str], tuple[str, str]]]] = {}
        inter_chain_bonds = [("A", "B")]

        chain_entity, pn_unit_entity, molecule_entity = compute_hash_entities(
            chain_order=chain_order,
            molecule_mapping=molecule_mapping,
            pn_unit_mapping=pn_unit_mapping,
            chain_residues=chain_residues,
            residue_names=residue_names,
            intra_chain_bonds=intra_chain_bonds,
            inter_chain_bonds=inter_chain_bonds,
        )

        # Different chains, different chain entities
        assert chain_entity["A"] != chain_entity["B"]
        # Different pn_units, different pn_unit entities
        assert pn_unit_entity["A"] != pn_unit_entity["B"]
        # Same molecule, same molecule entity
        assert molecule_entity["A"] == molecule_entity["B"]
