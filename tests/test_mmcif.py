"""Tests for gemmi_extra_id.mmcif module."""

import json
from pathlib import Path

import gemmi
import pytest

from gemmi_extra_id.mmcif import (
    VALID_SWAP_TARGETS,
    assign_extended_ids,
    assign_molecule_id,
    swap_auth_asym_id,
)

DATA_DIR = Path(__file__).parent / "data"
FROM_PDB_DIR = DATA_DIR / "from_pdb"
REF_DIR = DATA_DIR / "ref"


def _has_test_data(filename: str) -> bool:
    """Check if test data file exists."""
    return (FROM_PDB_DIR / filename).exists()


class TestAssignMoleculeId:
    """Tests for assign_molecule_id function."""

    def test_file_not_found(self) -> None:
        """Raises FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError):
            assign_molecule_id("nonexistent.cif")

    def test_directory_input(self, tmp_path: Path) -> None:
        """Raises ValueError for directory input."""
        with pytest.raises(ValueError, match="not a file"):
            assign_molecule_id(tmp_path)

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_returns_ordered_dict(self) -> None:
        """Returns OrderedDict with correct structure."""
        mapping = assign_molecule_id(FROM_PDB_DIR / "148L.cif")

        assert isinstance(mapping, dict)
        assert len(mapping) > 0
        assert all(isinstance(k, str) for k in mapping)
        assert all(isinstance(v, int) for v in mapping.values())

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_writes_output_file(self, tmp_path: Path) -> None:
        """Writes output file when path is provided."""
        output_file = tmp_path / "output.cif"

        assign_molecule_id(FROM_PDB_DIR / "148L.cif", output_file)

        assert output_file.exists()
        assert output_file.stat().st_size > 0

    @pytest.mark.skipif(not _has_test_data("148L.cif.gz"), reason="Test data not available")
    def test_reads_gzipped_cif(self) -> None:
        """Can read gzipped CIF files."""
        mapping = assign_molecule_id(FROM_PDB_DIR / "148L.cif.gz")

        assert isinstance(mapping, dict)
        assert len(mapping) > 0

    @pytest.mark.skipif(
        not (_has_test_data("148L.cif") and _has_test_data("148L.cif.gz")),
        reason="Test data not available",
    )
    def test_cif_and_gzip_produce_same_result(self) -> None:
        """CIF and gzipped CIF produce identical results."""
        mapping_cif = assign_molecule_id(FROM_PDB_DIR / "148L.cif")
        mapping_gz = assign_molecule_id(FROM_PDB_DIR / "148L.cif.gz")

        assert mapping_cif == mapping_gz


class TestAssignExtendedIds:
    """Tests for assign_extended_ids function."""

    def test_file_not_found(self) -> None:
        """Raises FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError):
            assign_extended_ids("nonexistent.cif")

    def test_directory_input(self, tmp_path: Path) -> None:
        """Raises ValueError for directory input."""
        with pytest.raises(ValueError, match="not a file"):
            assign_extended_ids(tmp_path)

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_returns_assignment_result(self) -> None:
        """Returns AssignmentResult with chain info."""
        result = assign_extended_ids(FROM_PDB_DIR / "148L.cif")

        assert hasattr(result, "chain_info")
        assert len(result.chain_info) > 0

        # Check ChainInfo structure
        for chain_id, info in result.chain_info.items():
            assert info.label_asym_id == chain_id
            assert isinstance(info.auth_asym_id, str)
            assert isinstance(info.entity_id, str)
            assert isinstance(info.entity_type, str)
            assert isinstance(info.molecule_id, int)
            assert isinstance(info.pn_unit_id, str)
            assert isinstance(info.pn_unit_entity, str)
            assert isinstance(info.molecule_entity, str)

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_molecule_id_mapping_property(self) -> None:
        """molecule_id_mapping property returns correct format."""
        result = assign_extended_ids(FROM_PDB_DIR / "148L.cif")

        mapping = result.molecule_id_mapping
        assert isinstance(mapping, dict)
        assert all(isinstance(v, int) for v in mapping.values())

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_writes_output_file(self, tmp_path: Path) -> None:
        """Writes output file when path is provided."""
        output_file = tmp_path / "output.cif"

        assign_extended_ids(FROM_PDB_DIR / "148L.cif", output_file)

        assert output_file.exists()
        assert output_file.stat().st_size > 0


class TestReferenceData:
    """Tests comparing results against reference data."""

    @pytest.mark.parametrize(
        "pdb_id,filename",
        [
            ("148L", "148L.cif"),
            ("148L", "148L.cif.gz"),
            ("1A0H", "1A0H.cif"),
            ("5FQD", "5FQD.cif.gz"),
        ],
    )
    def test_matches_reference(self, pdb_id: str, filename: str) -> None:
        """Results match AtomWorks reference data."""
        cif_path = FROM_PDB_DIR / filename
        ref_path = REF_DIR / f"{pdb_id}_ref.json"

        if not cif_path.exists():
            pytest.skip(f"Test data not available: {filename}")
        if not ref_path.exists():
            pytest.skip(f"Reference data not available: {pdb_id}_ref.json")

        with open(ref_path) as f:
            ref_data = json.load(f)

        result = assign_extended_ids(cif_path)

        # Check chain count
        assert len(result.chain_info) == ref_data["chain_count"]

        # Check molecule count
        molecule_ids = {info.molecule_id for info in result.chain_info.values()}
        assert len(molecule_ids) == ref_data["molecule_count"]

        # Check molecule_id mapping (must match AtomWorks exactly)
        for chain, mol_id in result.molecule_id_mapping.items():
            assert mol_id == ref_data["molecule_id_mapping"][chain], (
                f"Chain {chain}: got molecule_id={mol_id}, "
                f"expected {ref_data['molecule_id_mapping'][chain]}"
            )

        # Check pn_unit_id (AtomWorks reference uses same format)
        for chain, info in result.chain_info.items():
            ref_info = ref_data["chain_info"][chain]
            assert info.molecule_id == ref_info["molecule_id"]
            assert info.pn_unit_id == ref_info["pn_unit_id"]

    @pytest.mark.skipif(not _has_test_data("19HC.cif.gz"), reason="Test data not available")
    def test_19hc_grouping_consistency(self) -> None:
        """19HC: chains that should be grouped together are grouped together."""
        cif_path = FROM_PDB_DIR / "19HC.cif.gz"
        ref_path = REF_DIR / "19HC_ref.json"

        if not ref_path.exists():
            pytest.skip("Reference data not available")

        with open(ref_path) as f:
            ref_data = json.load(f)

        result = assign_extended_ids(cif_path)

        # Build groups from AtomWorks reference
        ref_groups: dict[int, set[str]] = {}
        for chain, mol_id in ref_data["molecule_id_mapping"].items():
            if mol_id not in ref_groups:
                ref_groups[mol_id] = set()
            ref_groups[mol_id].add(chain)

        # Check that chains in the same AtomWorks group have the same molecule_id in our result
        for _ref_mol_id, ref_chains in ref_groups.items():
            common_chains = [c for c in ref_chains if c in result.molecule_id_mapping]
            if len(common_chains) < 2:
                continue  # Can't verify grouping with less than 2 chains

            our_mol_ids = {result.molecule_id_mapping[c] for c in common_chains}
            assert len(our_mol_ids) == 1, (
                f"Chains {common_chains} should be in same molecule, "
                f"but got molecule_ids {our_mol_ids}"
            )


class TestSpecificStructures:
    """Tests for specific PDB structures with known characteristics."""

    @pytest.mark.skipif(not _has_test_data("5FQD.cif.gz"), reason="Test data not available")
    def test_5fqd_no_covalent_bonds(self) -> None:
        """5FQD has no covalent bonds - each chain is its own molecule."""
        result = assign_extended_ids(FROM_PDB_DIR / "5FQD.cif.gz")

        # Each chain should have a unique molecule_id
        molecule_ids = [info.molecule_id for info in result.chain_info.values()]
        assert len(molecule_ids) == len(set(molecule_ids)), "Each chain should be its own molecule"

        # Each chain's pn_unit should be itself (no grouping)
        for chain, info in result.chain_info.items():
            assert info.pn_unit_id == chain, f"Chain {chain} should be its own pn_unit"

    @pytest.mark.skipif(not _has_test_data("1A0H.cif"), reason="Test data not available")
    def test_1a0h_connectivity(self) -> None:
        """1A0H has chains connected via covale bonds (not disulf)."""
        result = assign_extended_ids(FROM_PDB_DIR / "1A0H.cif")

        # Disulfide bonds (A-B) are NOT counted as molecule connections
        # A is standalone (disulf to B doesn't connect them)
        # B, E, G are connected via covale bonds
        a_mol = result.chain_info["A"].molecule_id
        b_mol = result.chain_info["B"].molecule_id
        assert a_mol != b_mol, "A and B should be in different molecules (disulf ignored)"
        assert result.chain_info["E"].molecule_id == b_mol, "E should be in same molecule as B"
        assert result.chain_info["G"].molecule_id == b_mol, "G should be in same molecule as B"

        # Similar for the other complex: C is standalone, D,F,H are connected
        c_mol = result.chain_info["C"].molecule_id
        d_mol = result.chain_info["D"].molecule_id
        assert c_mol != d_mol, "C and D should be in different molecules (disulf ignored)"
        assert result.chain_info["F"].molecule_id == d_mol, "F should be in same molecule as D"
        assert result.chain_info["H"].molecule_id == d_mol, "H should be in same molecule as D"

    @pytest.mark.skipif(not _has_test_data("19HC.cif.gz"), reason="Test data not available")
    def test_19hc_many_ligands(self) -> None:
        """19HC has many ligands connected to polymer chains."""
        result = assign_extended_ids(FROM_PDB_DIR / "19HC.cif.gz")

        # Chain A (polymer) should have multiple ligands in the same molecule
        chain_a_mol = result.chain_info["A"].molecule_id
        mol_0_chains = [
            c for c, info in result.chain_info.items() if info.molecule_id == chain_a_mol
        ]
        assert len(mol_0_chains) > 1, "Chain A should share molecule with ligands"

    @pytest.mark.slow
    @pytest.mark.skipif(not _has_test_data("5Y6P.cif.gz"), reason="Test data not available")
    def test_5y6p_large_structure(self) -> None:
        """5Y6P is a large ribosome structure (performance test).

        This test verifies that the algorithm can handle large structures.
        No AtomWorks reference comparison is performed because:
        1. The reference file would be very large
        2. The primary purpose is to test performance, not correctness
        3. Correctness is validated by other structure tests
        """
        result = assign_extended_ids(FROM_PDB_DIR / "5Y6P.cif.gz")

        # Large structure with many chains (5Y6P has ~1100 chains)
        assert len(result.chain_info) > 1000


class TestSwapAuthAsymId:
    """Tests for swap_auth_asym_id function."""

    def test_valid_swap_targets_constant(self) -> None:
        """VALID_SWAP_TARGETS contains expected IDs."""
        assert "molecule_id" in VALID_SWAP_TARGETS
        assert "pn_unit_id" in VALID_SWAP_TARGETS
        assert "entity_id" in VALID_SWAP_TARGETS
        assert "label_asym_id" in VALID_SWAP_TARGETS

    def test_invalid_swap_target(self, tmp_path: Path) -> None:
        """Raises ValueError for invalid swap target."""
        if not (FROM_PDB_DIR / "148L.cif").exists():
            pytest.skip("Test data not available")

        with pytest.raises(ValueError, match="Invalid swap_with value"):
            swap_auth_asym_id(
                FROM_PDB_DIR / "148L.cif",
                tmp_path / "output.cif",
                swap_with="invalid_id",
            )

    def test_file_not_found(self, tmp_path: Path) -> None:
        """Raises FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError):
            swap_auth_asym_id(
                "nonexistent.cif",
                tmp_path / "output.cif",
            )

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_swap_with_molecule_id(self, tmp_path: Path) -> None:
        """Swaps auth_asym_id with molecule_id."""
        output_file = tmp_path / "output.cif"

        result = swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="molecule_id",
        )

        assert output_file.exists()
        assert len(result.chain_info) > 0

        # Verify the swap by reading the output
        doc = gemmi.cif.read(str(output_file))
        block = doc.sole_block()

        # Check that auth_asym_id contains molecule_id values (integer strings)
        auth_col = block.find_values("_atom_site.auth_asym_id")
        unique_auth = {str(v) for v in auth_col if str(v) not in ("?", ".")}
        assert all(v.isdigit() for v in unique_auth), f"Expected integers, got {unique_auth}"

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_swap_with_label_asym_id(self, tmp_path: Path) -> None:
        """Swaps auth_asym_id with label_asym_id."""
        output_file = tmp_path / "output.cif"

        swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="label_asym_id",
        )

        assert output_file.exists()

        # Verify auth_asym_id now contains label_asym_id values
        doc = gemmi.cif.read(str(output_file))
        block = doc.sole_block()

        auth_col = list(block.find_values("_atom_site.auth_asym_id"))
        label_col = list(block.find_values("_atom_site.label_asym_id"))
        assert auth_col == label_col

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_preserves_original_auth_asym_id(self, tmp_path: Path) -> None:
        """Original auth_asym_id is preserved in orig_auth_asym_id."""
        output_file = tmp_path / "output.cif"

        swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="molecule_id",
            preserve_original=True,
        )

        doc = gemmi.cif.read(str(output_file))
        block = doc.sole_block()

        # Check that orig_auth_asym_id column exists
        orig_col = block.find_values("_atom_site.orig_auth_asym_id")
        assert orig_col is not None, "orig_auth_asym_id column should exist"

        # Original values should include typical chain identifiers
        unique_orig = {str(v) for v in orig_col if str(v) not in ("?", ".")}
        assert len(unique_orig) > 0, "Should have preserved original auth values"

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_no_preserve_original(self, tmp_path: Path) -> None:
        """No orig_auth_asym_id when preserve_original=False."""
        output_file = tmp_path / "output.cif"

        swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="molecule_id",
            preserve_original=False,
        )

        doc = gemmi.cif.read(str(output_file))
        block = doc.sole_block()

        # orig_auth_asym_id column should not exist
        orig_col = block.find_values("_atom_site.orig_auth_asym_id")
        assert orig_col is None, "orig_auth_asym_id should not exist when preserve_original=False"

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_swap_with_pn_unit_id(self, tmp_path: Path) -> None:
        """Swaps auth_asym_id with pn_unit_id."""
        output_file = tmp_path / "output.cif"

        result = swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="pn_unit_id",
        )

        assert output_file.exists()
        assert len(result.chain_info) > 0

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_swap_with_entity_id(self, tmp_path: Path) -> None:
        """Swaps auth_asym_id with entity_id."""
        output_file = tmp_path / "output.cif"

        result = swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="entity_id",
        )

        assert output_file.exists()
        assert len(result.chain_info) > 0

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_returns_assignment_result(self, tmp_path: Path) -> None:
        """Returns AssignmentResult with correct structure."""
        output_file = tmp_path / "output.cif"

        result = swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
            swap_with="molecule_id",
        )

        # Check AssignmentResult structure
        assert hasattr(result, "chain_info")
        assert len(result.chain_info) > 0

        for chain_id, info in result.chain_info.items():
            assert info.label_asym_id == chain_id
            assert isinstance(info.molecule_id, int)
