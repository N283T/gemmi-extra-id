"""Tests for gemmi_extra_id.mmcif module."""

from collections import OrderedDict
from pathlib import Path

import gemmi
import pytest

from gemmi_extra_id.mmcif import (
    VALID_SWAP_TARGETS,
    assign_molecule_id,
    swap_auth_asym_id,
)

DATA_DIR = Path(__file__).parent / "data"
FROM_PDB_DIR = DATA_DIR / "from_pdb"


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

        assert isinstance(mapping, OrderedDict)
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

        assert isinstance(mapping, OrderedDict)
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

    @pytest.mark.skipif(not _has_test_data("1A0H.cif"), reason="Test data not available")
    def test_1a0h_connectivity(self) -> None:
        """1A0H has chains connected via covale bonds (not disulf)."""
        mapping = assign_molecule_id(FROM_PDB_DIR / "1A0H.cif")

        # Disulfide bonds (A-B) are NOT counted as molecule connections
        # A is standalone (disulf to B doesn't connect them)
        # B, E, G are connected via covale bonds
        a_mol = mapping["A"]
        b_mol = mapping["B"]
        assert a_mol != b_mol, "A and B should be in different molecules (disulf ignored)"
        assert mapping["E"] == b_mol, "E should be in same molecule as B"
        assert mapping["G"] == b_mol, "G should be in same molecule as B"

        # Similar for the other complex: C is standalone, D,F,H are connected
        c_mol = mapping["C"]
        d_mol = mapping["D"]
        assert c_mol != d_mol, "C and D should be in different molecules (disulf ignored)"
        assert mapping["F"] == d_mol, "F should be in same molecule as D"
        assert mapping["H"] == d_mol, "H should be in same molecule as D"


class TestSwapAuthAsymId:
    """Tests for swap_auth_asym_id function."""

    def test_valid_swap_targets_constant(self) -> None:
        """VALID_SWAP_TARGETS contains expected IDs."""
        assert "molecule_id" in VALID_SWAP_TARGETS

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
        )

        assert output_file.exists()
        assert isinstance(result, OrderedDict)
        assert len(result) > 0

        # Verify the swap by reading the output
        doc = gemmi.cif.read(str(output_file))
        block = doc.sole_block()

        # Check that auth_asym_id contains molecule_id values (integer strings)
        auth_col = block.find_values("_atom_site.auth_asym_id")
        unique_auth = {str(v) for v in auth_col if str(v) not in ("?", ".")}
        assert all(v.isdigit() for v in unique_auth), f"Expected integers, got {unique_auth}"

    @pytest.mark.skipif(not _has_test_data("148L.cif"), reason="Test data not available")
    def test_preserves_original_auth_asym_id(self, tmp_path: Path) -> None:
        """Original auth_asym_id is preserved in orig_auth_asym_id."""
        output_file = tmp_path / "output.cif"

        swap_auth_asym_id(
            FROM_PDB_DIR / "148L.cif",
            output_file,
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
            preserve_original=False,
        )

        doc = gemmi.cif.read(str(output_file))
        block = doc.sole_block()

        # orig_auth_asym_id column should not exist in atom_site loop
        atom_site_loop = block.find_loop("_atom_site.label_asym_id").get_loop()
        tags = list(atom_site_loop.tags)
        assert "_atom_site.orig_auth_asym_id" not in tags, (
            "orig_auth_asym_id should not exist when preserve_original=False"
        )
