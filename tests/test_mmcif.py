"""Tests for gemmi_extra_id.mmcif module."""

from pathlib import Path

import pytest

from gemmi_extra_id.mmcif import assign_molecule_id, assign_extended_ids

DATA_DIR = Path(__file__).parent.parent / "data"


class TestAssignMoleculeId:
    """Tests for assign_molecule_id function."""

    def test_file_not_found(self) -> None:
        """Raises FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError):
            assign_molecule_id("nonexistent.cif")

    def test_directory_input(self, tmp_path: Path) -> None:
        """Raises OSError for directory input."""
        with pytest.raises(OSError):
            assign_molecule_id(tmp_path)

    @pytest.mark.skipif(not (DATA_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_returns_ordered_dict(self) -> None:
        """Returns OrderedDict with correct structure."""
        mapping = assign_molecule_id(DATA_DIR / "148L.cif")

        assert isinstance(mapping, dict)
        assert len(mapping) > 0
        assert all(isinstance(k, str) for k in mapping.keys())
        assert all(isinstance(v, int) for v in mapping.values())

    @pytest.mark.skipif(not (DATA_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_writes_output_file(self, tmp_path: Path) -> None:
        """Writes output file when path is provided."""
        output_file = tmp_path / "output.cif"

        assign_molecule_id(DATA_DIR / "148L.cif", output_file)

        assert output_file.exists()
        assert output_file.stat().st_size > 0


class TestAssignExtendedIds:
    """Tests for assign_extended_ids function."""

    def test_file_not_found(self) -> None:
        """Raises FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError):
            assign_extended_ids("nonexistent.cif")

    @pytest.mark.skipif(not (DATA_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_returns_assignment_result(self) -> None:
        """Returns AssignmentResult with chain info."""
        result = assign_extended_ids(DATA_DIR / "148L.cif")

        assert hasattr(result, "chain_info")
        assert len(result.chain_info) > 0

        # Check ChainInfo structure
        for chain_id, info in result.chain_info.items():
            assert info.label_asym_id == chain_id
            assert isinstance(info.molecule_id, int)
            assert isinstance(info.pn_unit_id, str)
            assert isinstance(info.entity_type, str)

    @pytest.mark.skipif(not (DATA_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_molecule_id_mapping_property(self) -> None:
        """molecule_id_mapping property returns correct format."""
        result = assign_extended_ids(DATA_DIR / "148L.cif")

        mapping = result.molecule_id_mapping
        assert isinstance(mapping, dict)
        assert all(isinstance(v, int) for v in mapping.values())
