"""Tests for gemmi_extra_id.cli module."""

import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

from gemmi_extra_id.cli import app

DATA_DIR = Path(__file__).parent / "data"
FROM_PDB_DIR = DATA_DIR / "from_pdb"
runner = CliRunner()


class TestCLI:
    """Tests for CLI commands."""

    def test_version(self) -> None:
        """--version prints version and exits."""
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "gemmi-extra-id" in result.output
        assert "0.1.0" in result.output

    def test_no_args_shows_help(self) -> None:
        """No arguments shows help."""
        result = runner.invoke(app, [])
        # typer exits with 2 when no command provided but shows help
        assert "Usage" in result.output
        assert "assign" in result.output

    def test_assign_missing_file(self) -> None:
        """assign with non-existent file returns error."""
        result = runner.invoke(app, ["assign", "nonexistent.cif"])
        assert result.exit_code == 2  # typer exits with 2 for invalid path

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_assign_with_output(self, tmp_path: Path) -> None:
        """assign creates output file."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.cif"

        result = runner.invoke(app, ["assign", str(input_file), str(output_file)])

        assert result.exit_code == 0
        assert output_file.exists()
        assert "Wrote" in result.output

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_assign_default_output(self, tmp_path: Path) -> None:
        """assign without output uses default naming."""
        input_file = tmp_path / "test.cif"
        shutil.copy(FROM_PDB_DIR / "148L.cif", input_file)

        result = runner.invoke(app, ["assign", str(input_file)])

        assert result.exit_code == 0
        expected_output = tmp_path / "test_molid.cif"
        assert expected_output.exists()

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_assign_custom_conn_types(self, tmp_path: Path) -> None:
        """assign with custom conn-types option."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.cif"

        result = runner.invoke(
            app, ["assign", str(input_file), str(output_file), "--conn-types", "covale"]
        )

        assert result.exit_code == 0
        assert output_file.exists()


class TestSwapOption:
    """Tests for --swap CLI option."""

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_swap_molecule_id(self, tmp_path: Path) -> None:
        """--swap molecule_id replaces auth_asym_id."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.cif"

        result = runner.invoke(
            app, ["assign", str(input_file), str(output_file), "--swap", "molecule_id"]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert output_file.exists()

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_swap_label_asym_id(self, tmp_path: Path) -> None:
        """--swap label_asym_id replaces auth_asym_id with label_asym_id."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.cif"

        result = runner.invoke(
            app, ["assign", str(input_file), str(output_file), "--swap", "label_asym_id"]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert output_file.exists()

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_swap_invalid_target(self, tmp_path: Path) -> None:
        """--swap with invalid target shows error."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.cif"

        result = runner.invoke(
            app, ["assign", str(input_file), str(output_file), "--swap", "invalid"]
        )

        assert result.exit_code != 0
        assert "Invalid" in result.output or "Error" in result.output

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_swap_incompatible_with_json_format(self, tmp_path: Path) -> None:
        """--swap cannot be used with non-CIF output format."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.json"

        result = runner.invoke(
            app,
            ["assign", str(input_file), str(output_file), "-f", "json", "--swap", "molecule_id"],
        )

        assert result.exit_code != 0

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_swap_incompatible_with_stdout(self) -> None:
        """--swap cannot be used with stdout output."""
        input_file = FROM_PDB_DIR / "148L.cif"

        result = runner.invoke(app, ["assign", str(input_file), "-", "--swap", "molecule_id"])

        assert result.exit_code != 0

    @pytest.mark.skipif(not (FROM_PDB_DIR / "148L.cif").exists(), reason="Test data not available")
    def test_swap_with_quiet(self, tmp_path: Path) -> None:
        """--swap works with --quiet option."""
        input_file = FROM_PDB_DIR / "148L.cif"
        output_file = tmp_path / "output.cif"

        result = runner.invoke(
            app, ["assign", str(input_file), str(output_file), "--swap", "molecule_id", "-q"]
        )

        assert result.exit_code == 0
        assert output_file.exists()
        assert "Wrote" not in result.output  # quiet mode
