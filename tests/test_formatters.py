"""Tests for cifmolid.output module."""

import json
from collections import OrderedDict

from cifmolid.formatters import to_csv, to_json, to_table, to_tsv


class TestOutputFormats:
    """Tests for output format functions."""

    def test_to_json(self) -> None:
        """to_json produces valid JSON."""
        mapping = OrderedDict([("A", 0), ("B", 0), ("C", 1)])
        result = to_json(mapping)
        parsed = json.loads(result)
        assert parsed == {"A": 0, "B": 0, "C": 1}

    def test_to_csv(self) -> None:
        """to_csv produces valid CSV."""
        mapping = OrderedDict([("A", 0), ("B", 1)])
        result = to_csv(mapping)
        lines = [line.strip() for line in result.strip().split("\n")]
        assert lines[0] == "label_asym_id,molecule_id"
        assert lines[1] == "A,0"
        assert lines[2] == "B,1"

    def test_to_tsv(self) -> None:
        """to_tsv produces valid TSV."""
        mapping = OrderedDict([("A", 0), ("B", 1)])
        result = to_tsv(mapping)
        lines = [line.strip() for line in result.strip().split("\n")]
        assert lines[0] == "label_asym_id\tmolecule_id"
        assert lines[1] == "A\t0"
        assert lines[2] == "B\t1"

    def test_to_table(self) -> None:
        """to_table produces human-readable table."""
        mapping = OrderedDict([("A", 0), ("B", 1)])
        result = to_table(mapping)
        assert "label_asym_id" in result
        assert "molecule_id" in result
        assert "A" in result
        assert "B" in result

    def test_empty_mapping(self) -> None:
        """Empty mapping produces valid output."""
        mapping: OrderedDict[str, int] = OrderedDict()
        assert to_json(mapping) == "{}"
        assert "label_asym_id,molecule_id" in to_csv(mapping)
