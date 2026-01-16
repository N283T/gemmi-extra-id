"""Output formatters for molecule_id mappings."""

from __future__ import annotations

import csv
import json
import sys
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections import OrderedDict


def to_json(mapping: OrderedDict[str, int], indent: int = 2) -> str:
    """Convert mapping to JSON string."""
    return json.dumps(dict(mapping), indent=indent)


def to_csv(mapping: OrderedDict[str, int]) -> str:
    """Convert mapping to CSV string."""
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(["label_asym_id", "molecule_id"])
    for chain, mol_id in mapping.items():
        writer.writerow([chain, mol_id])
    return output.getvalue()


def to_tsv(mapping: OrderedDict[str, int]) -> str:
    """Convert mapping to TSV string."""
    output = StringIO()
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(["label_asym_id", "molecule_id"])
    for chain, mol_id in mapping.items():
        writer.writerow([chain, mol_id])
    return output.getvalue()


def to_table(mapping: OrderedDict[str, int]) -> str:
    """Convert mapping to human-readable table string."""
    lines = ["label_asym_id  molecule_id", "-" * 26]
    for chain, mol_id in mapping.items():
        lines.append(f"{chain:<14} {mol_id}")
    return "\n".join(lines)


def write_output(
    mapping: OrderedDict[str, int],
    output_path: str | Path | None,
    fmt: str = "table",
) -> None:
    """
    Write mapping to file or stdout in specified format.

    Args:
        mapping: The molecule_id mapping to write.
        output_path: Path to write to, or "-" for stdout, or None for no output.
        fmt: Output format ("json", "csv", "tsv", "table").
    """
    if output_path is None:
        return

    formatters = {
        "json": to_json,
        "csv": to_csv,
        "tsv": to_tsv,
        "table": to_table,
    }

    if fmt not in formatters:
        raise ValueError(f"Unknown format: {fmt}. Use one of: {', '.join(formatters)}")

    content = formatters[fmt](mapping)

    if str(output_path) == "-":
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")
    else:
        Path(output_path).write_text(content)
