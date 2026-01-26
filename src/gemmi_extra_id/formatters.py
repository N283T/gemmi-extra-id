"""Output formatters for molecule_id mappings."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections import OrderedDict


def to_json(mapping: OrderedDict[str, int], indent: int = 2) -> str:
    """Convert mapping to JSON string."""
    return json.dumps(dict(mapping), indent=indent)


def write_output(
    mapping: OrderedDict[str, int],
    output_path: str | Path | None,
    fmt: str = "json",
) -> None:
    """
    Write mapping to file or stdout in specified format.

    Args:
        mapping: The molecule_id mapping to write.
        output_path: Path to write to, or "-" for stdout, or None for no output.
        fmt: Output format ("json").
    """
    if output_path is None:
        return

    if fmt != "json":
        raise ValueError(f"Unknown format: {fmt}. Use 'json'")

    content = to_json(mapping)

    if str(output_path) == "-":
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")
    else:
        Path(output_path).write_text(content)
