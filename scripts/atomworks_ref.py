#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "atomworks>=2.2.0",
# ]
# ///
"""
Generate AtomWorks reference mmCIF with molecule_id.

This script uses AtomWorks to parse an mmCIF file and write a reference
output with molecule_id assignments for comparison with molid.

Usage:
    uv run scripts/atomworks_ref.py input.cif [output.cif]
"""

from __future__ import annotations

import sys
from pathlib import Path

from atomworks.io import parse_mmcif, write_mmcif


def main() -> int:
    """Generate AtomWorks reference mmCIF."""
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.cif> [output.cif]", file=sys.stderr)
        return 1

    input_path = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        output_path = Path(sys.argv[2])
    else:
        output_path = input_path.with_stem(f"{input_path.stem}_atomworks_ref")

    if not input_path.exists():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        return 1

    print(f"Parsing {input_path} with AtomWorks...")
    atom_array = parse_mmcif(input_path)

    print(f"Writing {output_path}...")
    write_mmcif(atom_array, output_path)

    # Print molecule_id summary
    if hasattr(atom_array, "molecule_id"):
        mol_ids = atom_array.molecule_id
        unique_ids = set(mol_ids)
        print(f"Found {len(unique_ids)} unique molecule_id(s)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
