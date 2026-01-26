#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "atomworks>=2.2.0",
# ]
# ///
"""
Generate AtomWorks reference data for testing.

This script uses AtomWorks to parse an mmCIF file and extract the
molecule_id mapping for comparison with gemmi-extra-id.

Usage:
    uv run scripts/atomworks_ref.py input.cif [output.json]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> int:
    """Generate AtomWorks reference data."""
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.cif> [output.json]", file=sys.stderr)
        return 1

    input_path = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        output_path = Path(sys.argv[2])
    else:
        output_path = input_path.with_stem(f"{input_path.stem}_ref").with_suffix(".json")

    if not input_path.exists():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        return 1

    # Import here to avoid slow import on --help
    from atomworks.io import parse

    print(f"Parsing {input_path} with AtomWorks...")
    result = parse(
        filename=str(input_path),
        build_assembly=None,  # ASU only
        add_missing_atoms=False,
        remove_waters=False,
    )

    asym = result["asym_unit"][0]

    # Extract unique chain -> molecule_id mapping
    chain_ids = asym.chain_id
    molecule_ids = asym.molecule_id

    # Build mapping (first occurrence of each chain)
    mapping: dict[str, int] = {}
    chain_order: list[str] = []
    for chain, mol_id in zip(chain_ids, molecule_ids, strict=True):
        if chain not in mapping:
            mapping[chain] = int(mol_id)
            chain_order.append(chain)

    # Also extract other IDs if available
    categories = asym.get_annotation_categories()
    chain_info: dict[str, dict] = {}

    for chain in chain_order:
        # Get first atom of this chain
        mask = asym.chain_id == chain
        first_idx = mask.argmax()

        info: dict[str, str | int] = {
            "molecule_id": int(asym.molecule_id[first_idx]),
        }

        if "chain_entity" in categories:
            info["chain_entity"] = int(asym.chain_entity[first_idx])
        if "pn_unit_id" in categories:
            info["pn_unit_id"] = str(asym.pn_unit_id[first_idx])
        if "pn_unit_entity" in categories:
            info["pn_unit_entity"] = int(asym.pn_unit_entity[first_idx])
        if "molecule_entity" in categories:
            info["molecule_entity"] = int(asym.molecule_entity[first_idx])
        if "chain_type" in categories:
            info["chain_type"] = str(asym.chain_type[first_idx])
        if "label_entity_id" in categories:
            info["label_entity_id"] = str(asym.label_entity_id[first_idx])

        chain_info[chain] = info

    # Build output
    output_data = {
        "source": "atomworks",
        "input_file": input_path.name,
        "chain_count": len(chain_order),
        "molecule_count": len(set(mapping.values())),
        "molecule_id_mapping": mapping,
        "chain_info": chain_info,
    }

    print(f"Writing {output_path}...")
    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)

    print(f"Found {len(chain_order)} chains, {len(set(mapping.values()))} molecules")
    return 0


if __name__ == "__main__":
    sys.exit(main())
