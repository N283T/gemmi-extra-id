#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "atomworks>=2.2.0",
#     "gemmi>=0.7.0",
#     "gemmi-extra-id",
# ]
# ///
"""
Test gemmi-extra-id against AtomWorks on all PDB files.

Usage:
    uv run scripts/test_all_pdb.py /path/to/pdb/mirror [--output results.json]
"""

import json
import sys
from collections.abc import Iterator
from pathlib import Path

import gemmi


def find_cif_files(mirror_path: Path) -> Iterator[Path]:
    """Find all CIF files in PDB mirror directory using gemmi.CifWalk."""
    for cif_str in gemmi.CifWalk(str(mirror_path)):
        yield Path(cif_str)


def extract_atomworks_ids(cif_path: Path) -> dict:
    """Extract molecule_id mapping from AtomWorks parse result."""
    from atomworks.io import parse

    result = parse(
        filename=cif_path,
        build_assembly=None,
        add_missing_atoms=False,
        remove_waters=False,
    )
    asym = result["asym_unit"][0]

    # Build chain -> molecule_id mapping
    mapping = {}
    for chain, mol_id in zip(asym.chain_id, asym.molecule_id, strict=True):
        if chain not in mapping:
            mapping[chain] = int(mol_id)

    return {"molecule_id_mapping": mapping}


def extract_gemmi_ids(cif_path: Path) -> dict:
    """Extract molecule_id mapping from gemmi-extra-id."""
    from gemmi_extra_id import assign_extended_ids

    result = assign_extended_ids(cif_path)
    return {"molecule_id_mapping": result.molecule_id_mapping}


def compare_results(expected: dict, actual: dict) -> tuple[bool, list[str]]:
    """Compare AtomWorks and gemmi-extra-id results."""
    errors = []

    exp_mapping = expected["molecule_id_mapping"]
    act_mapping = actual["molecule_id_mapping"]

    # Check all chains match
    if set(exp_mapping.keys()) != set(act_mapping.keys()):
        errors.append(
            f"Chain mismatch: expected {set(exp_mapping.keys())}, got {set(act_mapping.keys())}"
        )

    # Check molecule_id values
    for chain in exp_mapping:
        if chain in act_mapping and exp_mapping[chain] != act_mapping[chain]:
            errors.append(
                f"Chain {chain}: expected mol_id={exp_mapping[chain]}, got {act_mapping[chain]}"
            )

    return len(errors) == 0, errors


def main() -> int:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <pdb_mirror_path> [--output results.json]")
        return 1

    mirror_path = Path(sys.argv[1])
    output_path = None
    if "--output" in sys.argv:
        idx = sys.argv.index("--output")
        output_path = Path(sys.argv[idx + 1])

    results: dict[str, list] = {"passed": [], "failed": [], "errors": []}

    for cif_path in find_cif_files(mirror_path):
        pdb_id = cif_path.stem.replace(".cif", "")
        try:
            expected = extract_atomworks_ids(cif_path)
            actual = extract_gemmi_ids(cif_path)
            passed, errors = compare_results(expected, actual)

            if passed:
                results["passed"].append(pdb_id)
                print(f"PASS: {pdb_id}")
            else:
                results["failed"].append({"pdb_id": pdb_id, "errors": errors})
                print(f"FAIL: {pdb_id} - {errors}")
        except Exception as e:
            results["errors"].append({"pdb_id": pdb_id, "error": str(e)})
            print(f"ERROR: {pdb_id} - {e}")

    # Summary
    print(
        f"\nSummary: {len(results['passed'])} passed, "
        f"{len(results['failed'])} failed, {len(results['errors'])} errors"
    )

    if output_path:
        with open(output_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Results written to {output_path}")

    return 0 if not results["failed"] else 1


if __name__ == "__main__":
    sys.exit(main())
