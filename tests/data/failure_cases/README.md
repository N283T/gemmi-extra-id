# Historical Failure Cases

This directory contains PDB entries that previously failed molecule_id or pn_unit_id comparison with AtomWorks.

## Summary

| Category | Count | Root Cause | Fix |
|----------|-------|------------|-----|
| molecule_id | 43 | Various (struct_conn handling, etc.) | Earlier PRs |
| pn_unit_id (entity_type grouping) | 318 | Non-polymers grouped by entity_type | PR #12 |
| pn_unit_id (single-entity format) | 29 | CIF single-value format not parsed | PR #13 |
| **Total unique** | **72** | | |

## IGNORE_LIST (3 entries)

These entries have known PDB data quality issues (orphan chains not in `_struct_asym`):

| PDB ID | Issue |
|--------|-------|
| 2g10 | Chain F: 147 water atoms not in _struct_asym |
| 1ts6 | Chain C: 7 water atoms not in _struct_asym |
| 2k9y | Chains C, D: not in _struct_asym |

## pn_unit_id Failures (Fixed in PR #12)

**Root Cause**: `find_pn_units()` grouped non-polymer chains by `entity_type` first, then found connected components within each type. This broke connectivity when chains with different entity_types (e.g., branched vs non-polymer) were covalently linked.

**Example (146d)**: Chains C (branched), D (branched), H (non-polymer) connected via covale bonds were incorrectly separated.

## pn_unit_id Failures (Fixed in PR #13)

**Root Cause**: CIF files with single entity use non-loop format:
```
_entity.id   1
_entity.type polymer
```

`find_loop()` returns a nil Column object (not Python `None`), so the check `if entity_col is not None` passed but the loop was empty. This caused `entity_type` to default to "unknown", treating polymers as non-polymers.

**Affected entries**: 29 PDB entries with single entity (all chains are the same polymer type).

## Usage

Run regression tests:
```bash
# Test failure cases only
uv run python scripts/test_all_pdb.py subset tests/data/failure_cases -c both

# Test failure cases + 100 random entries
uv run python scripts/test_all_pdb.py subset tests/data/failure_cases -c both -r 100
```

## File List

See `failure_ids.txt` for complete list of 72 PDB IDs.
