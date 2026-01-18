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
# Test molecule_id and pn_unit_id
uv run python scripts/test_all_pdb.py subset tests/data/failure_cases -c both

# Test entity equivalence relations only
uv run python scripts/test_all_pdb.py subset tests/data/failure_cases -c entity

# Test all (molecule_id + pn_unit_id + entity)
uv run python scripts/test_all_pdb.py subset tests/data/failure_cases -c all

# Test with 100 random entries from PDB mirror
uv run python scripts/test_all_pdb.py subset tests/data/failure_cases -c all -r 100
```

### Comparison Modes

| Mode | Description |
|------|-------------|
| `molecule` | Compare molecule_id grouping |
| `pn_unit` | Compare pn_unit_id mapping |
| `both` | Compare molecule_id + pn_unit_id |
| `entity` | Compare entity equivalence relations (chain_entity, pn_unit_entity, molecule_entity) |
| `all` | Compare molecule_id + pn_unit_id + entity |

## Entity Field Comparison

### Why Values Differ

AtomWorks and gemmi-extra-id produce **different values** for entity fields but with **equivalent semantic meaning**:

| Field | AtomWorks | gemmi-extra-id |
|-------|-----------|----------------|
| chain_entity | Graph hash of residue structure (0-indexed) | CIF `_struct_asym.entity_id` |
| pn_unit_entity | Graph hash of chain_entity composition | Minimum entity_id in pn_unit |
| molecule_entity | Graph hash of pn_unit_entity composition | Minimum entity_id in molecule |

**AtomWorks approach**:
- Computes a hash of the graph structure at each level
- Chains with identical residue sequences get the same `chain_entity`
- Uses 0-indexed integers assigned in order of first occurrence

**gemmi-extra-id approach**:
- Uses PDB's official `entity_id` from CIF file
- Chains with the same `entity_id` belong to the same molecular entity
- Uses 1-indexed integers from CIF

### Equivalence Relation

Both approaches produce the **same groupings** (equivalence relations):

```
AtomWorks:  {A: 0, B: 0, C: 1}  →  equivalence groups: {{A, B}, {C}}
gemmi:      {A: 1, B: 1, C: 2}  →  equivalence groups: {{A, B}, {C}}
```

The test compares these equivalence groups, not the raw values.

## File List

See `failure_ids.txt` for complete list of 72 PDB IDs.
