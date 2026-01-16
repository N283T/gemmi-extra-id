# Molecule ID Implementation Notes

## Overview
This repository adds `molecule_id` to mmCIF files using a covalent connectivity
approximation based on `_struct_conn` links. The goal is to match AtomWorks
molecule grouping without assembly expansion.

Key points:
- `molecule_id` is assigned per connected component of `label_asym_id` chains.
- Connectivity comes from `_struct_conn.conn_type_id` filtered by covalent types.
- Output is written back into the mmCIF `_atom_site` loop.
- A `_molecule_id_map` loop is added to make the chain mapping explicit.

## Scripts

### Add molecule_id using struct_conn
`./scripts/add_molecule_id_struct_conn.py`
- Builds a chain graph using `_struct_conn` and assigns component IDs.
- Default covalent types: `covale`, `disulf`.
- Writes `_atom_site.molecule_id` and `_molecule_id_map.*`.

Example:
```
uv run python scripts/add_molecule_id_struct_conn.py \
  data/148L_no_waters.cif \
  data/148L_struct_conn_molecule_id_no_waters.cif
```

### Generate AtomWorks reference
`./scripts/gen_atomworks_ref.py`
- Uses AtomWorks parser to write a reference mmCIF with `molecule_id`.

Example:
```
uv run python scripts/gen_atomworks_ref.py \
  data/148L.cif \
  data/148L_atomworks_reference.cif
```

### Generate AtomWorks reference (no waters)
`./scripts/gen_atomworks_molecule_id_no_waters.py`
- Same as above but targeted to `data/148L_no_waters.cif`.

Example:
```
uv run python scripts/gen_atomworks_molecule_id_no_waters.py
```

### Compare molecule_id between files
`./scripts/compare_molecule_id.py`
- Compares `_atom_site.molecule_id` across two files.
- Defaults to `_atom_site.id` as key if unique; otherwise uses a compound key.
- Prints summary stats and samples of mismatches.

Example:
```
uv run python scripts/compare_molecule_id.py \
  data/148L_atomworks_molecule_id_no_waters.cif \
  data/148L_struct_conn_molecule_id_no_waters.cif
```

## Data Layout
- `data/148L.cif`: original example input
- `data/148L_no_waters.cif`: water-removed input
- `data/148L_atomworks_reference.cif`: AtomWorks output from full input
- `data/148L_atomworks_molecule_id_no_waters.cif`: AtomWorks output from no-waters input
- `data/148L_struct_conn_molecule_id.cif`: struct_conn-based output from full input
- `data/148L_struct_conn_molecule_id_no_waters.cif`: struct_conn-based output from no-waters input
