# molecule-id-gemmi

Utilities for assigning `molecule_id` to mmCIF files using Gemmi, with
AtomWorks references for comparison.

## Layout
- `scripts/`: executable helpers
- `data/`: example inputs/outputs
- `docs/IMPLEMENTATION.md`: implementation notes and usage

## Quick start
```
uv run python scripts/add_molecule_id_struct_conn.py
uv run python scripts/compare_molecule_id.py \
  data/148L_atomworks_molecule_id_no_waters.cif \
  data/148L_struct_conn_molecule_id_no_waters.cif
```
