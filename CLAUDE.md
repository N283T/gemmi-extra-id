# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

gemmi-extra-id is a Python library that assigns extra IDs (`molecule_id`, `pn_unit_id`, etc.) to mmCIF files based on covalent connectivity from `_struct_conn` links. It identifies connected components of chains and writes the grouping back to the mmCIF file.

## Commands

```bash
# Install in development mode
uv pip install -e ".[dev]"

# Run CLI
uv run gemmi-extra-id assign data/148L.cif output.cif

# Run tests
uv run pytest

# Run tests on multiple Python versions (3.11, 3.12, 3.13)
uv run nox

# Format and lint
uv run ruff format src tests
uv run ruff check src tests

# Type check
uv run ty check
```

## Architecture

### Module Structure
- `src/gemmi_extra_id/graph.py` - Graph algorithms (`find_components`, `find_pn_units`)
- `src/gemmi_extra_id/mmcif.py` - mmCIF I/O with Gemmi (`assign_molecule_id`, `assign_extended_ids`)
- `src/gemmi_extra_id/formatters.py` - Output formatters (JSON, CSV, TSV, table, tree)
- `src/gemmi_extra_id/cli.py` - CLI entry point (typer/rich)
- `src/gemmi_extra_id/__init__.py` - Public API exports

### Dev Scripts
- `scripts/atomworks_ref.py` - Generate AtomWorks reference (PEP 723, standalone)

### Processing Pipeline
1. Parse mmCIF with Gemmi
2. Build chain list from `_atom_site.label_asym_id`
3. Extract covalent edges from `_struct_conn` (types: `covale`, `disulf`)
4. Find connected components via iterative DFS
5. Write `_atom_site.molecule_id` and `_molecule_id_map` loop

### Public API
```python
from gemmi_extra_id import assign_molecule_id, assign_extended_ids, find_components
```

## Terminology
- **label_asym_id**: Chain identifier in mmCIF (used for connectivity)
- **molecule_id**: Connected component ID (chains linked by covalent bonds share the same ID)
- **pn_unit_id**: Same-type chains covalently linked (comma-separated chain list)
- **entity_id**: Entity ID from `_struct_asym.entity_id`
- **struct_conn**: mmCIF category storing inter-residue connections
