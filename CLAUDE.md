# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cifmolid is a Python library that assigns `molecule_id` to mmCIF files based on covalent connectivity from `_struct_conn` links. It identifies connected components of chains and writes the grouping back to the mmCIF file.

## Commands

```bash
# Install in development mode
uv pip install -e ".[dev]"

# Run CLI
uv run cifmolid assign data/148L.cif output.cif

# Run tests
uv run pytest

# Run tests on multiple Python versions (3.10, 3.11, 3.12)
uv run nox

# Format and lint
uv run ruff format src tests
uv run ruff check src tests
```

## Architecture

### Module Structure
- `src/cifmolid/core.py` - Graph algorithms (`find_components`)
- `src/cifmolid/cif.py` - mmCIF I/O with Gemmi (`assign_molecule_id`)
- `src/cifmolid/cli.py` - CLI entry point (typer/rich)
- `src/cifmolid/__init__.py` - Public API exports

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
from cifmolid import assign_molecule_id, find_components
```

## Terminology
- **label_asym_id**: Chain identifier in mmCIF (used for connectivity)
- **molecule_id**: Connected component ID (chains linked by covalent bonds share the same ID)
- **struct_conn**: mmCIF category storing inter-residue connections
