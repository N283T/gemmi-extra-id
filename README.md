# gemmi-extra-id

[![CI](https://github.com/N283T/gemmi-extra-id/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/gemmi-extra-id/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Assign `molecule_id` to mmCIF files based on covalent bond connectivity.

## What is molecule_id?

In mmCIF files, chains (identified by `label_asym_id`) can be covalently linked together. `molecule_id` groups these connected chains into logical molecules.

**Example:** A protein-ligand complex where chain A (protein) and chain B (covalently attached ligand) would share the same `molecule_id`, while chain C (a separate small molecule) would have a different `molecule_id`.

```
Chain A (protein) ─── covale ─── Chain B (ligand)  →  molecule_id = 0
Chain C (separate molecule)                        →  molecule_id = 1
```

## Installation

```bash
# Library only
pip install gemmi-extra-id

# With CLI
pip install gemmi-extra-id[cli]
```

## Quick Start

### Command Line

```bash
# Basic usage - adds molecule_id to CIF
gemmi-extra-id assign input.cif output.cif

# Output as JSON mapping
gemmi-extra-id assign input.cif -f json -
# {"A": 0, "B": 0, "C": 1}
```

### Python API

```python
from gemmi_extra_id import assign_molecule_id

mapping = assign_molecule_id("input.cif", "output.cif")
print(mapping)
# OrderedDict([('A', 0), ('B', 0), ('C', 1)])
```

## Output Format

The tool adds two elements to the mmCIF file:

1. **`_atom_site.molecule_id`** - Component ID for each atom
2. **`_molecule_id_map`** - Summary mapping table

```
loop_
_molecule_id_map.label_asym_id
_molecule_id_map.auth_asym_id
_molecule_id_map.molecule_id
A A 0
B B 0
C C 1
```

## CLI Options

```bash
gemmi-extra-id assign [OPTIONS] INPUT_FILE [OUTPUT_FILE]

Options:
  -f, --format [cif|json]   Output format (default: cif)
  -c, --conn-types TEXT     Bond types to consider (default: covale)
  -s, --swap molecule_id    Replace auth_asym_id with molecule_id
  -q, --quiet               Suppress status messages
```

### Swap Mode

Replace `auth_asym_id` with `molecule_id` for applications that only read `auth_asym_id`:

```bash
gemmi-extra-id assign input.cif output.cif --swap molecule_id
```

Original values are preserved in `_atom_site.orig_auth_asym_id`.

## Bond Type Configuration

By default, only `covale` (covalent) bonds are considered. Disulfide bonds (`disulf`) are excluded.

```bash
# Include disulfide bonds
gemmi-extra-id assign input.cif output.cif --conn-types covale,disulf
```

```python
mapping = assign_molecule_id("input.cif", covalent_types={"covale", "disulf"})
```

## API Reference

### assign_molecule_id

```python
def assign_molecule_id(
    input_path: str | Path,
    output_path: str | Path | None = None,
    covalent_types: set[str] | None = None,  # default: {"covale"}
) -> OrderedDict[str, int]:
    """Assign molecule_id based on covalent connectivity."""
```

### swap_auth_asym_id

```python
def swap_auth_asym_id(
    input_path: str | Path,
    output_path: str | Path,
    preserve_original: bool = True,
    covalent_types: set[str] | None = None,
) -> OrderedDict[str, int]:
    """Replace auth_asym_id with molecule_id."""
```

## Development

```bash
# Clone and install
git clone https://github.com/N283T/gemmi-extra-id.git
cd gemmi-extra-id
uv sync --group dev

# Run tests
uv run pytest

# Lint
uv run ruff check src tests
uv run ruff format src tests
```

## License

MIT
