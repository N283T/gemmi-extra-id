# cifmolid

Assign `molecule_id` to mmCIF files based on covalent connectivity.

cifmolid identifies connected components of chains using covalent bond information from `_struct_conn` and writes the grouping back to the mmCIF file.

## Installation

```bash
pip install cifmolid
```

Or with uv:

```bash
uv add cifmolid
```

## Usage

### Command Line

```bash
# Assign molecule_id to an mmCIF file
cifmolid assign input.cif output.cif

# Use default output name (input_molid.cif)
cifmolid assign input.cif

# Specify custom covalent bond types
cifmolid assign input.cif --conn-types covale,disulf,metalc
```

### Python API

```python
from cifmolid import assign_molecule_id

# Assign and write to file
mapping = assign_molecule_id("input.cif", "output.cif")
print(mapping)
# OrderedDict([('A', 0), ('B', 0), ('C', 1), ...])

# Get mapping without writing
mapping = assign_molecule_id("input.cif")

# Custom covalent types
mapping = assign_molecule_id("input.cif", "output.cif", covalent_types={"covale", "disulf", "metalc"})
```

### Lower-level API

```python
from cifmolid import find_components

# Find connected components from nodes and edges
nodes = ["A", "B", "C", "D"]
edges = [("A", "B"), ("C", "D")]
mapping = find_components(nodes, edges)
# OrderedDict([('A', 0), ('B', 0), ('C', 1), ('D', 1)])
```

## Output

cifmolid adds two elements to the output mmCIF:

1. `_atom_site.molecule_id` column - component ID for each atom
2. `_molecule_id_map` loop - mapping of `label_asym_id` to `molecule_id`

## Development

```bash
# Install with dev dependencies
uv pip install -e ".[dev]"

# Run tests
uv run pytest

# Run tests on multiple Python versions
uv run nox

# Lint and format
uv run ruff check src tests
uv run ruff format src tests
```

## License

MIT
