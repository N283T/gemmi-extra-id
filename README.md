# gemmi-extra-id

Assign molecule_id to mmCIF files based on covalent connectivity.

This is a [Gemmi](https://gemmi.readthedocs.io/)-based tool that identifies connected components of chains using covalent bond information from `_struct_conn` and writes the grouping back to the mmCIF file.

## Installation

### Library only

```bash
pip install gemmi-extra-id
```

Or with uv:

```bash
uv add gemmi-extra-id
```

### With CLI

```bash
pip install gemmi-extra-id[cli]
```

Or with uv:

```bash
uv add gemmi-extra-id[cli]
```

## Usage

### Command Line

> **Note:** CLI requires installation with `[cli]` extra.

```bash
# Assign molecule_id to an mmCIF file
gemmi-extra-id assign input.cif output.cif

# Use default output name (input_molid.cif)
gemmi-extra-id assign input.cif

# Output as JSON
gemmi-extra-id assign input.cif -f json output.json

# Output to stdout (for piping)
gemmi-extra-id assign input.cif -f json -

# Specify custom covalent bond types
gemmi-extra-id assign input.cif --conn-types covale,disulf,metalc
```

### Swap Mode

Use `--swap` to replace `auth_asym_id` with `molecule_id`. This is useful for applications that only read `auth_asym_id`:

```bash
# Replace auth_asym_id with molecule_id
gemmi-extra-id assign input.cif output.cif --swap molecule_id
```

The original `auth_asym_id` values are preserved in `_atom_site.orig_auth_asym_id`.

### Python API

```python
from gemmi_extra_id import assign_molecule_id

# Assign and write to file
mapping = assign_molecule_id("input.cif", "output.cif")
print(mapping)
# OrderedDict([('A', 0), ('B', 0), ('C', 1), ...])

# Get mapping without writing
mapping = assign_molecule_id("input.cif")

# Custom covalent types
mapping = assign_molecule_id("input.cif", "output.cif", covalent_types={"covale", "disulf"})
```

### Swap API

```python
from gemmi_extra_id import swap_auth_asym_id

# Swap auth_asym_id with molecule_id
mapping = swap_auth_asym_id("input.cif", "output.cif")

# Don't preserve original auth_asym_id
mapping = swap_auth_asym_id("input.cif", "output.cif", preserve_original=False)
```

### Output Formatters

```python
from gemmi_extra_id import assign_molecule_id, to_json

mapping = assign_molecule_id("input.cif")

# Convert to JSON
print(to_json(mapping))   # {"A": 0, "B": 0, "C": 1}
```

## Output

gemmi-extra-id adds two elements to the output mmCIF:

1. `_atom_site.molecule_id` column - component ID for each atom
2. `_molecule_id_map` loop - mapping of `label_asym_id` to `molecule_id`

## Covalent Bond Types

By default, gemmi-extra-id uses only `covale` (covalent) bond types for molecule grouping. Disulfide bonds (`disulf`) are **not** considered as intermolecular connections.

This means chains connected only by disulfide bonds will have different `molecule_id` values:

```python
# 1A0H: Chain A and B are connected by disulfide bond
# But they are assigned different molecule_ids
mapping = assign_molecule_id("1A0H.cif")
assert mapping["A"] != mapping["B"]
```

To include disulfide bonds in molecule grouping, specify `--conn-types`:

```bash
gemmi-extra-id assign input.cif --conn-types covale,disulf
```

Or in Python:

```python
mapping = assign_molecule_id("input.cif", covalent_types={"covale", "disulf"})
```

## Development

```bash
# Install with CLI and dev dependencies
uv pip install -e ".[cli]" --group dev

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
