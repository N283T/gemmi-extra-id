# gemmi-extra-id

Assign extra IDs to mmCIF files based on covalent connectivity.

This is a [Gemmi](https://gemmi.readthedocs.io/)-based reimplementation of the ID system defined in the [AtomWorks Glossary](https://rosettacommons.github.io/atomworks/latest/glossary.html). It identifies connected components of chains using covalent bond information from `_struct_conn` and writes the grouping back to the mmCIF file.

## Installation

```bash
pip install gemmi-extra-id
```

Or with uv:

```bash
uv add gemmi-extra-id
```

## Usage

### Command Line

```bash
# Assign molecule_id to an mmCIF file
gemmi-extra-id assign input.cif output.cif

# Use default output name (input_molid.cif)
gemmi-extra-id assign input.cif

# Output as JSON
gemmi-extra-id assign input.cif -f json output.json

# Output to stdout (for piping)
gemmi-extra-id assign input.cif -f json -

# Available formats: cif, json, csv, tsv, table, tree
gemmi-extra-id assign input.cif -f csv output.csv

# Specify custom covalent bond types
gemmi-extra-id assign input.cif --conn-types covale,disulf,metalc
```

### Extended Mode

Use `--extended` flag to include all IDs (chain_entity, pn_unit_id, pn_unit_entity, molecule_entity):

```bash
# Extended output as table
gemmi-extra-id assign input.cif -f table --extended

# Tree view shows hierarchical structure
gemmi-extra-id assign input.cif -f tree
```

Tree output displays the ID hierarchy:

```
148L (6 chains)
├── molecule_id=0 (entity=1)
│   ├── pn_unit A (entity=1)
│   │   └── A (polymer, entity=1)
│   ├── pn_unit B (entity=2)
│   │   └── B (polymer, entity=2)
│   └── pn_unit C (entity=3)
│       └── C (branched, entity=3)
├── molecule_id=1 (entity=4)
│   └── pn_unit D (entity=4)
│       └── D (non-polymer, entity=4)
...
```

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
mapping = assign_molecule_id("input.cif", "output.cif", covalent_types={"covale", "disulf", "metalc"})
```

### Extended API

```python
from gemmi_extra_id import assign_extended_ids

# Get all IDs
result = assign_extended_ids("input.cif", "output.cif")

for chain, info in result.chain_info.items():
    print(f"{chain}: molecule={info.molecule_id}, pn_unit={info.pn_unit_id}")
```

### Output Formatters

```python
from gemmi_extra_id import assign_molecule_id, to_json, to_csv

mapping = assign_molecule_id("input.cif")

# Convert to different formats
print(to_json(mapping))   # {"A": 0, "B": 0, "C": 1}
print(to_csv(mapping))    # label_asym_id,molecule_id\nA,0\n...
```

## Output

gemmi-extra-id adds two elements to the output mmCIF:

1. `_atom_site.molecule_id` column - component ID for each atom
2. `_molecule_id_map` loop - mapping of `label_asym_id` to `molecule_id`

With `--extended`, additional IDs are computed:

| ID | Description |
|----|-------------|
| `chain_entity` | Entity ID from `_struct_asym.entity_id` |
| `pn_unit_id` | Same-type chains covalently linked |
| `pn_unit_entity` | Entity ID for the pn_unit |
| `molecule_entity` | Entity ID for the molecule (minimum entity_id when mixed) |

## AtomWorks Compatibility

By default, gemmi-extra-id uses only `covale` (covalent) bond types for molecule grouping, matching [AtomWorks](https://github.com/rosettacommons/atomworks) behavior. Disulfide bonds (`disulf`) are **not** considered as intermolecular connections.

This means chains connected only by disulfide bonds will have different `molecule_id` values:

```python
# 1A0H: Chain A and B are connected by disulfide bond
# But they are assigned different molecule_ids
result = assign_extended_ids("1A0H.cif")
assert result.chain_info["A"].molecule_id != result.chain_info["B"].molecule_id
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
