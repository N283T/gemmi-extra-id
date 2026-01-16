---
name: gemmi-improve
description: Improve gemmi-extra-id code using Gemmi documentation
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
timeout: 120000
---

# Gemmi Improve

Improve gemmi-extra-id code by referencing the local Gemmi documentation.

## Documentation Files

Read these files from `gemmi/docs/` as needed:

| File | Content |
|------|---------|
| `mol.rst` | Structure, Model, Chain, Residue, Atom hierarchy |
| `cif.rst` | CIF parsing, reading/writing mmCIF files |
| `chemistry.rst` | Chemical components, bonds, restraints |
| `analysis.rst` | Structure analysis functions |

## Source Files

Current implementation in `src/gemmi_extra_id/`:

| File | Purpose |
|------|---------|
| `mmcif.py` | mmCIF I/O with Gemmi |
| `graph.py` | Graph algorithms (find_components, find_pn_units) |
| `formatters.py` | Output formatters (JSON, CSV, table, tree) |
| `cli.py` | CLI interface |

## Workflow

1. **Ask what to improve**: If user didn't specify, ask what aspect they want to improve.

2. **Read relevant gemmi docs**:
   ```
   Read gemmi/docs/mol.rst   # For structure handling
   Read gemmi/docs/cif.rst   # For CIF I/O
   ```

3. **Review current implementation**: Read the relevant source file.

4. **Suggest improvements**: Based on documentation, suggest:
   - Better Gemmi API usage
   - Missing features Gemmi supports
   - Performance optimizations
   - Code simplifications

5. **Implement if approved**: Make the changes with proper testing.

## Common Improvement Areas

### Entity Handling
- `_entity` table - entity types and descriptions
- `_struct_asym` - asymmetric unit components
- `_entity_poly` - polymer sequences

### Connection Data
- `_struct_conn` - inter-residue connections
- Bond types: covale, disulf, hydrog, metalc
- Partner atoms and distances

### CIF Output
- Adding custom categories
- Writing loops efficiently
- Preserving original data

## Output Format

When suggesting improvements:

```markdown
## Current Implementation
[Description of current approach]

## Gemmi Documentation
[Relevant quotes from gemmi/docs]

## Suggested Improvement
[Proposed code changes]

## Benefits
- [Benefit 1]
- [Benefit 2]
```

## Example Session

User: "How can I improve entity extraction?"

1. Read `gemmi/docs/mol.rst` for entity-related sections
2. Read `src/gemmi_extra_id/mmcif.py` current implementation
3. Compare and suggest improvements
4. Implement if user approves
