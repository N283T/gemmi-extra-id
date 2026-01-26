# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - Unreleased

Major simplification release focusing on core molecule_id functionality.

### Removed

- Complete mode and `assign_extended_ids()` function
- Extended IDs: `chain_entity`, `pn_unit_id`, `pn_unit_entity`, `molecule_entity`
- Graph hashing functionality (`graph_hash.py`, `--hash-entities`)
- `find_pn_units()` function
- CSV, TSV, table, and tree output formats (only CIF and JSON remain)
- `--swap label_asym_id` option (only `molecule_id` remains)
- `--extended` CLI flag
- `--mode` CLI flag
- `hash` and `complete` optional dependencies

### Changed

- CLI is now an optional dependency (`pip install gemmi-extra-id[cli]`)
- `--swap` option now only accepts `molecule_id`
- Output formats reduced to `cif` (default) and `json`

## [0.1.0] - Unreleased

Initial release.

### Added

- Core functionality for assigning molecule IDs to mmCIF files based on covalent connectivity
- `assign_molecule_id()` - Assign molecule IDs using `_struct_conn` links
- `find_components()` - Find connected components in chain graphs
- CLI tool (`gemmi-extra-id assign`) with CIF and JSON output
- `swap_auth_asym_id()` - Swap `auth_asym_id` with `molecule_id`
- `--swap molecule_id` CLI option
- Support for Python 3.11, 3.12, and 3.13
- Type annotations throughout the codebase

### Changed

- Exclude disulfide bonds (`disulf`) from molecule_id calculation to align with AtomWorks (#5)

[Unreleased]: https://github.com/N283T/gemmi-extra-id/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/N283T/gemmi-extra-id/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/N283T/gemmi-extra-id/releases/tag/v0.1.0
