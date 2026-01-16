# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- CLI is now an optional dependency (`pip install gemmi-extra-id[cli]`)

## [0.1.0] - Unreleased

Initial release.

### Added

- Core functionality for assigning molecule IDs to mmCIF files based on covalent connectivity
- `assign_molecule_id()` - Assign molecule IDs using `_struct_conn` links
- `assign_extended_ids()` - Assign extended IDs (molecule_id, pn_unit_id, entity_id)
- `find_components()` - Find connected components in chain graphs
- `find_pn_units()` - Group chains by entity into PN units
- CLI tool (`gemmi-extra-id assign`) with multiple output formats
- Output formatters: JSON, CSV, TSV, table, tree
- `swap_auth_asym_id()` - Swap `auth_asym_id` with assigned IDs for legacy application compatibility
- `--swap` CLI option to replace `auth_asym_id` with `molecule_id`, `pn_unit_id`, `entity_id`, or `label_asym_id`
- Support for Python 3.11, 3.12, and 3.13
- Type annotations throughout the codebase

### Changed

- Exclude disulfide bonds (`disulf`) from molecule_id calculation to align with AtomWorks (#5)

[Unreleased]: https://github.com/N283T/gemmi-extra-id/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/N283T/gemmi-extra-id/releases/tag/v0.1.0
