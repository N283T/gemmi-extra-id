"""mmCIF I/O operations using Gemmi."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi

from gemmi_extra_id.graph import DEFAULT_COVALENT_TYPES, find_components, find_pn_units

if TYPE_CHECKING:
    from collections.abc import Set as AbstractSet


@dataclass
class ChainInfo:
    """Information about a single chain with all assigned IDs."""

    label_asym_id: str
    auth_asym_id: str
    entity_id: str
    entity_type: str
    molecule_id: int
    pn_unit_id: str
    pn_unit_entity: str
    molecule_entity: str


@dataclass
class AssignmentResult:
    """Result container for extended ID assignment."""

    chain_info: OrderedDict[str, ChainInfo]

    @property
    def molecule_id_mapping(self) -> OrderedDict[str, int]:
        """Backward-compatible molecule_id mapping."""
        return OrderedDict((k, v.molecule_id) for k, v in self.chain_info.items())


# mmCIF tags
_ATOM_SITE_LABEL = "_atom_site.label_asym_id"
_ATOM_SITE_AUTH = "_atom_site.auth_asym_id"
_ATOM_SITE_MOLID = "_atom_site.molecule_id"
_STRUCT_CONN_TYPE = "_struct_conn.conn_type_id"
_STRUCT_CONN_P1 = "_struct_conn.ptnr1_label_asym_id"
_STRUCT_CONN_P2 = "_struct_conn.ptnr2_label_asym_id"
_MOLID_MAP_CATEGORY = "_molecule_id_map."

# Entity tags
_STRUCT_ASYM_ID = "_struct_asym.id"
_STRUCT_ASYM_ENTITY = "_struct_asym.entity_id"
_ENTITY_ID = "_entity.id"
_ENTITY_TYPE = "_entity.type"

# Extended ID tags
_EXTENDED_MAP_CATEGORY = "_extended_id_map."

# Swap feature
_ATOM_SITE_ORIG_AUTH = "_atom_site.orig_auth_asym_id"
VALID_SWAP_TARGETS = frozenset({"molecule_id", "pn_unit_id", "entity_id", "label_asym_id"})


def _get_chain_info(
    loop: gemmi.cif.Loop,
) -> tuple[list[str], dict[str, str], int]:
    """Extract chain order and label-to-auth mapping from atom_site loop."""
    tags = list(loop.tags)
    if _ATOM_SITE_LABEL not in tags:
        raise ValueError(f"Missing required tag: {_ATOM_SITE_LABEL}")

    label_idx = tags.index(_ATOM_SITE_LABEL)
    auth_idx = tags.index(_ATOM_SITE_AUTH) if _ATOM_SITE_AUTH in tags else label_idx

    order: list[str] = []
    label_to_auth: dict[str, str] = {}

    for row_idx in range(loop.length()):
        label = loop[row_idx, label_idx]
        if label not in label_to_auth:
            order.append(label)
            label_to_auth[label] = loop[row_idx, auth_idx]

    return order, label_to_auth, label_idx


def _get_covalent_edges(
    block: gemmi.cif.Block,
    covalent_types: AbstractSet[str],
) -> list[tuple[str, str]]:
    """Extract covalent bond edges from struct_conn.

    Handles both loop format (multiple rows) and key-value format (single row).
    """
    # Use find_values which works for both loops and key-value pairs
    conn_col = block.find_values(_STRUCT_CONN_TYPE)
    if not conn_col:
        return []

    p1_col = block.find_values(_STRUCT_CONN_P1)
    p2_col = block.find_values(_STRUCT_CONN_P2)

    if not p1_col or not p2_col:
        return []

    # All columns should have the same length
    n_rows = len(conn_col)
    if len(p1_col) != n_rows or len(p2_col) != n_rows:
        return []

    edges: list[tuple[str, str]] = []
    for i in range(n_rows):
        conn_type = conn_col[i].lower()
        if conn_type not in covalent_types:
            continue
        a = p1_col[i]
        b = p2_col[i]
        if a in ("?", ".") or b in ("?", "."):
            continue
        edges.append((a, b))

    return edges


def _get_entity_mapping(
    block: gemmi.cif.Block,
) -> dict[str, tuple[str, str]]:
    """
    Extract entity_id and entity_type for each chain.

    Reads _struct_asym and _entity to build:
    label_asym_id -> (entity_id, entity_type)

    Args:
        block: mmCIF data block.

    Returns:
        Dict mapping label_asym_id to (entity_id, entity_type) tuple.
        Returns empty dict if required loops are missing.
    """
    # First, build entity_id -> entity_type mapping from _entity loop
    entity_types: dict[str, str] = {}
    entity_col = block.find_loop(_ENTITY_ID)
    if entity_col is not None:
        loop = entity_col.get_loop()
        if loop is not None:
            tags = list(loop.tags)
            if _ENTITY_TYPE in tags:
                id_idx = tags.index(_ENTITY_ID)
                type_idx = tags.index(_ENTITY_TYPE)
                for row_idx in range(loop.length()):
                    eid = loop[row_idx, id_idx]
                    etype = loop[row_idx, type_idx]
                    if eid not in ("?", "."):
                        entity_types[eid] = etype if etype not in ("?", ".") else "unknown"

    # Then, build label_asym_id -> (entity_id, entity_type) from _struct_asym
    result: dict[str, tuple[str, str]] = {}
    asym_col = block.find_loop(_STRUCT_ASYM_ID)
    if asym_col is not None:
        loop = asym_col.get_loop()
        if loop is not None:
            tags = list(loop.tags)
            if _STRUCT_ASYM_ENTITY in tags:
                id_idx = tags.index(_STRUCT_ASYM_ID)
                entity_idx = tags.index(_STRUCT_ASYM_ENTITY)
                for row_idx in range(loop.length()):
                    asym_id = loop[row_idx, id_idx]
                    entity_id = loop[row_idx, entity_idx]
                    if asym_id not in ("?", ".") and entity_id not in ("?", "."):
                        entity_type = entity_types.get(entity_id, "unknown")
                        result[asym_id] = (entity_id, entity_type)

    return result


def _compute_extended_chain_info(
    chain_order: list[str],
    label_to_auth: dict[str, str],
    edges: list[tuple[str, str]],
    entity_mapping: dict[str, tuple[str, str]],
) -> tuple[OrderedDict[str, ChainInfo], dict[str, int]]:
    """
    Compute all extended IDs for chains.

    This is the shared computation logic used by both assign_extended_ids
    and swap_auth_asym_id.

    Args:
        chain_order: List of label_asym_id in order of appearance.
        label_to_auth: Mapping from label_asym_id to auth_asym_id.
        edges: List of covalent bond edges (chain pairs).
        entity_mapping: Mapping from label_asym_id to (entity_id, entity_type).

    Returns:
        Tuple of (chain_info_dict, molecule_mapping).
    """
    # Compute molecule_id (connected components)
    molecule_mapping = find_components(chain_order, edges)

    # Extract entity types for pn_unit calculation
    entity_types = {chain: entity_mapping.get(chain, ("?", "unknown"))[1] for chain in chain_order}

    # Determine which chains are polymers (entity_type == "polymer")
    is_polymer = {chain: entity_types.get(chain, "unknown") == "polymer" for chain in chain_order}

    # Compute pn_unit_id (AtomWorks-compatible: polymers never grouped)
    pn_unit_mapping = find_pn_units(chain_order, edges, entity_types, is_polymer)

    # Group chains by molecule_id to compute molecule_entity
    molecule_to_chains: dict[int, list[str]] = {}
    for chain, mol_id in molecule_mapping.items():
        if mol_id not in molecule_to_chains:
            molecule_to_chains[mol_id] = []
        molecule_to_chains[mol_id].append(chain)

    # Compute molecule_entity (smallest entity_id in molecule)
    molecule_entity_map: dict[int, str] = {}
    for mol_id, chains in molecule_to_chains.items():
        entity_ids = [
            entity_mapping.get(c, ("?", "unknown"))[0]
            for c in chains
            if entity_mapping.get(c, ("?", "unknown"))[0] != "?"
        ]
        if entity_ids:
            # Use smallest entity_id (numerically if possible, lexically otherwise)
            try:
                molecule_entity_map[mol_id] = str(min(int(e) for e in entity_ids))
            except ValueError:
                molecule_entity_map[mol_id] = min(entity_ids)
        else:
            molecule_entity_map[mol_id] = "?"

    # Build ChainInfo for each chain
    chain_info_dict: OrderedDict[str, ChainInfo] = OrderedDict()
    for chain in chain_order:
        entity_id, entity_type = entity_mapping.get(chain, ("?", "unknown"))
        pn_unit_id = pn_unit_mapping.get(chain, chain)

        # pn_unit_entity: entity_id of the pn_unit (should be same for all chains)
        pn_members = pn_unit_id.split(",")
        pn_entity_ids = [
            entity_mapping.get(c, ("?", "unknown"))[0]
            for c in pn_members
            if entity_mapping.get(c, ("?", "unknown"))[0] != "?"
        ]
        if pn_entity_ids:
            try:
                pn_unit_entity = str(min(int(e) for e in pn_entity_ids))
            except ValueError:
                pn_unit_entity = min(pn_entity_ids)
        else:
            pn_unit_entity = "?"

        mol_id = molecule_mapping[chain]
        chain_info_dict[chain] = ChainInfo(
            label_asym_id=chain,
            auth_asym_id=label_to_auth.get(chain, "?"),
            entity_id=entity_id,
            entity_type=entity_type,
            molecule_id=mol_id,
            pn_unit_id=pn_unit_id,
            pn_unit_entity=pn_unit_entity,
            molecule_entity=molecule_entity_map[mol_id],
        )

    return chain_info_dict, molecule_mapping


def _add_molecule_id_column(
    loop: gemmi.cif.Loop,
    mapping: dict[str, int],
    label_idx: int,
) -> None:
    """Add or update molecule_id column in atom_site loop."""
    tags = list(loop.tags)
    if _ATOM_SITE_MOLID not in tags:
        loop.add_columns([_ATOM_SITE_MOLID], value="?")
        tags = list(loop.tags)

    mol_idx = tags.index(_ATOM_SITE_MOLID)

    for row_idx in range(loop.length()):
        label = loop[row_idx, label_idx]
        loop[row_idx, mol_idx] = str(mapping[label])


def _add_mapping_loop(
    block: gemmi.cif.Block,
    chain_order: list[str],
    label_to_auth: dict[str, str],
    mapping: dict[str, int],
) -> None:
    """Add molecule_id_map loop to block."""
    loop = block.init_mmcif_loop(
        _MOLID_MAP_CATEGORY, ["label_asym_id", "auth_asym_id", "molecule_id"]
    )
    for label in chain_order:
        loop.add_row([label, label_to_auth.get(label, "?"), str(mapping[label])])


def assign_molecule_id(
    input_path: str | Path,
    output_path: str | Path | None = None,
    covalent_types: AbstractSet[str] | None = None,
) -> OrderedDict[str, int]:
    """
    Assign molecule_id to mmCIF file based on covalent connectivity.

    Reads an mmCIF file, identifies connected components of chains based on
    covalent bonds in _struct_conn, and writes the result with molecule_id
    assignments.

    Args:
        input_path: Path to input mmCIF file.
        output_path: Path to output mmCIF file. If None, no file is written.
        covalent_types: Set of conn_type_id values to treat as covalent bonds.
            Defaults to {"covale", "disulf"}.

    Returns:
        OrderedDict mapping label_asym_id to molecule_id (0-indexed).

    Raises:
        ValueError: If required mmCIF tags are missing.
        FileNotFoundError: If input file does not exist.
    """
    if covalent_types is None:
        covalent_types = DEFAULT_COVALENT_TYPES
    else:
        covalent_types = frozenset(t.lower() for t in covalent_types)

    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")

    doc = gemmi.cif.read(str(input_path))
    block = doc.sole_block()

    atom_site_col = block.find_loop(_ATOM_SITE_LABEL)
    if atom_site_col is None:
        raise ValueError(f"Missing {_ATOM_SITE_LABEL} in {input_path}")

    atom_site_loop = atom_site_col.get_loop()
    chain_order, label_to_auth, label_idx = _get_chain_info(atom_site_loop)
    edges = _get_covalent_edges(block, covalent_types)
    mapping = find_components(chain_order, edges)

    if output_path is not None:
        _add_molecule_id_column(atom_site_loop, mapping, label_idx)
        _add_mapping_loop(block, chain_order, label_to_auth, mapping)
        doc.write_file(str(output_path))

    return mapping


def _add_extended_mapping_loop(
    block: gemmi.cif.Block,
    chain_order: list[str],
    chain_info_dict: dict[str, ChainInfo],
) -> None:
    """Add extended_id_map loop to block."""
    loop = block.init_mmcif_loop(
        _EXTENDED_MAP_CATEGORY,
        [
            "label_asym_id",
            "auth_asym_id",
            "entity_id",
            "entity_type",
            "molecule_id",
            "pn_unit_id",
            "pn_unit_entity",
            "molecule_entity",
        ],
    )
    for label in chain_order:
        info = chain_info_dict[label]
        loop.add_row(
            [
                info.label_asym_id,
                info.auth_asym_id,
                info.entity_id,
                info.entity_type,
                str(info.molecule_id),
                info.pn_unit_id,
                info.pn_unit_entity,
                info.molecule_entity,
            ]
        )


def assign_extended_ids(
    input_path: str | Path,
    output_path: str | Path | None = None,
    covalent_types: AbstractSet[str] | None = None,
) -> AssignmentResult:
    """
    Assign all extended IDs to mmCIF file based on covalent connectivity.

    Computes and assigns:
    - molecule_id: Connected component ID (existing behavior)
    - chain_entity: Entity ID for each chain
    - pn_unit_id: Comma-separated chain IDs of same-type chains covalently linked
    - pn_unit_entity: Entity ID for each pn_unit
    - molecule_entity: Entity ID for each molecule (smallest if mixed)

    Args:
        input_path: Path to input mmCIF file.
        output_path: Path to output mmCIF file. If None, no file is written.
        covalent_types: Set of conn_type_id values to treat as covalent bonds.
            Defaults to {"covale", "disulf"}.

    Returns:
        AssignmentResult containing ChainInfo for each chain.

    Raises:
        ValueError: If required mmCIF tags are missing.
        FileNotFoundError: If input file does not exist.
    """
    if covalent_types is None:
        covalent_types = DEFAULT_COVALENT_TYPES
    else:
        covalent_types = frozenset(t.lower() for t in covalent_types)

    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")

    doc = gemmi.cif.read(str(input_path))
    block = doc.sole_block()

    atom_site_col = block.find_loop(_ATOM_SITE_LABEL)
    if atom_site_col is None:
        raise ValueError(f"Missing {_ATOM_SITE_LABEL} in {input_path}")

    atom_site_loop = atom_site_col.get_loop()
    chain_order, label_to_auth, label_idx = _get_chain_info(atom_site_loop)
    edges = _get_covalent_edges(block, covalent_types)
    entity_mapping = _get_entity_mapping(block)

    # Compute all extended IDs using shared helper
    chain_info_dict, molecule_mapping = _compute_extended_chain_info(
        chain_order, label_to_auth, edges, entity_mapping
    )

    if output_path is not None:
        _add_molecule_id_column(atom_site_loop, molecule_mapping, label_idx)
        _add_extended_mapping_loop(block, chain_order, chain_info_dict)
        doc.write_file(str(output_path))

    return AssignmentResult(chain_info=chain_info_dict)


def _swap_auth_asym_id(
    loop: gemmi.cif.Loop,
    chain_info_dict: dict[str, ChainInfo],
    swap_with: str,
    label_idx: int,
    preserve_original: bool = True,
) -> None:
    """
    Swap auth_asym_id with specified ID in atom_site loop.

    Args:
        loop: The _atom_site loop to modify.
        chain_info_dict: ChainInfo for each chain.
        swap_with: Which ID to swap with.
        label_idx: Index of label_asym_id column in loop.
        preserve_original: Whether to preserve original auth_asym_id.
    """
    tags = list(loop.tags)

    # Get or create auth_asym_id index
    if _ATOM_SITE_AUTH in tags:
        auth_idx = tags.index(_ATOM_SITE_AUTH)
    else:
        loop.add_columns([_ATOM_SITE_AUTH], value="?")
        tags = list(loop.tags)
        auth_idx = tags.index(_ATOM_SITE_AUTH)

    # Optionally preserve original auth_asym_id
    if preserve_original:
        if _ATOM_SITE_ORIG_AUTH not in tags:
            loop.add_columns([_ATOM_SITE_ORIG_AUTH], value="?")
            tags = list(loop.tags)
        orig_auth_idx = tags.index(_ATOM_SITE_ORIG_AUTH)

        # Copy current auth values to orig column
        for row_idx in range(loop.length()):
            loop[row_idx, orig_auth_idx] = loop[row_idx, auth_idx]

    # Perform swap
    for row_idx in range(loop.length()):
        label = loop[row_idx, label_idx]
        info = chain_info_dict[label]

        if swap_with == "molecule_id":
            new_value = str(info.molecule_id)
        elif swap_with == "pn_unit_id":
            new_value = info.pn_unit_id
        elif swap_with == "entity_id":
            new_value = info.entity_id
        elif swap_with == "label_asym_id":
            new_value = info.label_asym_id
        else:
            raise ValueError(f"Unknown swap_with value: {swap_with}")

        loop[row_idx, auth_idx] = new_value


def swap_auth_asym_id(
    input_path: str | Path,
    output_path: str | Path,
    swap_with: str = "molecule_id",
    preserve_original: bool = True,
    covalent_types: AbstractSet[str] | None = None,
) -> AssignmentResult:
    """
    Assign IDs and swap auth_asym_id with the specified extra-id.

    This function computes all extended IDs and then replaces the
    auth_asym_id column with the specified ID. The original auth_asym_id
    values are optionally preserved in a new column.

    Args:
        input_path: Path to input mmCIF file.
        output_path: Path to output mmCIF file.
        swap_with: The ID to swap into auth_asym_id:
            - "molecule_id": Connected component ID (integer as string)
            - "pn_unit_id": Same-type covalent unit ID
            - "entity_id": Entity ID from _struct_asym
            - "label_asym_id": Label asym ID
        preserve_original: If True, store original auth_asym_id in
            _atom_site.orig_auth_asym_id.
        covalent_types: Set of conn_type_id values to treat as covalent bonds.

    Returns:
        AssignmentResult containing ChainInfo for each chain.

    Raises:
        ValueError: If swap_with is not a valid target.
        FileNotFoundError: If input file does not exist.
    """
    if swap_with not in VALID_SWAP_TARGETS:
        raise ValueError(
            f"Invalid swap_with value: {swap_with}. "
            f"Valid options: {', '.join(sorted(VALID_SWAP_TARGETS))}"
        )

    if covalent_types is None:
        covalent_types = DEFAULT_COVALENT_TYPES
    else:
        covalent_types = frozenset(t.lower() for t in covalent_types)

    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")

    doc = gemmi.cif.read(str(input_path))
    block = doc.sole_block()

    atom_site_col = block.find_loop(_ATOM_SITE_LABEL)
    if atom_site_col is None:
        raise ValueError(f"Missing {_ATOM_SITE_LABEL} in {input_path}")

    atom_site_loop = atom_site_col.get_loop()
    chain_order, label_to_auth, label_idx = _get_chain_info(atom_site_loop)
    edges = _get_covalent_edges(block, covalent_types)
    entity_mapping = _get_entity_mapping(block)

    # Compute all extended IDs using shared helper
    chain_info_dict, molecule_mapping = _compute_extended_chain_info(
        chain_order, label_to_auth, edges, entity_mapping
    )

    # Perform the swap
    _swap_auth_asym_id(atom_site_loop, chain_info_dict, swap_with, label_idx, preserve_original)

    # Also add molecule_id column and mapping loop
    _add_molecule_id_column(atom_site_loop, molecule_mapping, label_idx)
    _add_extended_mapping_loop(block, chain_order, chain_info_dict)

    doc.write_file(str(output_path))

    return AssignmentResult(chain_info=chain_info_dict)
