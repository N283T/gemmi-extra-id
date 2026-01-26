"""mmCIF I/O operations using Gemmi."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi

from gemmi_extra_id.graph import DEFAULT_COVALENT_TYPES, find_components

if TYPE_CHECKING:
    from collections.abc import Set as AbstractSet


# mmCIF tags
_ATOM_SITE_LABEL = "_atom_site.label_asym_id"
_ATOM_SITE_AUTH = "_atom_site.auth_asym_id"
_ATOM_SITE_MOLID = "_atom_site.molecule_id"
_STRUCT_CONN_TYPE = "_struct_conn.conn_type_id"
_STRUCT_CONN_P1 = "_struct_conn.ptnr1_label_asym_id"
_STRUCT_CONN_P2 = "_struct_conn.ptnr2_label_asym_id"
_MOLID_MAP_CATEGORY = "_molecule_id_map."

# Swap feature
_ATOM_SITE_ORIG_AUTH = "_atom_site.orig_auth_asym_id"
VALID_SWAP_TARGETS = frozenset({"molecule_id"})


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
    """Extract covalent bond edges from struct_conn."""
    conn_col = block.find_values(_STRUCT_CONN_TYPE)
    if not conn_col:
        return []

    p1_col = block.find_values(_STRUCT_CONN_P1)
    p2_col = block.find_values(_STRUCT_CONN_P2)

    if not p1_col or not p2_col:
        return []

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
            Defaults to {"covale"}.

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


def swap_auth_asym_id(
    input_path: str | Path,
    output_path: str | Path,
    preserve_original: bool = True,
    covalent_types: AbstractSet[str] | None = None,
) -> OrderedDict[str, int]:
    """
    Assign molecule_id and swap auth_asym_id with it.

    This function computes molecule_id and replaces the auth_asym_id column
    with it. The original auth_asym_id values are optionally preserved.

    Args:
        input_path: Path to input mmCIF file.
        output_path: Path to output mmCIF file.
        preserve_original: If True, store original auth_asym_id in
            _atom_site.orig_auth_asym_id.
        covalent_types: Set of conn_type_id values to treat as covalent bonds.

    Returns:
        OrderedDict mapping label_asym_id to molecule_id.

    Raises:
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

    # Perform the swap
    tags = list(atom_site_loop.tags)

    # Get or create auth_asym_id index
    if _ATOM_SITE_AUTH in tags:
        auth_idx = tags.index(_ATOM_SITE_AUTH)
    else:
        atom_site_loop.add_columns([_ATOM_SITE_AUTH], value="?")
        tags = list(atom_site_loop.tags)
        auth_idx = tags.index(_ATOM_SITE_AUTH)

    # Optionally preserve original auth_asym_id
    if preserve_original:
        if _ATOM_SITE_ORIG_AUTH not in tags:
            atom_site_loop.add_columns([_ATOM_SITE_ORIG_AUTH], value="?")
            tags = list(atom_site_loop.tags)
        orig_auth_idx = tags.index(_ATOM_SITE_ORIG_AUTH)

        for row_idx in range(atom_site_loop.length()):
            atom_site_loop[row_idx, orig_auth_idx] = atom_site_loop[row_idx, auth_idx]

    # Swap auth_asym_id with molecule_id
    for row_idx in range(atom_site_loop.length()):
        label = atom_site_loop[row_idx, label_idx]
        atom_site_loop[row_idx, auth_idx] = str(mapping[label])

    # Add molecule_id column and mapping loop
    _add_molecule_id_column(atom_site_loop, mapping, label_idx)
    _add_mapping_loop(block, chain_order, label_to_auth, mapping)

    doc.write_file(str(output_path))

    return mapping
