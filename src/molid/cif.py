"""mmCIF I/O operations using Gemmi."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi

from molid.core import DEFAULT_COVALENT_TYPES, find_components

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
    """Extract covalent bond edges from struct_conn loop."""
    col = block.find_loop(_STRUCT_CONN_TYPE)
    if col is None:
        return []

    loop = col.get_loop()
    tags = list(loop.tags)

    if _STRUCT_CONN_P1 not in tags or _STRUCT_CONN_P2 not in tags:
        return []

    conn_idx = tags.index(_STRUCT_CONN_TYPE)
    p1_idx = tags.index(_STRUCT_CONN_P1)
    p2_idx = tags.index(_STRUCT_CONN_P2)

    edges: list[tuple[str, str]] = []
    for row_idx in range(loop.length()):
        conn_type = loop[row_idx, conn_idx].lower()
        if conn_type not in covalent_types:
            continue
        a = loop[row_idx, p1_idx]
        b = loop[row_idx, p2_idx]
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
