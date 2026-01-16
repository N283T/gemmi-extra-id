from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path

import gemmi


ATOM_SITE_LABEL_TAG = "_atom_site.label_asym_id"
ATOM_SITE_AUTH_TAG = "_atom_site.auth_asym_id"
ATOM_SITE_MOL_TAG = "_atom_site.molecule_id"
MAP_CATEGORY = "_molecule_id_map."

STRUCT_CONN_CONN_TAG = "_struct_conn.conn_type_id"
STRUCT_CONN_P1_TAG = "_struct_conn.ptnr1_label_asym_id"
STRUCT_CONN_P2_TAG = "_struct_conn.ptnr2_label_asym_id"

DEFAULT_COVALENT_TYPES = {"covale", "disulf"}


def build_chain_list(
    loop: gemmi.cif.Loop,
) -> tuple[list[str], dict[str, str], int]:
    tags = list(loop.tags)
    if ATOM_SITE_LABEL_TAG not in tags:
        raise ValueError(f"Missing required tag: {ATOM_SITE_LABEL_TAG}")
    label_idx = tags.index(ATOM_SITE_LABEL_TAG)
    auth_idx = tags.index(ATOM_SITE_AUTH_TAG) if ATOM_SITE_AUTH_TAG in tags else label_idx

    order: list[str] = []
    label_to_auth: dict[str, str] = {}
    for row_idx in range(loop.length()):
        label = loop[row_idx, label_idx]
        if label not in label_to_auth:
            order.append(label)
            label_to_auth[label] = loop[row_idx, auth_idx]
    return order, label_to_auth, label_idx


def build_struct_conn_edges(
    block: gemmi.cif.Block,
    covalent_types: set[str],
) -> list[tuple[str, str]]:
    col = block.find_loop(STRUCT_CONN_CONN_TAG)
    if col is None:
        return []
    loop = col.get_loop()

    tags = list(loop.tags)
    if STRUCT_CONN_P1_TAG not in tags or STRUCT_CONN_P2_TAG not in tags:
        return []
    conn_idx = tags.index(STRUCT_CONN_CONN_TAG)
    p1_idx = tags.index(STRUCT_CONN_P1_TAG)
    p2_idx = tags.index(STRUCT_CONN_P2_TAG)

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


def build_components(
    chain_order: list[str],
    edges: list[tuple[str, str]],
) -> OrderedDict[str, str]:
    adjacency: dict[str, set[str]] = {chain: set() for chain in chain_order}
    for a, b in edges:
        if a not in adjacency:
            adjacency[a] = set()
            chain_order.append(a)
        if b not in adjacency:
            adjacency[b] = set()
            chain_order.append(b)
        adjacency[a].add(b)
        adjacency[b].add(a)

    mapping: OrderedDict[str, str] = OrderedDict()
    next_id = 0
    for chain in chain_order:
        if chain in mapping:
            continue
        stack = [chain]
        while stack:
            current = stack.pop()
            if current in mapping:
                continue
            mapping[current] = str(next_id)
            stack.extend(adjacency.get(current, []))
        next_id += 1
    return mapping


def add_molecule_id_to_atom_site(
    loop: gemmi.cif.Loop,
    mapping: dict[str, str],
    label_idx: int,
) -> None:
    tags = list(loop.tags)
    if ATOM_SITE_MOL_TAG not in tags:
        loop.add_columns([ATOM_SITE_MOL_TAG], value="?")
        tags = list(loop.tags)
    mol_idx = tags.index(ATOM_SITE_MOL_TAG)

    for row_idx in range(loop.length()):
        label = loop[row_idx, label_idx]
        loop[row_idx, mol_idx] = mapping[label]


def add_mapping_category(
    block: gemmi.cif.Block,
    chain_order: list[str],
    label_to_auth: dict[str, str],
    mapping: dict[str, str],
) -> None:
    loop = block.init_mmcif_loop(MAP_CATEGORY, ["label_asym_id", "auth_asym_id", "molecule_id"])
    for label in chain_order:
        loop.add_row([label, label_to_auth.get(label, "?"), mapping[label]])


def run(input_path: Path, output_path: Path, covalent_types: set[str]) -> None:
    doc = gemmi.cif.read(str(input_path))
    block = doc.sole_block()
    atom_site_loop = block.find_loop(ATOM_SITE_LABEL_TAG).get_loop()

    chain_order, label_to_auth, label_idx = build_chain_list(atom_site_loop)
    edges = build_struct_conn_edges(block, covalent_types)
    mapping = build_components(chain_order, edges)

    add_molecule_id_to_atom_site(atom_site_loop, mapping, label_idx)
    add_mapping_category(block, chain_order, label_to_auth, mapping)

    doc.write_file(str(output_path))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Append molecule_id based on struct_conn covalent links and write mmCIF."
    )
    parser.add_argument("input", nargs="?", default="data/148L.cif")
    parser.add_argument("output", nargs="?", default="data/148L_struct_conn_molecule_id.cif")
    parser.add_argument(
        "--conn-types",
        default=",".join(sorted(DEFAULT_COVALENT_TYPES)),
        help="Comma-separated conn_type_id values to treat as covalent.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    covalent_types = {t.strip().lower() for t in args.conn_types.split(",") if t.strip()}
    run(input_path, output_path, covalent_types)


if __name__ == "__main__":
    main()
