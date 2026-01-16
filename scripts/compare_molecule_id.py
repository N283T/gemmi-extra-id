from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import gemmi

ATOM_SITE_ID_TAG = "_atom_site.id"
ATOM_SITE_LABEL_TAG = "_atom_site.label_asym_id"
ATOM_SITE_AUTH_TAG = "_atom_site.auth_asym_id"
ATOM_SITE_MOL_TAG = "_atom_site.molecule_id"

DEFAULT_KEY_FIELDS = [
    "_atom_site.label_asym_id",
    "_atom_site.label_seq_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_atom_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_atom_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.label_alt_id",
    "_atom_site.pdbx_PDB_model_num",
]

MISSING_VALUES = {None, "?", "."}


@dataclass(frozen=True)
class IndexResult:
    mapping: dict[tuple[str, ...], set[str | None]]
    missing: int
    rows: int


def find_atom_site_loop(block: gemmi.cif.Block) -> gemmi.cif.Loop:
    for tag in (ATOM_SITE_ID_TAG, "_atom_site.group_PDB", ATOM_SITE_LABEL_TAG):
        col = block.find_loop(tag)
        if col is not None:
            return col.get_loop()
    raise ValueError("Could not find _atom_site loop in CIF.")


def is_unique(loop: gemmi.cif.Loop, idx: int) -> bool:
    seen: set[str] = set()
    for row_idx in range(loop.length()):
        val = loop[row_idx, idx]
        if val in seen:
            return False
        seen.add(val)
    return True


def choose_key_fields(
    loop_a: gemmi.cif.Loop, loop_b: gemmi.cif.Loop
) -> list[str]:
    tags_a = set(loop_a.tags)
    tags_b = set(loop_b.tags)
    if ATOM_SITE_ID_TAG in tags_a and ATOM_SITE_ID_TAG in tags_b:
        idx_a = list(loop_a.tags).index(ATOM_SITE_ID_TAG)
        idx_b = list(loop_b.tags).index(ATOM_SITE_ID_TAG)
        if is_unique(loop_a, idx_a) and is_unique(loop_b, idx_b):
            return [ATOM_SITE_ID_TAG]

    fields = [f for f in DEFAULT_KEY_FIELDS if f in tags_a and f in tags_b]
    if not fields:
        raise ValueError("No shared key fields found for comparison.")
    return fields


def build_index(
    loop: gemmi.cif.Loop, key_fields: Iterable[str]
) -> IndexResult:
    tags = list(loop.tags)
    key_idxs = [tags.index(f) for f in key_fields]
    mol_idx = tags.index(ATOM_SITE_MOL_TAG) if ATOM_SITE_MOL_TAG in tags else None

    mapping: dict[tuple[str, ...], set[str | None]] = defaultdict(set)
    missing = 0
    for row_idx in range(loop.length()):
        key = tuple(loop[row_idx, idx] for idx in key_idxs)
        mol = loop[row_idx, mol_idx] if mol_idx is not None else None
        if mol in MISSING_VALUES:
            missing += 1
        mapping[key].add(mol)
    return IndexResult(mapping=mapping, missing=missing, rows=loop.length())


def build_pair_map(loop: gemmi.cif.Loop) -> dict[tuple[str, str], set[str | None]] | None:
    tags = list(loop.tags)
    if ATOM_SITE_LABEL_TAG not in tags:
        return None
    label_idx = tags.index(ATOM_SITE_LABEL_TAG)
    auth_idx = tags.index(ATOM_SITE_AUTH_TAG) if ATOM_SITE_AUTH_TAG in tags else label_idx
    mol_idx = tags.index(ATOM_SITE_MOL_TAG) if ATOM_SITE_MOL_TAG in tags else None

    mapping: dict[tuple[str, str], set[str | None]] = defaultdict(set)
    for row_idx in range(loop.length()):
        pair = (loop[row_idx, label_idx], loop[row_idx, auth_idx])
        mol = loop[row_idx, mol_idx] if mol_idx is not None else None
        mapping[pair].add(mol)
    return mapping


def describe_duplicates(mapping: dict[tuple[str, ...], set[str | None]]) -> list[tuple[tuple[str, ...], set[str | None]]]:
    return [(k, v) for k, v in mapping.items() if len(v) > 1]


def format_key(key: tuple[str, ...]) -> str:
    return " | ".join(key)


def format_values(values: set[str | None]) -> str:
    def sort_key(val: str | None) -> tuple[int, str]:
        if val is None:
            return (0, "")
        return (1, str(val))

    return str(sorted(values, key=sort_key))


def print_samples(
    title: str, items: list, max_items: int, formatter
) -> None:
    if not items:
        return
    print(f"{title} (up to {max_items})")
    for item in items[:max_items]:
        print(f"- {formatter(item)}")
    print("")


def compare_files(ref_path: Path, tgt_path: Path, key_fields: list[str] | None, max_items: int) -> None:
    ref_doc = gemmi.cif.read(str(ref_path))
    tgt_doc = gemmi.cif.read(str(tgt_path))
    ref_loop = find_atom_site_loop(ref_doc.sole_block())
    tgt_loop = find_atom_site_loop(tgt_doc.sole_block())

    if key_fields is None:
        key_fields = choose_key_fields(ref_loop, tgt_loop)

    ref = build_index(ref_loop, key_fields)
    tgt = build_index(tgt_loop, key_fields)

    ref_keys = set(ref.mapping.keys())
    tgt_keys = set(tgt.mapping.keys())

    only_ref = ref_keys - tgt_keys
    only_tgt = tgt_keys - ref_keys
    common = ref_keys & tgt_keys

    mismatches = []
    for key in common:
        if ref.mapping[key] != tgt.mapping[key]:
            mismatches.append((key, ref.mapping[key], tgt.mapping[key]))

    ref_dupes = describe_duplicates(ref.mapping)
    tgt_dupes = describe_duplicates(tgt.mapping)

    print("Files")
    print(f"- reference: {ref_path}")
    print(f"- target:    {tgt_path}")
    print("")
    print("Atom site")
    print(f"- key fields: {', '.join(key_fields)}")
    print(f"- rows: ref={ref.rows} tgt={tgt.rows}")
    print(f"- molecule_id missing: ref={ref.missing} tgt={tgt.missing}")
    print(f"- keys only in ref: {len(only_ref)}")
    print(f"- keys only in tgt: {len(only_tgt)}")
    print(f"- mismatched molecule_id: {len(mismatches)}")
    print(f"- duplicate keys (ref): {len(ref_dupes)}")
    print(f"- duplicate keys (tgt): {len(tgt_dupes)}")
    print("")

    print_samples(
        "Sample keys only in ref",
        list(only_ref),
        max_items,
        lambda key: format_key(key),
    )
    print_samples(
        "Sample keys only in tgt",
        list(only_tgt),
        max_items,
        lambda key: format_key(key),
    )
    print_samples(
        "Sample mismatches",
        mismatches,
        max_items,
        lambda item: f"{format_key(item[0])} :: ref={format_values(item[1])} tgt={format_values(item[2])}",
    )
    print_samples(
        "Sample duplicate keys in ref",
        ref_dupes,
        max_items,
        lambda item: f"{format_key(item[0])} :: {format_values(item[1])}",
    )
    print_samples(
        "Sample duplicate keys in tgt",
        tgt_dupes,
        max_items,
        lambda item: f"{format_key(item[0])} :: {format_values(item[1])}",
    )

    ref_pair = build_pair_map(ref_loop)
    tgt_pair = build_pair_map(tgt_loop)
    if ref_pair is not None and tgt_pair is not None:
        pair_keys = set(ref_pair) | set(tgt_pair)
        pair_mismatches = []
        for pair in pair_keys:
            if ref_pair.get(pair, set()) != tgt_pair.get(pair, set()):
                pair_mismatches.append((pair, ref_pair.get(pair, set()), tgt_pair.get(pair, set())))

        print("label/auth asym_id -> molecule_id mapping")
        print(f"- pair mismatches: {len(pair_mismatches)}")
        print_samples(
            "Sample pair mismatches",
            pair_mismatches,
            max_items,
            lambda item: f"{item[0][0]} / {item[0][1]} :: ref={format_values(item[1])} tgt={format_values(item[2])}",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare molecule_id assignments between two mmCIF files."
    )
    parser.add_argument("reference", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument(
        "--key",
        help="Comma-separated list of _atom_site tags to use as the key.",
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=20,
        help="Maximum number of sample rows to show.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    key_fields = [k.strip() for k in args.key.split(",")] if args.key else None
    compare_files(args.reference, args.target, key_fields, args.max_items)


if __name__ == "__main__":
    main()
