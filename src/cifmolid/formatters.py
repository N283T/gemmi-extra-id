"""Output formatters for molecule_id mappings."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import asdict
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections import OrderedDict

    from cifmolid.mmcif import AssignmentResult


def to_json(mapping: OrderedDict[str, int], indent: int = 2) -> str:
    """Convert mapping to JSON string."""
    return json.dumps(dict(mapping), indent=indent)


def to_csv(mapping: OrderedDict[str, int]) -> str:
    """Convert mapping to CSV string."""
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(["label_asym_id", "molecule_id"])
    for chain, mol_id in mapping.items():
        writer.writerow([chain, mol_id])
    return output.getvalue()


def to_tsv(mapping: OrderedDict[str, int]) -> str:
    """Convert mapping to TSV string."""
    output = StringIO()
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(["label_asym_id", "molecule_id"])
    for chain, mol_id in mapping.items():
        writer.writerow([chain, mol_id])
    return output.getvalue()


def to_table(mapping: OrderedDict[str, int]) -> str:
    """Convert mapping to human-readable table string."""
    lines = ["label_asym_id  molecule_id", "-" * 26]
    for chain, mol_id in mapping.items():
        lines.append(f"{chain:<14} {mol_id}")
    return "\n".join(lines)


def write_output(
    mapping: OrderedDict[str, int],
    output_path: str | Path | None,
    fmt: str = "table",
) -> None:
    """
    Write mapping to file or stdout in specified format.

    Args:
        mapping: The molecule_id mapping to write.
        output_path: Path to write to, or "-" for stdout, or None for no output.
        fmt: Output format ("json", "csv", "tsv", "table").
    """
    if output_path is None:
        return

    formatters = {
        "json": to_json,
        "csv": to_csv,
        "tsv": to_tsv,
        "table": to_table,
    }

    if fmt not in formatters:
        raise ValueError(f"Unknown format: {fmt}. Use one of: {', '.join(formatters)}")

    content = formatters[fmt](mapping)

    if str(output_path) == "-":
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")
    else:
        Path(output_path).write_text(content)


# Extended formatters for AssignmentResult
_EXTENDED_HEADERS = [
    "label_asym_id",
    "auth_asym_id",
    "entity_id",
    "entity_type",
    "molecule_id",
    "pn_unit_id",
    "pn_unit_entity",
    "molecule_entity",
]


def to_extended_json(result: AssignmentResult, indent: int = 2) -> str:
    """Convert full assignment result to JSON."""
    data = {chain: asdict(info) for chain, info in result.chain_info.items()}
    return json.dumps(data, indent=indent)


def to_extended_csv(result: AssignmentResult) -> str:
    """Convert full assignment result to CSV with all columns."""
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(_EXTENDED_HEADERS)
    for info in result.chain_info.values():
        writer.writerow(
            [
                info.label_asym_id,
                info.auth_asym_id,
                info.entity_id,
                info.entity_type,
                info.molecule_id,
                info.pn_unit_id,
                info.pn_unit_entity,
                info.molecule_entity,
            ]
        )
    return output.getvalue()


def to_extended_tsv(result: AssignmentResult) -> str:
    """Convert full assignment result to TSV with all columns."""
    output = StringIO()
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(_EXTENDED_HEADERS)
    for info in result.chain_info.values():
        writer.writerow(
            [
                info.label_asym_id,
                info.auth_asym_id,
                info.entity_id,
                info.entity_type,
                info.molecule_id,
                info.pn_unit_id,
                info.pn_unit_entity,
                info.molecule_entity,
            ]
        )
    return output.getvalue()


def to_extended_table(result: AssignmentResult) -> str:
    """Convert full assignment result to human-readable table string."""
    # Calculate column widths
    widths = {h: len(h) for h in _EXTENDED_HEADERS}
    for info in result.chain_info.values():
        widths["label_asym_id"] = max(widths["label_asym_id"], len(info.label_asym_id))
        widths["auth_asym_id"] = max(widths["auth_asym_id"], len(info.auth_asym_id))
        widths["entity_id"] = max(widths["entity_id"], len(info.entity_id))
        widths["entity_type"] = max(widths["entity_type"], len(info.entity_type))
        widths["molecule_id"] = max(widths["molecule_id"], len(str(info.molecule_id)))
        widths["pn_unit_id"] = max(widths["pn_unit_id"], len(info.pn_unit_id))
        widths["pn_unit_entity"] = max(widths["pn_unit_entity"], len(info.pn_unit_entity))
        widths["molecule_entity"] = max(widths["molecule_entity"], len(info.molecule_entity))

    # Build format string
    fmt_parts = [f"{{:<{widths[h]}}}" for h in _EXTENDED_HEADERS]
    fmt_str = "  ".join(fmt_parts)

    lines = [fmt_str.format(*_EXTENDED_HEADERS)]
    total_width = sum(widths.values()) + 2 * (len(_EXTENDED_HEADERS) - 1)
    lines.append("-" * total_width)

    for info in result.chain_info.values():
        lines.append(
            fmt_str.format(
                info.label_asym_id,
                info.auth_asym_id,
                info.entity_id,
                info.entity_type,
                str(info.molecule_id),
                info.pn_unit_id,
                info.pn_unit_entity,
                info.molecule_entity,
            )
        )

    return "\n".join(lines)


def to_extended_tree(result: AssignmentResult) -> str:
    """Convert full assignment result to hierarchical tree string.

    Shows the relationship: molecule → pn_unit → chain
    """
    from collections import defaultdict

    from rich.console import Console
    from rich.tree import Tree

    # Group chains by molecule_id, then by pn_unit_id
    molecules: dict[int, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    molecule_entities: dict[int, str] = {}

    for info in result.chain_info.values():
        mol_id = info.molecule_id
        molecules[mol_id][info.pn_unit_id].append(info)
        molecule_entities[mol_id] = info.molecule_entity

    # Build tree
    root = Tree("[bold]ID Hierarchy[/bold]")

    for mol_id in sorted(molecules.keys()):
        mol_entity = molecule_entities[mol_id]
        mol_branch = root.add(
            f"[bold cyan]molecule {mol_id}[/bold cyan] [dim](entity: {mol_entity})[/dim]"
        )

        pn_units = molecules[mol_id]
        for pn_unit_id in sorted(pn_units.keys()):
            chains = pn_units[pn_unit_id]
            # Get pn_unit info from first chain
            first_chain = chains[0]
            etype = first_chain.entity_type
            pn_ent = first_chain.pn_unit_entity
            pn_label = f"[green]pn_unit {pn_unit_id}[/green] [dim]({etype}, entity: {pn_ent})[/dim]"
            pn_branch = mol_branch.add(pn_label)

            for info in chains:
                pn_branch.add(
                    f"[yellow]{info.label_asym_id}[/yellow] "
                    f"[dim](auth: {info.auth_asym_id}, entity: {info.entity_id})[/dim]"
                )

    # Render to string
    console = Console(force_terminal=False, no_color=False, width=120)
    with console.capture() as capture:
        console.print(root)
    return capture.get()


def print_extended_tree(result: AssignmentResult) -> None:
    """Print hierarchical tree directly to console with colors."""
    from collections import defaultdict

    from rich.console import Console
    from rich.tree import Tree

    # Group chains by molecule_id, then by pn_unit_id
    molecules: dict[int, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    molecule_entities: dict[int, str] = {}

    for info in result.chain_info.values():
        mol_id = info.molecule_id
        molecules[mol_id][info.pn_unit_id].append(info)
        molecule_entities[mol_id] = info.molecule_entity

    # Build tree
    root = Tree("[bold]ID Hierarchy[/bold]")

    for mol_id in sorted(molecules.keys()):
        mol_entity = molecule_entities[mol_id]
        mol_branch = root.add(
            f"[bold cyan]molecule {mol_id}[/bold cyan] [dim](entity: {mol_entity})[/dim]"
        )

        pn_units = molecules[mol_id]
        for pn_unit_id in sorted(pn_units.keys()):
            chains = pn_units[pn_unit_id]
            first_chain = chains[0]
            etype = first_chain.entity_type
            pn_ent = first_chain.pn_unit_entity
            pn_label = f"[green]pn_unit {pn_unit_id}[/green] [dim]({etype}, entity: {pn_ent})[/dim]"
            pn_branch = mol_branch.add(pn_label)

            for info in chains:
                pn_branch.add(
                    f"[yellow]{info.label_asym_id}[/yellow] "
                    f"[dim](auth: {info.auth_asym_id}, entity: {info.entity_id})[/dim]"
                )

    Console().print(root)


def write_extended_output(
    result: AssignmentResult,
    output_path: str | Path | None,
    fmt: str = "table",
) -> None:
    """
    Write extended result to file or stdout in specified format.

    Args:
        result: The AssignmentResult to write.
        output_path: Path to write to, or "-" for stdout, or None for no output.
        fmt: Output format ("json", "csv", "tsv", "table").
    """
    if output_path is None:
        return

    # Tree format is special - print directly for colors
    if fmt == "tree":
        if str(output_path) == "-":
            print_extended_tree(result)
        else:
            content = to_extended_tree(result)
            Path(output_path).write_text(content)
        return

    formatters = {
        "json": to_extended_json,
        "csv": to_extended_csv,
        "tsv": to_extended_tsv,
        "table": to_extended_table,
    }

    if fmt not in formatters:
        raise ValueError(f"Unknown format: {fmt}. Use one of: {', '.join(formatters)}, tree")

    content = formatters[fmt](result)

    if str(output_path) == "-":
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")
    else:
        Path(output_path).write_text(content)
