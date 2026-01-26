"""Command-line interface for gemmi_extra_id."""

from __future__ import annotations

import sys

try:
    import typer
    from rich.console import Console
except ImportError:
    print(
        "CLI dependencies (typer, rich) not installed.\n"
        "Install with: pip install gemmi-extra-id[cli]",
        file=sys.stderr,
    )
    sys.exit(1)

from enum import Enum
from pathlib import Path
from typing import Annotated

from gemmi_extra_id import __version__
from gemmi_extra_id.formatters import write_extended_output, write_output
from gemmi_extra_id.graph import DEFAULT_COVALENT_TYPES
from gemmi_extra_id.mmcif import (
    VALID_SWAP_TARGETS,
    assign_extended_ids,
    assign_molecule_id,
    swap_auth_asym_id,
)

app = typer.Typer(
    name="gemmi-extra-id",
    help="Assign extra IDs to mmCIF files based on covalent connectivity.",
    no_args_is_help=True,
)
console = Console()
err_console = Console(stderr=True)


class OutputFormat(str, Enum):
    """Output format options."""

    cif = "cif"
    json = "json"
    csv = "csv"
    tsv = "tsv"
    table = "table"
    tree = "tree"


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        console.print(f"gemmi-extra-id {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool | None,
        typer.Option("--version", "-V", callback=version_callback, is_eager=True),
    ] = None,
) -> None:
    """Assign molecule_id to mmCIF files based on covalent connectivity."""


@app.command()
def assign(
    input_file: Annotated[
        Path,
        typer.Argument(help="Input mmCIF file", exists=True, dir_okay=False),
    ],
    output_file: Annotated[
        str | None,
        typer.Argument(help="Output file (default: <input>_molid.cif, use '-' for stdout)"),
    ] = None,
    fmt: Annotated[
        OutputFormat,
        typer.Option(
            "--format",
            "-f",
            help="Output format",
        ),
    ] = OutputFormat.cif,
    conn_types: Annotated[
        str,
        typer.Option(
            "--conn-types",
            "-c",
            help="Comma-separated conn_type_id values to treat as covalent",
        ),
    ] = ",".join(sorted(DEFAULT_COVALENT_TYPES)),
    quiet: Annotated[
        bool,
        typer.Option(
            "--quiet",
            "-q",
            help="Suppress status messages (useful with stdout output)",
        ),
    ] = False,
    extended: Annotated[
        bool,
        typer.Option(
            "--extended",
            "-e",
            help="Include all IDs: chain_entity, pn_unit_id, pn_unit_entity, molecule_entity",
        ),
    ] = False,
    swap: Annotated[
        str | None,
        typer.Option(
            "--swap",
            "-s",
            help=(
                "Swap auth_asym_id with specified ID. "
                f"Options: {', '.join(sorted(VALID_SWAP_TARGETS))}"
            ),
        ),
    ] = None,
    hash_entities: Annotated[
        bool,
        typer.Option(
            "--hash-entities",
            "-H",
            help="Use Weisfeiler-Lehman graph hashing for entity IDs. "
            "Requires networkx: pip install gemmi-extra-id[hash]",
        ),
    ] = False,
) -> None:
    """Assign molecule_id to an mmCIF file based on covalent connectivity."""
    covalent_types = {t.strip().lower() for t in conn_types.split(",") if t.strip()}

    # Tree format requires extended mode
    if fmt == OutputFormat.tree:
        extended = True

    # --hash-entities implies extended mode
    if hash_entities:
        extended = True

    # Determine output path
    is_stdout = output_file == "-"
    if output_file is None:
        ext = ".cif" if fmt == OutputFormat.cif else f".{fmt.value}"
        output_path: Path | str = input_file.with_stem(f"{input_file.stem}_molid").with_suffix(ext)
    else:
        output_path = "-" if is_stdout else Path(output_file)

    # Validate --swap option
    if swap is not None:
        if swap not in VALID_SWAP_TARGETS:
            err_console.print(
                f"[red]Error:[/red] Invalid --swap value: {swap}. "
                f"Valid options: {', '.join(sorted(VALID_SWAP_TARGETS))}"
            )
            raise typer.Exit(1)

        if fmt != OutputFormat.cif:
            err_console.print("[red]Error:[/red] --swap option only works with CIF output format")
            raise typer.Exit(1)

        if is_stdout:
            err_console.print("[red]Error:[/red] --swap option cannot write to stdout")
            raise typer.Exit(1)

    try:
        if swap is not None:
            # Swap mode: replace auth_asym_id with specified ID
            result = swap_auth_asym_id(
                input_file,
                output_path,
                swap_with=swap,
                covalent_types=covalent_types,
            )
            mapping = result.molecule_id_mapping
        elif extended:
            # Extended mode: include all IDs
            if fmt == OutputFormat.cif:
                result = assign_extended_ids(
                    input_file,
                    None if is_stdout else output_path,
                    covalent_types,
                    use_hash=hash_entities,
                )
                if is_stdout:
                    err_console.print("[red]Error:[/red] CIF format cannot be written to stdout")
                    raise typer.Exit(1)
            else:
                result = assign_extended_ids(
                    input_file, None, covalent_types, use_hash=hash_entities
                )
                write_extended_output(result, output_path, fmt.value)
            mapping = result.molecule_id_mapping
        else:
            # Standard mode: molecule_id only
            if fmt == OutputFormat.cif:
                mapping = assign_molecule_id(
                    input_file, None if is_stdout else output_path, covalent_types
                )
                if is_stdout:
                    err_console.print("[red]Error:[/red] CIF format cannot be written to stdout")
                    raise typer.Exit(1)
            else:
                mapping = assign_molecule_id(input_file, None, covalent_types)
                write_output(mapping, output_path, fmt.value)
    except ImportError as e:
        err_console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1) from None
    except ValueError as e:
        err_console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1) from None

    if not quiet and not is_stdout:
        n_molecules = len(set(mapping.values()))
        n_chains = len(mapping)
        console.print(f"[green]Wrote[/green] {output_path}")
        console.print(
            f"Assigned [cyan]{n_molecules}[/cyan] molecule_id(s) "
            f"to [cyan]{n_chains}[/cyan] chain(s)"
        )
        if swap is not None:
            console.print(f"Swapped auth_asym_id with [cyan]{swap}[/cyan]")
        if hash_entities:
            console.print("Used [cyan]graph hashing[/cyan] for entity IDs")


if __name__ == "__main__":
    app()
