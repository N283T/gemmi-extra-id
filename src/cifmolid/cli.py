"""Command-line interface for cifmolid."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from cifmolid import __version__
from cifmolid.formatters import write_extended_output, write_output
from cifmolid.graph import DEFAULT_COVALENT_TYPES
from cifmolid.mmcif import assign_extended_ids, assign_molecule_id

app = typer.Typer(
    name="cifmolid",
    help="Assign molecule_id to mmCIF files based on covalent connectivity.",
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
        console.print(f"cifmolid {__version__}")
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
) -> None:
    """Assign molecule_id to an mmCIF file based on covalent connectivity."""
    covalent_types = {t.strip().lower() for t in conn_types.split(",") if t.strip()}

    # Tree format requires extended mode
    if fmt == OutputFormat.tree:
        extended = True

    # Determine output path
    is_stdout = output_file == "-"
    if output_file is None:
        ext = ".cif" if fmt == OutputFormat.cif else f".{fmt.value}"
        output_path: Path | str = input_file.with_stem(f"{input_file.stem}_molid").with_suffix(ext)
    else:
        output_path = "-" if is_stdout else Path(output_file)

    try:
        if extended:
            # Extended mode: include all IDs
            if fmt == OutputFormat.cif:
                result = assign_extended_ids(
                    input_file, None if is_stdout else output_path, covalent_types
                )
                if is_stdout:
                    err_console.print("[red]Error:[/red] CIF format cannot be written to stdout")
                    raise typer.Exit(1)
            else:
                result = assign_extended_ids(input_file, None, covalent_types)
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


if __name__ == "__main__":
    app()
