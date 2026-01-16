"""Command-line interface for molid."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from molid import __version__
from molid.cif import assign_molecule_id
from molid.core import DEFAULT_COVALENT_TYPES

app = typer.Typer(
    name="molid",
    help="Assign molecule_id to mmCIF files based on covalent connectivity.",
    no_args_is_help=True,
)
console = Console()
err_console = Console(stderr=True)


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        console.print(f"molid {__version__}")
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
        Path | None,
        typer.Argument(help="Output mmCIF file (default: <input>_molid.cif)"),
    ] = None,
    conn_types: Annotated[
        str,
        typer.Option(
            "--conn-types",
            "-c",
            help="Comma-separated conn_type_id values to treat as covalent",
        ),
    ] = ",".join(sorted(DEFAULT_COVALENT_TYPES)),
) -> None:
    """Assign molecule_id to an mmCIF file based on covalent connectivity."""
    if output_file is None:
        output_file = input_file.with_stem(f"{input_file.stem}_molid")

    covalent_types = {t.strip().lower() for t in conn_types.split(",") if t.strip()}

    try:
        mapping = assign_molecule_id(input_file, output_file, covalent_types)
    except ValueError as e:
        err_console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1) from None

    n_molecules = len(set(mapping.values()))
    n_chains = len(mapping)
    console.print(f"[green]Wrote[/green] {output_file}")
    console.print(
        f"Assigned [cyan]{n_molecules}[/cyan] molecule_id(s) to [cyan]{n_chains}[/cyan] chain(s)"
    )


if __name__ == "__main__":
    app()
