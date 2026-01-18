"""Compare gemmi-extra-id with AtomWorks extra-id assignment."""

import contextlib
import io
import json
import logging
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

import gemmi
import typer
from rich.console import Console
from rich.table import Table
from tqdm import tqdm

logger = logging.getLogger(__name__)
console = Console()

# PDB entries with known data quality issues (orphan chains not in _struct_asym)
# These are excluded from comparison as the difference is due to PDB data issues,
# not gemmi-extra-id or AtomWorks behavior.
# See: plans/orphan_chain_detection.md
IGNORE_LIST: set[str] = {
    "2g10",  # Chain F: 147 water atoms not in _struct_asym
    "1ts6",  # Chain C: 7 water atoms not in _struct_asym
    "2k9y",  # Chains C, D: not in _struct_asym
}


def setup_logging(log_file: Path | None = None) -> None:
    """Configure logging to file only (to avoid interfering with tqdm)."""
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    # Suppress all console logging (AtomWorks, Biotite, etc.)
    # Only file handler will receive logs
    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            root_logger.removeHandler(handler)

    # Suppress atomworks/biotite loggers from console
    for lib in ["atomworks", "biotite"]:
        lib_logger = logging.getLogger(lib)
        lib_logger.setLevel(logging.ERROR)
        lib_logger.propagate = False

    logging.captureWarnings(True)


@dataclass
class CompareResult:
    """Result of comparing a single CIF file."""

    pdb_id: str
    status: str  # "pass", "fail", "error"
    gemmi_groups: set[frozenset[str]] = field(default_factory=set)
    atomworks_groups: set[frozenset[str]] = field(default_factory=set)
    error: str | None = None


def find_cif_files(mirror_path: Path) -> list[Path]:
    """Find all CIF files in directory using gemmi.CifWalk."""
    return [Path(cif_str) for cif_str in gemmi.CifWalk(str(mirror_path))]


@contextlib.contextmanager
def suppress_output():
    """Suppress stdout and stderr (for AtomWorks print statements)."""
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        yield


def initialize_chain_info_lite(cif_block, atom_array) -> dict:
    """Lightweight chain_info initialization using only _entity.type.

    This is an alternative to initialize_chain_info_from_category() that doesn't
    require _entity_poly category. For polymer-less entries (branched/non-polymer),
    all chains will have is_polymer=False.
    """
    from atomworks.enums import ChainType

    chain_info_dict = {}

    # Build chain_id -> entity_id mapping from atom_array
    chain_ids = atom_array.get_annotation("chain_id")
    entity_ids = atom_array.get_annotation("label_entity_id").astype(str)
    chain_to_entity = dict(zip(chain_ids, entity_ids, strict=True))

    # Get entity_id -> entity_type mapping from _entity
    entity_types = {}
    if "entity" in cif_block:
        entity_cat = cif_block["entity"]
        ids = entity_cat["id"].as_array()
        types = entity_cat["type"].as_array()
        entity_types = dict(zip(ids.astype(str), types, strict=True))

    # Build chain_info_dict
    for chain_id, entity_id in chain_to_entity.items():
        entity_type = entity_types.get(entity_id, "non-polymer")
        chain_type = ChainType.as_enum(entity_type)
        chain_info_dict[chain_id] = {
            "rcsb_entity": entity_id,
            "chain_type": chain_type,
            "is_polymer": chain_type.is_polymer(),
        }

    return chain_info_dict


def load_with_fixed_nonpoly_res_id(cif_path: Path):
    """Load atom array with non-polymer res_id fixed before adding bonds.

    Flow:
    1. load_any(include_bonds=False) - load without bonds
    2. update_nonpoly_seq_ids() - fix res_id for non-polymers
    3. _add_bonds() - add bonds with fixed res_id
    """
    from atomworks.io.transforms.atom_array import add_polymer_annotation
    from atomworks.io.transforms.categories import initialize_chain_info_from_category
    from atomworks.io.utils.io_utils import _add_bonds, get_structure, read_any
    from biotite.structure import AtomArrayStack

    # Read the CIF file
    cif_file = read_any(cif_path)

    # Load structure WITHOUT bonds
    try:
        atom_array_or_stack = get_structure(
            cif_file,
            extra_fields=["label_entity_id", "auth_seq_id", "occupancy"],
            include_bonds=False,  # Don't add bonds yet!
            model=None,
            altloc="all",
        )
    except Exception:
        # Fall back to first model if models have unequal atoms
        atom_array_or_stack = get_structure(
            cif_file,
            extra_fields=["label_entity_id", "auth_seq_id", "occupancy"],
            include_bonds=False,
            model=1,
            altloc="all",
        )

    # Handle AtomArrayStack (multiple models)
    if isinstance(atom_array_or_stack, AtomArrayStack):
        atom_array = atom_array_or_stack[0]
    else:
        atom_array = atom_array_or_stack

    # Initialize chain_info to get is_polymer annotation
    # Use lightweight version for polymer-less entries (no entity_poly category)
    if "entity_poly" in cif_file.block:
        chain_info = initialize_chain_info_from_category(cif_file.block, atom_array)
    else:
        chain_info = initialize_chain_info_lite(cif_file.block, atom_array)
    atom_array = add_polymer_annotation(atom_array, chain_info)

    # Fix: Update non-polymer res_id to auth_seq_id BEFORE adding bonds
    # This is what parse() does via update_nonpoly_seq_ids()
    if hasattr(atom_array, "auth_seq_id") and hasattr(atom_array, "is_polymer"):
        import numpy as np

        non_polymer_mask = ~atom_array.is_polymer
        auth_seq_ids = atom_array.auth_seq_id
        # Handle string auth_seq_id
        if auth_seq_ids.dtype.kind in ("U", "S", "O"):
            auth_seq_ids = np.array(
                [int(x) if str(x) not in (".", "?", "") else -1 for x in auth_seq_ids]
            )
        atom_array.res_id[non_polymer_mask] = auth_seq_ids[non_polymer_mask]

    # Now add bonds with the fixed res_id
    atom_array = _add_bonds(
        atom_array,
        cif_file.block,
        add_bond_types_from_struct_conn=["covale"],
        fix_bond_types=False,
    )

    return atom_array


def extract_atomworks_molecule_groups(cif_path: Path) -> set[frozenset[str]]:
    """Extract molecule grouping from AtomWorks with fixed non-polymer handling."""
    from atomworks.io.transforms.atom_array import (
        add_molecule_id_annotation,
        add_pn_unit_id_annotation,
    )

    with warnings.catch_warnings(), suppress_output():
        warnings.simplefilter("ignore")
        atom_array = load_with_fixed_nonpoly_res_id(cif_path)
        atom_array = add_pn_unit_id_annotation(atom_array)
        atom_array = add_molecule_id_annotation(atom_array)

    # chain_id -> molecule_id mapping
    mapping: dict[str, int] = {}
    for i in range(len(atom_array)):
        chain = str(atom_array.chain_id[i])
        mol_id = int(atom_array.molecule_id[i])
        if chain not in mapping:
            mapping[chain] = mol_id

    # molecule_id -> set of chains
    groups: dict[int, set[str]] = {}
    for chain, mol_id in mapping.items():
        if mol_id not in groups:
            groups[mol_id] = set()
        groups[mol_id].add(chain)

    return {frozenset(g) for g in groups.values()}


def extract_gemmi_molecule_groups(cif_path: Path) -> set[frozenset[str]]:
    """Extract molecule grouping from gemmi-extra-id."""
    from gemmi_extra_id import assign_extended_ids

    result = assign_extended_ids(cif_path)

    # molecule_id -> set of chains
    groups: dict[int, set[str]] = {}
    for chain, info in result.chain_info.items():
        mol_id = info.molecule_id
        if mol_id not in groups:
            groups[mol_id] = set()
        groups[mol_id].add(chain)

    return {frozenset(g) for g in groups.values()}


def compare_one(cif_path: Path) -> CompareResult:
    """Compare molecule groupings for a single CIF file."""
    pdb_id = cif_path.stem.replace(".cif", "")
    try:
        gemmi_groups = extract_gemmi_molecule_groups(cif_path)
        atomworks_groups = extract_atomworks_molecule_groups(cif_path)

        if gemmi_groups == atomworks_groups:
            return CompareResult(pdb_id, "pass")
        else:
            return CompareResult(
                pdb_id,
                "fail",
                gemmi_groups=gemmi_groups,
                atomworks_groups=atomworks_groups,
            )
    except Exception as e:
        return CompareResult(pdb_id, "error", error=str(e))


def format_groups(groups: set[frozenset[str]]) -> str:
    """Format molecule groups for display."""
    sorted_groups = sorted([sorted(g) for g in groups], key=lambda x: x[0] if x else "")
    return ", ".join("{" + ",".join(g) + "}" for g in sorted_groups)


app = typer.Typer(help="Compare gemmi-extra-id with AtomWorks extra-id assignment.")


@app.command()
def main(
    pdb_mirror: Annotated[
        Path, typer.Argument(help="Path to PDB mmCIF mirror directory or test data")
    ],
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Output JSON file for results")
    ] = None,
    workers: Annotated[
        int, typer.Option("--workers", "-w", help="Number of parallel workers")
    ] = os.cpu_count() or 4,
    log: Annotated[Path | None, typer.Option("--log", "-l", help="Log file path")] = None,
    limit: Annotated[
        int | None, typer.Option("--limit", "-n", help="Limit number of files to process")
    ] = None,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Show details of failed cases")
    ] = False,
    include_ignored: Annotated[
        bool,
        typer.Option("--include-ignored", help="Include entries in IGNORE_LIST"),
    ] = False,
) -> None:
    """Compare gemmi-extra-id with AtomWorks extra-id assignment.

    Uses AtomWorks' lightweight load_any() function instead of parse()
    to avoid heavy preprocessing while still getting molecule_id assignment.
    """
    setup_logging(log)

    # Find all CIF files
    console.print(f"Scanning [cyan]{pdb_mirror}[/cyan] for CIF files...")
    cif_files = find_cif_files(pdb_mirror)

    # Filter out entries in IGNORE_LIST (known PDB data quality issues)
    ignored_count = 0
    if not include_ignored:
        original_count = len(cif_files)
        cif_files = [f for f in cif_files if f.stem.split(".")[0] not in IGNORE_LIST]
        ignored_count = original_count - len(cif_files)

    if limit:
        cif_files = cif_files[:limit]

    total = len(cif_files)
    ignore_msg = f" (ignored {ignored_count})" if ignored_count > 0 else ""
    console.print(
        f"Found [green]{total}[/green] CIF files{ignore_msg}, "
        f"processing with [yellow]{workers}[/yellow] workers..."
    )

    logger.info(f"Found {total} CIF files{ignore_msg}, processing with {workers} workers...")

    passed: list[str] = []
    failed: list[CompareResult] = []
    errors: list[CompareResult] = []

    # Process files in parallel with tqdm progress bar
    with (
        ProcessPoolExecutor(max_workers=workers) as executor,
        tqdm(total=total, desc="Processing", unit="file") as pbar,
    ):
        for result in executor.map(compare_one, cif_files):
            if result.status == "pass":
                passed.append(result.pdb_id)
                logger.info(f"PASS: {result.pdb_id}")
            elif result.status == "fail":
                failed.append(result)
                logger.warning(
                    f"FAIL: {result.pdb_id} - "
                    f"gemmi={format_groups(result.gemmi_groups)} vs "
                    f"atomworks={format_groups(result.atomworks_groups)}"
                )
                if verbose:
                    tqdm.write(
                        f"FAIL: {result.pdb_id}: "
                        f"gemmi={format_groups(result.gemmi_groups)} vs "
                        f"atomworks={format_groups(result.atomworks_groups)}"
                    )
            else:
                errors.append(result)
                logger.error(f"ERROR: {result.pdb_id} - {result.error}")
                if verbose:
                    tqdm.write(f"ERROR: {result.pdb_id} - {result.error}")
            pbar.update(1)

    # Summary table
    console.print()
    table = Table(title="Summary")
    table.add_column("Status", style="bold")
    table.add_column("Count", justify="right")
    table.add_row("[green]Pass[/green]", str(len(passed)))
    table.add_row("[yellow]Fail[/yellow]", str(len(failed)))
    table.add_row("[red]Error[/red]", str(len(errors)))
    console.print(table)

    # Show failed cases (first 10) if not in verbose mode
    if failed and not verbose:
        console.print()
        console.print("[yellow]Failed cases (showing first 10):[/yellow]")
        for result in failed[:10]:
            console.print(
                f"  {result.pdb_id}: "
                f"gemmi={format_groups(result.gemmi_groups)} vs "
                f"atomworks={format_groups(result.atomworks_groups)}"
            )
        if len(failed) > 10:
            console.print(f"  ... and {len(failed) - 10} more")

    # Show errors (first 5) if not in verbose mode
    if errors and not verbose:
        console.print()
        console.print("[red]Errors (showing first 5):[/red]")
        for result in errors[:5]:
            console.print(f"  {result.pdb_id}: {result.error}")
        if len(errors) > 5:
            console.print(f"  ... and {len(errors) - 5} more")

    logger.info(f"Summary: {len(passed)} passed, {len(failed)} failed, {len(errors)} errors")

    if output:
        output_data = {
            "passed": passed,
            "failed": [
                {
                    "pdb_id": r.pdb_id,
                    "gemmi_groups": [list(g) for g in r.gemmi_groups],
                    "atomworks_groups": [list(g) for g in r.atomworks_groups],
                }
                for r in failed
            ],
            "errors": [{"pdb_id": r.pdb_id, "error": r.error} for r in errors],
        }
        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)
        console.print(f"\nResults written to [cyan]{output}[/cyan]")
        logger.info(f"Results written to {output}")

    if failed or errors:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
