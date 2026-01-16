"""Test gemmi-extra-id against AtomWorks on original PDB files."""

import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Optional

import gemmi
import typer
from rich.console import Console
from rich.table import Table
from tqdm import tqdm

logger = logging.getLogger(__name__)
console = Console()


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

    logging.captureWarnings(True)


@dataclass
class TestResult:
    """Result of testing a single CIF file."""

    pdb_id: str
    status: str  # "passed", "failed", "error"
    errors: list[str] | None = None
    exception: str | None = None


def find_cif_files(mirror_path: Path) -> list[Path]:
    """Find all CIF files in directory using gemmi.CifWalk."""
    return [Path(cif_str) for cif_str in gemmi.CifWalk(str(mirror_path))]


def extract_atomworks_mapping(cif_path: Path) -> dict[str, int]:
    """Extract molecule_id mapping from AtomWorks."""
    from atomworks.io import parse

    result = parse(
        filename=cif_path,
        build_assembly=None,
        add_missing_atoms=False,
        remove_waters=False,
    )

    asym = result["asym_unit"][0]

    mapping: dict[str, int] = {}
    for i in range(len(asym)):
        chain = str(asym.chain_id[i])
        if chain not in mapping:
            mapping[chain] = int(asym.molecule_id[i])

    return mapping


def extract_gemmi_ids(cif_path: Path) -> dict[str, int]:
    """Extract molecule_id mapping from gemmi-extra-id."""
    from gemmi_extra_id import assign_extended_ids

    result = assign_extended_ids(cif_path)
    return dict(result.molecule_id_mapping)


def compare_results(
    expected: dict[str, int], actual: dict[str, int]
) -> tuple[bool, list[str]]:
    """Compare molecule groupings between AtomWorks and gemmi-extra-id.

    Compares which chains are grouped together, not the exact molecule_id values.
    AtomWorks may filter out certain chains (e.g., some ions), so we only
    compare chains that exist in both outputs.
    """
    errors = []

    # Get common chains
    common_chains = set(expected.keys()) & set(actual.keys())
    if not common_chains:
        errors.append("No common chains between AtomWorks and gemmi-extra-id")
        return False, errors

    # Build groupings: mol_id -> set of chains
    def build_groups(mapping: dict[str, int], chains: set[str]) -> dict[int, set[str]]:
        groups: dict[int, set[str]] = {}
        for chain in chains:
            if chain in mapping:
                mol_id = mapping[chain]
                if mol_id not in groups:
                    groups[mol_id] = set()
                groups[mol_id].add(chain)
        return groups

    expected_groups = build_groups(expected, common_chains)
    actual_groups = build_groups(actual, common_chains)

    # Convert to sets of frozensets for comparison
    expected_sets = {frozenset(g) for g in expected_groups.values()}
    actual_sets = {frozenset(g) for g in actual_groups.values()}

    if expected_sets != actual_sets:
        errors.append(
            f"Molecule groupings differ: expected {expected_sets}, got {actual_sets}"
        )

    return len(errors) == 0, errors


def process_one_file(cif_path: Path) -> TestResult:
    """Process a single CIF file and return the test result."""
    pdb_id = cif_path.stem.replace(".cif", "")
    try:
        # Get molecule_id from AtomWorks
        expected = extract_atomworks_mapping(cif_path)

        # Get molecule_id from gemmi-extra-id
        actual = extract_gemmi_ids(cif_path)

        # Compare results
        passed, errors = compare_results(expected, actual)

        if passed:
            return TestResult(pdb_id=pdb_id, status="passed")
        else:
            return TestResult(pdb_id=pdb_id, status="failed", errors=errors)
    except Exception as e:
        return TestResult(pdb_id=pdb_id, status="error", exception=str(e))


app = typer.Typer(help="Test gemmi-extra-id against AtomWorks on PDB files.")


@app.command()
def main(
    pdb_mirror: Annotated[
        Path, typer.Argument(help="Path to PDB mmCIF mirror directory")
    ],
    output: Annotated[
        Optional[Path], typer.Option("--output", "-o", help="Output JSON file for results")
    ] = None,
    workers: Annotated[
        int, typer.Option("--workers", "-w", help="Number of parallel workers")
    ] = os.cpu_count() or 4,
    log: Annotated[
        Optional[Path], typer.Option("--log", "-l", help="Log file path")
    ] = None,
    limit: Annotated[
        Optional[int], typer.Option("--limit", "-n", help="Limit number of files to process")
    ] = None,
) -> None:
    """Test gemmi-extra-id against AtomWorks on original PDB mmCIF files.

    Both AtomWorks and gemmi-extra-id process the same original PDB files,
    ensuring they have access to the same struct_conn data.
    """
    setup_logging(log)

    # Find all CIF files
    console.print(f"Scanning [cyan]{pdb_mirror}[/cyan] for CIF files...")
    cif_files = find_cif_files(pdb_mirror)

    if limit:
        cif_files = cif_files[:limit]

    total = len(cif_files)
    console.print(
        f"Found [green]{total}[/green] CIF files, processing with [yellow]{workers}[/yellow] workers..."
    )

    logger.info(f"Found {total} CIF files, processing with {workers} workers...")

    results: dict[str, list] = {"passed": [], "failed": [], "errors": []}

    # Process files in parallel with tqdm progress bar
    with ProcessPoolExecutor(max_workers=workers) as executor:
        with tqdm(total=total, desc="Processing", unit="file") as pbar:
            for result in executor.map(process_one_file, cif_files):
                if result.status == "passed":
                    results["passed"].append(result.pdb_id)
                    logger.info(f"PASS: {result.pdb_id}")
                elif result.status == "failed":
                    results["failed"].append(
                        {"pdb_id": result.pdb_id, "errors": result.errors}
                    )
                    logger.warning(f"FAIL: {result.pdb_id} - {result.errors}")
                    tqdm.write(f"FAIL: {result.pdb_id}")
                else:
                    results["errors"].append(
                        {"pdb_id": result.pdb_id, "error": result.exception}
                    )
                    logger.error(f"ERROR: {result.pdb_id} - {result.exception}")
                    tqdm.write(f"ERROR: {result.pdb_id} - {result.exception}")
                pbar.update(1)

    # Summary table
    table = Table(title="Test Summary")
    table.add_column("Status", style="bold")
    table.add_column("Count", justify="right")
    table.add_row("[green]Passed[/green]", str(len(results["passed"])))
    table.add_row("[yellow]Failed[/yellow]", str(len(results["failed"])))
    table.add_row("[red]Errors[/red]", str(len(results["errors"])))
    console.print(table)

    logger.info(
        f"Summary: {len(results['passed'])} passed, "
        f"{len(results['failed'])} failed, {len(results['errors'])} errors"
    )

    if output:
        with open(output, "w") as f:
            json.dump(results, f, indent=2)
        console.print(f"Results written to [cyan]{output}[/cyan]")
        logger.info(f"Results written to {output}")

    if results["failed"] or results["errors"]:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
