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
    # molecule_id groups
    gemmi_mol_groups: set[frozenset[str]] = field(default_factory=set)
    atomworks_mol_groups: set[frozenset[str]] = field(default_factory=set)
    # pn_unit_id mappings (chain -> pn_unit_id)
    gemmi_pn_units: dict[str, str] = field(default_factory=dict)
    atomworks_pn_units: dict[str, str] = field(default_factory=dict)
    # Entity equivalence groups (set of chains with same entity value)
    gemmi_chain_entity_groups: set[frozenset[str]] = field(default_factory=set)
    atomworks_chain_entity_groups: set[frozenset[str]] = field(default_factory=set)
    gemmi_pn_unit_entity_groups: set[frozenset[str]] = field(default_factory=set)
    atomworks_pn_unit_entity_groups: set[frozenset[str]] = field(default_factory=set)
    gemmi_mol_entity_groups: set[frozenset[str]] = field(default_factory=set)
    atomworks_mol_entity_groups: set[frozenset[str]] = field(default_factory=set)
    # Failure details
    mol_match: bool = True
    pn_match: bool = True
    chain_entity_match: bool = True
    pn_unit_entity_match: bool = True
    mol_entity_match: bool = True
    error: str | None = None


def find_cif_files(mirror_path: Path) -> list[Path]:
    """Find all CIF files in directory using gemmi.CifWalk."""
    return [Path(cif_str) for cif_str in gemmi.CifWalk(str(mirror_path))]


def to_equivalence_groups(mapping: dict[str, int | str]) -> set[frozenset[str]]:
    """Convert a value mapping to equivalence groups.

    Example:
        {"A": 1, "B": 1, "C": 2} -> {frozenset({"A", "B"}), frozenset({"C"})}
    """
    value_to_keys: dict[int | str, set[str]] = {}
    for key, value in mapping.items():
        if value not in value_to_keys:
            value_to_keys[value] = set()
        value_to_keys[value].add(key)
    return {frozenset(keys) for keys in value_to_keys.values()}


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


@dataclass
class AtomWorksResult:
    """Result from AtomWorks extraction."""

    mol_groups: set[frozenset[str]]
    pn_mapping: dict[str, str]
    chain_entity_groups: set[frozenset[str]]
    pn_unit_entity_groups: set[frozenset[str]]
    mol_entity_groups: set[frozenset[str]]


def extract_atomworks_ids(cif_path: Path) -> AtomWorksResult:
    """Extract molecule grouping, pn_unit_id, and entity fields from AtomWorks.

    Returns:
        AtomWorksResult with molecule groups, pn_unit mapping, and entity groups.
    """
    from atomworks.io.transforms.atom_array import add_id_and_entity_annotations

    with warnings.catch_warnings(), suppress_output():
        warnings.simplefilter("ignore")
        atom_array = load_with_fixed_nonpoly_res_id(cif_path)
        # Use add_id_and_entity_annotations which adds all 6 annotations
        atom_array = add_id_and_entity_annotations(atom_array)

    # chain_id -> molecule_id mapping
    mol_mapping: dict[str, int] = {}
    # chain_id -> pn_unit_id mapping
    pn_mapping: dict[str, str] = {}
    # Entity mappings (chain_id -> entity value)
    chain_entity_mapping: dict[str, int] = {}
    pn_unit_entity_mapping: dict[str, int] = {}
    mol_entity_mapping: dict[str, int] = {}

    for i in range(len(atom_array)):
        chain = str(atom_array.chain_id[i])
        if chain not in mol_mapping:
            mol_mapping[chain] = int(atom_array.molecule_id[i])
            pn_mapping[chain] = str(atom_array.pn_unit_id[i])
            chain_entity_mapping[chain] = int(atom_array.chain_entity[i])
            pn_unit_entity_mapping[chain] = int(atom_array.pn_unit_entity[i])
            mol_entity_mapping[chain] = int(atom_array.molecule_entity[i])

    # molecule_id -> set of chains
    groups: dict[int, set[str]] = {}
    for chain, mol_id in mol_mapping.items():
        if mol_id not in groups:
            groups[mol_id] = set()
        groups[mol_id].add(chain)

    return AtomWorksResult(
        mol_groups={frozenset(g) for g in groups.values()},
        pn_mapping=pn_mapping,
        chain_entity_groups=to_equivalence_groups(chain_entity_mapping),
        pn_unit_entity_groups=to_equivalence_groups(pn_unit_entity_mapping),
        mol_entity_groups=to_equivalence_groups(mol_entity_mapping),
    )


@dataclass
class GemmiResult:
    """Result from gemmi-extra-id extraction."""

    mol_groups: set[frozenset[str]]
    pn_mapping: dict[str, str]
    chain_entity_groups: set[frozenset[str]]
    pn_unit_entity_groups: set[frozenset[str]]
    mol_entity_groups: set[frozenset[str]]


def extract_gemmi_ids(cif_path: Path) -> GemmiResult:
    """Extract molecule grouping, pn_unit_id, and entity fields from gemmi-extra-id.

    Returns:
        GemmiResult with molecule groups, pn_unit mapping, and entity groups.
    """
    from gemmi_extra_id import assign_extended_ids

    result = assign_extended_ids(cif_path)

    # molecule_id -> set of chains
    groups: dict[int, set[str]] = {}
    # chain_id -> pn_unit_id
    pn_mapping: dict[str, str] = {}
    # Entity mappings (use entity_id for chain_entity, pn_unit_entity, molecule_entity)
    chain_entity_mapping: dict[str, str] = {}
    pn_unit_entity_mapping: dict[str, str] = {}
    mol_entity_mapping: dict[str, str] = {}

    for chain, info in result.chain_info.items():
        mol_id = info.molecule_id
        if mol_id not in groups:
            groups[mol_id] = set()
        groups[mol_id].add(chain)
        pn_mapping[chain] = info.pn_unit_id
        # Use entity_id as chain_entity (CIF entity)
        chain_entity_mapping[chain] = info.entity_id
        pn_unit_entity_mapping[chain] = info.pn_unit_entity
        mol_entity_mapping[chain] = info.molecule_entity

    return GemmiResult(
        mol_groups={frozenset(g) for g in groups.values()},
        pn_mapping=pn_mapping,
        chain_entity_groups=to_equivalence_groups(chain_entity_mapping),
        pn_unit_entity_groups=to_equivalence_groups(pn_unit_entity_mapping),
        mol_entity_groups=to_equivalence_groups(mol_entity_mapping),
    )


def compare_one(cif_path: Path) -> CompareResult:
    """Compare molecule_id, pn_unit_id, and entity fields for a single CIF file."""
    pdb_id = cif_path.stem.replace(".cif", "")
    try:
        gemmi_result = extract_gemmi_ids(cif_path)
        atomworks_result = extract_atomworks_ids(cif_path)

        # Initialize all matches as True
        mol_match = True
        pn_match = True
        chain_entity_match = True
        pn_unit_entity_match = True
        mol_entity_match = True

        # Determine what to compare based on global mode
        if _compare_mode == "molecule":
            mol_match = gemmi_result.mol_groups == atomworks_result.mol_groups
        elif _compare_mode == "pn_unit":
            pn_match = gemmi_result.pn_mapping == atomworks_result.pn_mapping
        elif _compare_mode == "entity":
            # Compare entity equivalence relations only
            chain_entity_match = (
                gemmi_result.chain_entity_groups == atomworks_result.chain_entity_groups
            )
            pn_unit_entity_match = (
                gemmi_result.pn_unit_entity_groups == atomworks_result.pn_unit_entity_groups
            )
            mol_entity_match = gemmi_result.mol_entity_groups == atomworks_result.mol_entity_groups
        elif _compare_mode == "all":
            mol_match = gemmi_result.mol_groups == atomworks_result.mol_groups
            pn_match = gemmi_result.pn_mapping == atomworks_result.pn_mapping
            chain_entity_match = (
                gemmi_result.chain_entity_groups == atomworks_result.chain_entity_groups
            )
            pn_unit_entity_match = (
                gemmi_result.pn_unit_entity_groups == atomworks_result.pn_unit_entity_groups
            )
            mol_entity_match = gemmi_result.mol_entity_groups == atomworks_result.mol_entity_groups
        else:  # "both" (molecule + pn_unit)
            mol_match = gemmi_result.mol_groups == atomworks_result.mol_groups
            pn_match = gemmi_result.pn_mapping == atomworks_result.pn_mapping

        all_match = (
            mol_match
            and pn_match
            and chain_entity_match
            and pn_unit_entity_match
            and mol_entity_match
        )

        if all_match:
            return CompareResult(pdb_id, "pass")
        else:
            return CompareResult(
                pdb_id,
                "fail",
                gemmi_mol_groups=gemmi_result.mol_groups,
                atomworks_mol_groups=atomworks_result.mol_groups,
                gemmi_pn_units=gemmi_result.pn_mapping,
                atomworks_pn_units=atomworks_result.pn_mapping,
                gemmi_chain_entity_groups=gemmi_result.chain_entity_groups,
                atomworks_chain_entity_groups=atomworks_result.chain_entity_groups,
                gemmi_pn_unit_entity_groups=gemmi_result.pn_unit_entity_groups,
                atomworks_pn_unit_entity_groups=atomworks_result.pn_unit_entity_groups,
                gemmi_mol_entity_groups=gemmi_result.mol_entity_groups,
                atomworks_mol_entity_groups=atomworks_result.mol_entity_groups,
                mol_match=mol_match,
                pn_match=pn_match,
                chain_entity_match=chain_entity_match,
                pn_unit_entity_match=pn_unit_entity_match,
                mol_entity_match=mol_entity_match,
            )
    except Exception as e:
        return CompareResult(pdb_id, "error", error=str(e))


def format_groups(groups: set[frozenset[str]]) -> str:
    """Format molecule groups for display."""
    sorted_groups = sorted([sorted(g) for g in groups], key=lambda x: x[0] if x else "")
    return ", ".join("{" + ",".join(g) + "}" for g in sorted_groups)


def format_pn_diff(gemmi: dict[str, str], atomworks: dict[str, str]) -> str:
    """Format pn_unit_id differences for display."""
    diffs = []
    all_chains = set(gemmi.keys()) | set(atomworks.keys())
    for chain in sorted(all_chains):
        g_pn = gemmi.get(chain, "?")
        a_pn = atomworks.get(chain, "?")
        if g_pn != a_pn:
            diffs.append(f"{chain}:[{g_pn}!={a_pn}]")
    return ", ".join(diffs) if diffs else "(no diff)"


def format_entity_diff(
    gemmi_groups: set[frozenset[str]], atomworks_groups: set[frozenset[str]]
) -> str:
    """Format entity equivalence group differences for display."""
    if gemmi_groups == atomworks_groups:
        return "(match)"
    return f"gemmi={format_groups(gemmi_groups)} vs aw={format_groups(atomworks_groups)}"


def format_failure(result: CompareResult) -> str:
    """Format failure details for display."""
    parts = []
    if not result.mol_match:
        parts.append(
            f"mol: gemmi={format_groups(result.gemmi_mol_groups)} "
            f"vs aw={format_groups(result.atomworks_mol_groups)}"
        )
    if not result.pn_match:
        parts.append(f"pn: {format_pn_diff(result.gemmi_pn_units, result.atomworks_pn_units)}")
    if not result.chain_entity_match:
        diff = format_entity_diff(
            result.gemmi_chain_entity_groups, result.atomworks_chain_entity_groups
        )
        parts.append(f"chain_entity: {diff}")
    if not result.pn_unit_entity_match:
        diff = format_entity_diff(
            result.gemmi_pn_unit_entity_groups, result.atomworks_pn_unit_entity_groups
        )
        parts.append(f"pn_unit_entity: {diff}")
    if not result.mol_entity_match:
        diff = format_entity_diff(
            result.gemmi_mol_entity_groups, result.atomworks_mol_entity_groups
        )
        parts.append(f"mol_entity: {diff}")
    return "; ".join(parts)


app = typer.Typer(help="Compare gemmi-extra-id with AtomWorks extra-id assignment.")


# Global variable to control comparison mode (set by CLI)
_compare_mode: str = "both"


def _run_comparison(
    cif_files: list[Path],
    workers: int,
    verbose: bool,
    output: Path | None,
) -> bool:
    """Run comparison on a list of CIF files.

    Returns True if all tests pass, False otherwise.
    """
    total = len(cif_files)
    console.print(
        f"Processing [green]{total}[/green] files with [yellow]{workers}[/yellow] workers..."
    )
    logger.info(f"Processing {total} files with {workers} workers...")

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
                failure_msg = format_failure(result)
                logger.warning(f"FAIL: {result.pdb_id} - {failure_msg}")
                if verbose:
                    tqdm.write(f"FAIL: {result.pdb_id}: {failure_msg}")
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
            console.print(f"  {result.pdb_id}: {format_failure(result)}")
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
                    "mol_match": r.mol_match,
                    "pn_match": r.pn_match,
                    "chain_entity_match": r.chain_entity_match,
                    "pn_unit_entity_match": r.pn_unit_entity_match,
                    "mol_entity_match": r.mol_entity_match,
                    "gemmi_mol_groups": [list(g) for g in r.gemmi_mol_groups],
                    "atomworks_mol_groups": [list(g) for g in r.atomworks_mol_groups],
                    "gemmi_pn_units": r.gemmi_pn_units,
                    "atomworks_pn_units": r.atomworks_pn_units,
                    "gemmi_chain_entity_groups": [list(g) for g in r.gemmi_chain_entity_groups],
                    "atomworks_chain_entity_groups": [
                        list(g) for g in r.atomworks_chain_entity_groups
                    ],
                    "gemmi_pn_unit_entity_groups": [list(g) for g in r.gemmi_pn_unit_entity_groups],
                    "atomworks_pn_unit_entity_groups": [
                        list(g) for g in r.atomworks_pn_unit_entity_groups
                    ],
                    "gemmi_mol_entity_groups": [list(g) for g in r.gemmi_mol_entity_groups],
                    "atomworks_mol_entity_groups": [list(g) for g in r.atomworks_mol_entity_groups],
                }
                for r in failed
            ],
            "errors": [{"pdb_id": r.pdb_id, "error": r.error} for r in errors],
        }
        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)
        console.print(f"\nResults written to [cyan]{output}[/cyan]")
        logger.info(f"Results written to {output}")

    return not (failed or errors)


def _validate_compare_mode(compare: str) -> None:
    """Validate and set compare mode."""
    global _compare_mode
    valid_modes = ("molecule", "pn_unit", "both", "entity", "all")
    if compare not in valid_modes:
        console.print(f"[red]Invalid --compare value: {compare}[/red]")
        console.print(f"Valid options: {', '.join(valid_modes)}")
        raise typer.Exit(code=1)
    _compare_mode = compare


@app.command("all")
def cmd_all(
    pdb_mirror: Annotated[Path, typer.Argument(help="Path to PDB mmCIF mirror directory")],
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
    compare: Annotated[
        str,
        typer.Option("--compare", "-c", help="molecule|pn_unit|both|entity|all"),
    ] = "molecule",
) -> None:
    """Run comparison on entire PDB mirror."""
    _validate_compare_mode(compare)
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
    console.print(f"Found [green]{total}[/green] CIF files{ignore_msg}")

    if not _run_comparison(cif_files, workers, verbose, output):
        raise typer.Exit(code=1)


@app.command("subset")
def cmd_subset(
    test_dir: Annotated[Path, typer.Argument(help="Directory containing test CIF files")],
    random: Annotated[
        int, typer.Option("--random", "-r", help="Add N random entries from PDB mirror")
    ] = 0,
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Output JSON file for results")
    ] = None,
    workers: Annotated[
        int, typer.Option("--workers", "-w", help="Number of parallel workers")
    ] = os.cpu_count() or 4,
    log: Annotated[Path | None, typer.Option("--log", "-l", help="Log file path")] = None,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Show details of failed cases")
    ] = False,
    include_ignored: Annotated[
        bool,
        typer.Option("--include-ignored", help="Include entries in IGNORE_LIST"),
    ] = False,
    compare: Annotated[
        str,
        typer.Option("--compare", "-c", help="molecule|pn_unit|both|entity|all"),
    ] = "molecule",
) -> None:
    """Run comparison on a subset of files with optional random additions.

    Takes a directory of specific test cases. Optionally adds random entries
    from the PDB mirror (requires GEMMI_PDB_DIR environment variable or
    ~/.local/share/pdb).
    """
    import random as rand_module

    _validate_compare_mode(compare)
    setup_logging(log)

    # Find CIF files in test directory
    console.print(f"Scanning [cyan]{test_dir}[/cyan] for CIF files...")
    cif_files = find_cif_files(test_dir)

    # Filter out IGNORE_LIST if not included
    if not include_ignored:
        original_count = len(cif_files)
        cif_files = [f for f in cif_files if f.stem.split(".")[0] not in IGNORE_LIST]
        ignored_count = original_count - len(cif_files)
        if ignored_count > 0:
            console.print(f"[dim]Ignored {ignored_count} entries in IGNORE_LIST[/dim]")

    test_pdb_ids = {f.stem.split(".")[0] for f in cif_files}
    console.print(f"Found [green]{len(cif_files)}[/green] test files")

    # Add random entries if requested
    if random > 0:
        # Get all PDB entries from mirror
        all_pdb_files: list[Path] = []
        for pdb_id in test_pdb_ids:
            # Use first entry to find PDB mirror location
            expanded = gemmi.expand_if_pdb_code(pdb_id)
            if expanded:
                pdb_mirror = Path(expanded).parent.parent.parent
                console.print(f"Using PDB mirror: [cyan]{pdb_mirror}[/cyan]")
                all_pdb_files = find_cif_files(pdb_mirror)
                break

        if not all_pdb_files:
            console.print("[red]Could not find PDB mirror location[/red]")
            raise typer.Exit(code=1)

        # Filter out test cases and IGNORE_LIST
        available = [
            f
            for f in all_pdb_files
            if f.stem.split(".")[0] not in test_pdb_ids and f.stem.split(".")[0] not in IGNORE_LIST
        ]

        # Sample random entries
        sample_size = min(random, len(available))
        rand_module.seed(42)  # Reproducible sampling
        random_files = rand_module.sample(available, sample_size)
        console.print(f"Adding [green]{len(random_files)}[/green] random entries")

        cif_files.extend(random_files)

    if not _run_comparison(cif_files, workers, verbose, output):
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
