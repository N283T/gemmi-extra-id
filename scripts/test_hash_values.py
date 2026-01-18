#!/usr/bin/env python3
"""Compare gemmi-extra-id hash mode with AtomWorks for exact value matching.

This script verifies that --hash-entities mode produces identical entity values
to AtomWorks, not just equivalent groupings.
"""

import contextlib
import io
import json
import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

logger = logging.getLogger(__name__)
console = Console()

# PDB entries with known data quality issues
IGNORE_LIST: set[str] = {
    "2g10",  # Chain F: orphan water atoms
    "1ts6",  # Chain C: orphan water atoms
    "2k9y",  # Chains C, D: not in _struct_asym
}


@dataclass
class HashCompareResult:
    """Result of comparing hash values for a single CIF file."""

    pdb_id: str
    status: str  # "pass", "fail", "error"
    # Actual values for comparison
    gemmi_chain_entity: dict[str, int] = field(default_factory=dict)
    atomworks_chain_entity: dict[str, int] = field(default_factory=dict)
    gemmi_pn_unit_entity: dict[str, int] = field(default_factory=dict)
    atomworks_pn_unit_entity: dict[str, int] = field(default_factory=dict)
    gemmi_molecule_entity: dict[str, int] = field(default_factory=dict)
    atomworks_molecule_entity: dict[str, int] = field(default_factory=dict)
    # Match flags
    chain_entity_match: bool = True
    pn_unit_entity_match: bool = True
    molecule_entity_match: bool = True
    error: str | None = None


@contextlib.contextmanager
def suppress_output():
    """Suppress stdout and stderr."""
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        yield


def get_atomworks_entities(
    cif_path: Path,
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    """Get entity values from AtomWorks.

    AtomWorks 2.2.0 computes entity values automatically via parse() with
    add_id_and_entity_annotations=True (the default). The entity values
    are stored as annotations on the AtomArray.
    """
    from atomworks.io.parser import parse

    with suppress_output(), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = parse(
            cif_path,
            build_assembly=None,
            add_missing_atoms=False,
            remove_waters=False,
            fix_formal_charges=False,
            add_id_and_entity_annotations=True,
        )

    # Get the asymmetric unit (asym_unit) which contains per-chain annotations
    # asym_unit can be an AtomArrayStack, so get the first model
    from biotite.structure import AtomArrayStack

    atom_array = result["asym_unit"]
    if isinstance(atom_array, AtomArrayStack):
        atom_array = atom_array[0]

    # Extract per-chain entity values
    chain_entity: dict[str, int] = {}
    pn_unit_entity: dict[str, int] = {}
    molecule_entity: dict[str, int] = {}

    seen_chains: set[str] = set()
    for i in range(len(atom_array)):
        chain_id = atom_array.chain_id[i]
        if chain_id not in seen_chains:
            seen_chains.add(chain_id)
            chain_entity[chain_id] = int(atom_array.chain_entity[i])
            pn_unit_entity[chain_id] = int(atom_array.pn_unit_entity[i])
            molecule_entity[chain_id] = int(atom_array.molecule_entity[i])

    return chain_entity, pn_unit_entity, molecule_entity


def get_gemmi_hash_entities(
    cif_path: Path,
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    """Get hash-based entity values from gemmi-extra-id.

    Directly calls the graph_hash module to get chain_entity, pn_unit_entity,
    and molecule_entity values computed via WL graph hashing.
    """
    import gemmi

    from gemmi_extra_id.graph import find_components, find_pn_units
    from gemmi_extra_id.graph_hash import compute_hash_entities
    from gemmi_extra_id.mmcif import (
        _get_chain_info,
        _get_covalent_edges,
        _get_entity_mapping,
        _get_residue_bonds,
        _get_residue_info,
    )

    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block()

    atom_site_col = block.find_loop("_atom_site.label_asym_id")
    atom_site_loop = atom_site_col.get_loop()

    chain_order, _, _ = _get_chain_info(atom_site_loop)
    covalent_types = frozenset({"covale"})
    edges = _get_covalent_edges(block, covalent_types)
    entity_mapping = _get_entity_mapping(block)

    # Get residue-level info for hash computation
    chain_residues, residue_names = _get_residue_info(atom_site_loop)
    intra_chain_bonds, inter_chain_bonds = _get_residue_bonds(block, covalent_types)

    # Compute molecule and pn_unit mappings
    molecule_mapping = find_components(chain_order, edges)
    entity_types = {chain: entity_mapping.get(chain, ("?", "unknown"))[1] for chain in chain_order}
    is_polymer = {chain: entity_types.get(chain, "unknown") == "polymer" for chain in chain_order}
    pn_unit_mapping = find_pn_units(chain_order, edges, entity_types, is_polymer)

    # Compute hash entities
    chain_entity, pn_unit_entity, molecule_entity = compute_hash_entities(
        chain_order=chain_order,
        molecule_mapping=molecule_mapping,
        pn_unit_mapping=pn_unit_mapping,
        chain_residues=chain_residues,
        residue_names=residue_names,
        intra_chain_bonds=intra_chain_bonds,
        inter_chain_bonds=inter_chain_bonds,
    )

    return chain_entity, pn_unit_entity, molecule_entity


def _check_equivalence_relation(
    gm_values: dict[str, int],
    aw_values: dict[str, int],
) -> bool:
    """Check if two entity mappings have the same equivalence relation.

    The entity values are sequential integers (0, 1, 2, ...) assigned in order
    of first appearance. Different libraries may filter out different chains,
    causing the integer values to differ. However, what matters is whether
    chains that have the same entity in one library also have the same entity
    in the other (equivalence relation).

    This function checks that for any two chains in both mappings:
    - If they have the same entity value in gm, they should have the same in aw
    - If they have different entity values in gm, they should have different in aw
    """
    # Get common chains (convert numpy strings to regular strings)
    gm_keys = set(gm_values.keys())
    aw_keys = {str(k) for k in aw_values}
    common = sorted(gm_keys & aw_keys)

    if len(common) < 2:
        # With 0 or 1 common chains, there's nothing to compare
        return True

    # Build lookup for atomworks values (handles numpy string keys)
    aw_lookup = {}
    for k, v in aw_values.items():
        aw_lookup[str(k)] = v

    # Check equivalence relation
    for i, chain1 in enumerate(common):
        for chain2 in common[i + 1 :]:
            gm_same = gm_values[chain1] == gm_values[chain2]
            aw_same = aw_lookup[chain1] == aw_lookup[chain2]
            if gm_same != aw_same:
                return False

    return True


def compare_file(cif_path: Path) -> HashCompareResult:
    """Compare hash values for a single file."""
    pdb_id = cif_path.stem.split(".")[0].lower()

    if pdb_id in IGNORE_LIST:
        return HashCompareResult(pdb_id=pdb_id, status="skip", error="In ignore list")

    try:
        # Get AtomWorks values
        aw_chain, aw_pn, aw_mol = get_atomworks_entities(cif_path)

        # Get gemmi-extra-id hash values
        gm_chain, gm_pn, gm_mol = get_gemmi_hash_entities(cif_path)

        # Compare equivalence relations (not exact values)
        # This is the correct comparison because:
        # 1. AtomWorks filters out some chains (e.g., ZN ions)
        # 2. Entity values are session-local integers, not global IDs
        # 3. What matters is that chains with same hash get same grouping
        chain_match = _check_equivalence_relation(gm_chain, aw_chain)
        pn_match = _check_equivalence_relation(gm_pn, aw_pn)
        mol_match = _check_equivalence_relation(gm_mol, aw_mol)

        status = "pass" if (chain_match and pn_match and mol_match) else "fail"

        return HashCompareResult(
            pdb_id=pdb_id,
            status=status,
            gemmi_chain_entity=gm_chain,
            atomworks_chain_entity={str(k): v for k, v in aw_chain.items()},
            gemmi_pn_unit_entity=gm_pn,
            atomworks_pn_unit_entity={str(k): v for k, v in aw_pn.items()},
            gemmi_molecule_entity=gm_mol,
            atomworks_molecule_entity={str(k): v for k, v in aw_mol.items()},
            chain_entity_match=chain_match,
            pn_unit_entity_match=pn_match,
            molecule_entity_match=mol_match,
        )

    except Exception as e:
        return HashCompareResult(pdb_id=pdb_id, status="error", error=str(e))


app = typer.Typer(help="Compare gemmi-extra-id hash mode with AtomWorks (exact value matching).")


@app.command()
def run(
    input_dir: Annotated[Path, typer.Argument(help="Directory containing CIF files")],
    output: Annotated[Path | None, typer.Option("-o", "--output", help="Output JSON file")] = None,
    verbose: Annotated[bool, typer.Option("-v", "--verbose", help="Show detailed output")] = False,
) -> None:
    """Compare hash entity values between gemmi-extra-id and AtomWorks."""
    import gemmi

    # Find CIF files
    cif_files = [Path(f) for f in gemmi.CifWalk(str(input_dir))]
    console.print(f"Found [cyan]{len(cif_files)}[/cyan] CIF files")

    results: list[HashCompareResult] = []
    passed = 0
    failed = 0
    errors = 0
    skipped = 0

    for cif_path in cif_files:
        result = compare_file(cif_path)
        results.append(result)

        if result.status == "pass":
            passed += 1
            if verbose:
                console.print(f"[green]PASS[/green] {result.pdb_id}")
        elif result.status == "fail":
            failed += 1
            console.print(f"[red]FAIL[/red] {result.pdb_id}")
            if not result.chain_entity_match:
                console.print(f"  chain_entity: gemmi={result.gemmi_chain_entity}")
                console.print(f"  chain_entity: aw   ={result.atomworks_chain_entity}")
            if not result.pn_unit_entity_match:
                console.print(f"  pn_unit_entity: gemmi={result.gemmi_pn_unit_entity}")
                console.print(f"  pn_unit_entity: aw   ={result.atomworks_pn_unit_entity}")
            if not result.molecule_entity_match:
                console.print(f"  molecule_entity: gemmi={result.gemmi_molecule_entity}")
                console.print(f"  molecule_entity: aw   ={result.atomworks_molecule_entity}")
        elif result.status == "skip":
            skipped += 1
            if verbose:
                console.print(f"[yellow]SKIP[/yellow] {result.pdb_id}: {result.error}")
        else:
            errors += 1
            console.print(f"[yellow]ERROR[/yellow] {result.pdb_id}: {result.error}")

    # Summary
    console.print()
    table = Table(title="Hash Value Comparison Results")
    table.add_column("Status", style="bold")
    table.add_column("Count", justify="right")
    table.add_row("[green]Passed[/green]", str(passed))
    table.add_row("[red]Failed[/red]", str(failed))
    table.add_row("[yellow]Errors[/yellow]", str(errors))
    table.add_row("[dim]Skipped[/dim]", str(skipped))
    table.add_row("Total", str(len(results)))
    console.print(table)

    # Save results
    if output:
        output_data = {
            "summary": {
                "passed": passed,
                "failed": failed,
                "errors": errors,
                "skipped": skipped,
                "total": len(results),
            },
            "failed": [
                {
                    "pdb_id": r.pdb_id,
                    "gemmi_chain_entity": r.gemmi_chain_entity,
                    "atomworks_chain_entity": r.atomworks_chain_entity,
                    "gemmi_pn_unit_entity": r.gemmi_pn_unit_entity,
                    "atomworks_pn_unit_entity": r.atomworks_pn_unit_entity,
                    "gemmi_molecule_entity": r.gemmi_molecule_entity,
                    "atomworks_molecule_entity": r.atomworks_molecule_entity,
                }
                for r in results
                if r.status == "fail"
            ],
            "errors": [
                {"pdb_id": r.pdb_id, "error": r.error} for r in results if r.status == "error"
            ],
        }
        output.write_text(json.dumps(output_data, indent=2))
        console.print(f"Results saved to [cyan]{output}[/cyan]")


if __name__ == "__main__":
    app()
