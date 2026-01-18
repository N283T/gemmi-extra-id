#!/usr/bin/env python3
"""Test exact value matching with AtomWorks.

This script tests whether gemmi-extra-id can produce EXACTLY the same
entity values as AtomWorks (not just equivalent groupings).

This is separate from test_hash_values.py which tests equivalence relations.
"""

import contextlib
import io
import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

console = Console()

# PDB entries with known data quality issues
IGNORE_LIST: set[str] = {
    "2g10",  # Chain F: orphan water atoms
    "1ts6",  # Chain C: orphan water atoms
    "2k9y",  # Chains C, D: not in _struct_asym
}


@dataclass
class ExactMatchResult:
    """Result of exact value matching for a single CIF file."""

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
    """Get entity values from AtomWorks."""
    from atomworks.io.parser import parse
    from biotite.structure import AtomArrayStack

    with suppress_output(), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = parse(
            cif_path,
            build_assembly=None,
            add_missing_atoms=True,
            remove_waters=False,
            fix_formal_charges=False,
            add_id_and_entity_annotations=True,
        )

    atom_array = result["asym_unit"]
    if isinstance(atom_array, AtomArrayStack):
        atom_array = atom_array[0]

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

    Uses canonical sequences for polymer chains (current implementation).
    """
    import gemmi

    from gemmi_extra_id.graph import find_components, find_pn_units
    from gemmi_extra_id.graph_hash import compute_hash_entities
    from gemmi_extra_id.mmcif import (
        _get_canonical_sequences,
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

    chain_residues, residue_names = _get_residue_info(atom_site_loop)
    intra_chain_bonds, inter_chain_bonds = _get_residue_bonds(block, covalent_types)

    # Get canonical sequences for polymer chains
    canonical_sequences = _get_canonical_sequences(block)

    # For polymer chains, use canonical sequence
    entity_types = {chain: entity_mapping.get(chain, ("?", "unknown"))[1] for chain in chain_order}
    for chain in chain_order:
        entity_id, entity_type = entity_mapping.get(chain, ("?", "unknown"))
        if entity_type == "polymer" and entity_id in canonical_sequences:
            canon_seq = canonical_sequences[entity_id]
            chain_residues[chain] = [res_key for res_key, _ in canon_seq]
            residue_names[chain] = {res_key: mon_id for res_key, mon_id in canon_seq}
            if len(canon_seq) > 1:
                polymer_bonds = []
                for i in range(len(canon_seq) - 1):
                    res_key_a = canon_seq[i][0]
                    res_key_b = canon_seq[i + 1][0]
                    polymer_bonds.append((res_key_a, res_key_b))
                intra_chain_bonds[chain] = polymer_bonds

    molecule_mapping = find_components(chain_order, edges)
    is_polymer = {chain: entity_types.get(chain, "unknown") == "polymer" for chain in chain_order}
    pn_unit_mapping = find_pn_units(chain_order, edges, entity_types, is_polymer)

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


def _check_exact_match(
    gm_values: dict[str, int],
    aw_values: dict[str, int],
) -> bool:
    """Check if two entity mappings have exactly the same values.

    Compares values for common chains only (AtomWorks may filter some chains).
    """
    gm_keys = set(gm_values.keys())
    aw_keys = {str(k) for k in aw_values}
    common = sorted(gm_keys & aw_keys)

    if not common:
        return True

    aw_lookup = {str(k): v for k, v in aw_values.items()}

    for chain in common:
        if gm_values[chain] != aw_lookup[chain]:
            return False

    return True


def compare_file(cif_path: Path) -> ExactMatchResult:
    """Compare exact values for a single file."""
    pdb_id = cif_path.stem.split(".")[0].lower()

    if pdb_id in IGNORE_LIST:
        return ExactMatchResult(pdb_id=pdb_id, status="skip", error="In ignore list")

    try:
        aw_chain, aw_pn, aw_mol = get_atomworks_entities(cif_path)
        gm_chain, gm_pn, gm_mol = get_gemmi_hash_entities(cif_path)

        # Check exact value match
        chain_match = _check_exact_match(gm_chain, aw_chain)
        pn_match = _check_exact_match(gm_pn, aw_pn)
        mol_match = _check_exact_match(gm_mol, aw_mol)

        status = "pass" if (chain_match and pn_match and mol_match) else "fail"

        return ExactMatchResult(
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
        return ExactMatchResult(pdb_id=pdb_id, status="error", error=str(e))


app = typer.Typer(help="Test exact value matching with AtomWorks.")


@app.command()
def run(
    input_dir: Annotated[Path, typer.Argument(help="Directory containing CIF files")],
    output: Annotated[Path | None, typer.Option("-o", "--output", help="Output JSON file")] = None,
    verbose: Annotated[bool, typer.Option("-v", "--verbose", help="Show detailed output")] = False,
) -> None:
    """Test exact entity value matching between gemmi-extra-id and AtomWorks."""
    import gemmi

    cif_files = [Path(f) for f in gemmi.CifWalk(str(input_dir))]
    console.print(f"Found [cyan]{len(cif_files)}[/cyan] CIF files")

    results: list[ExactMatchResult] = []
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

    console.print()
    table = Table(title="Exact Value Matching Results")
    table.add_column("Status", style="bold")
    table.add_column("Count", justify="right")
    table.add_row("[green]Passed[/green]", str(passed))
    table.add_row("[red]Failed[/red]", str(failed))
    table.add_row("[yellow]Errors[/yellow]", str(errors))
    table.add_row("[dim]Skipped[/dim]", str(skipped))
    table.add_row("Total", str(len(results)))
    console.print(table)

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
