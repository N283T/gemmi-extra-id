from __future__ import annotations

import argparse
from pathlib import Path

from atomworks.io import parse
from atomworks.io.utils.io_utils import to_cif_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate AtomWorks reference mmCIF with molecule_id."
    )
    parser.add_argument("input", nargs="?", default="data/148L.cif")
    parser.add_argument("output", nargs="?", default="data/148L_atomworks_reference.cif")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = parse(str(Path(args.input)))
    to_cif_file(
        result["asym_unit"],
        str(Path(args.output)),
        include_nan_coords=False,
        extra_fields=["molecule_id"],
    )


if __name__ == "__main__":
    main()
