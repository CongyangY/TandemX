#!/usr/bin/env python3
"""Compare two TandemX run directories for parameter consistency."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from tandemx.run_compare import compare_run_directories


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare two TandemX discover or pipeline output directories and write "
            "compare_runs.tsv plus compare_runs.md. This is a post hoc consistency "
            "check; it does not rerun discovery."
        )
    )
    parser.add_argument("--run-a", required=True, type=Path, help="First TandemX run or discover output directory.")
    parser.add_argument("--run-b", required=True, type=Path, help="Second TandemX run or discover output directory.")
    parser.add_argument("--outdir", required=True, type=Path, help="Directory for compare_runs.tsv and compare_runs.md.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = compare_run_directories(args.run_a, args.run_b, args.outdir)
    direct_row = next(row for row in rows if row["item"] == "directly_comparable")
    print(f"wrote {args.outdir / 'compare_runs.tsv'}")
    print(f"wrote {args.outdir / 'compare_runs.md'}")
    print(f"directly_comparable={direct_row['run_a_value']}")
    if direct_row["run_a_value"] == "false":
        print(direct_row["notes"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
