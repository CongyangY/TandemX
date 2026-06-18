#!/usr/bin/env python3
"""Compare known repeats post hoc against a de novo TandemX monomer catalog."""

from __future__ import annotations

import argparse
from pathlib import Path

from tandemx.sensitivity import check_known_repeats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path, help="monomers.fa produced by tandemx discover.")
    parser.add_argument("--known", required=True, type=Path, help="Known repeat FASTA used only after discovery.")
    parser.add_argument("--out", required=True, type=Path, help="Output known_repeat_matches.tsv path.")
    parser.add_argument("--kmer-size", type=int, default=11, help="Exact orientation-aware k-mer size.")
    args = parser.parse_args()
    for path in (args.catalog, args.known):
        if not path.is_file():
            parser.error(f"input file does not exist: {path}")
    try:
        matches = check_known_repeats(
            args.catalog,
            args.known,
            args.out,
            kmer_size=args.kmer_size,
        )
    except ValueError as exc:
        parser.error(str(exc))
    print(f"known-repeat check: wrote {len(matches)} post hoc matches to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
