#!/usr/bin/env python3
"""Render static figures for TandemX discover outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from tandemx.visualize.discover_catalog import DiscoverCatalogConfig, render_discover_catalog_report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate publication-oriented static figures from tandemx discover "
            "outputs: families.tsv, family_similarity.tsv, and monomers.fa."
        )
    )
    parser.add_argument("--families", required=True, type=Path, help="families.tsv from tandemx discover.")
    parser.add_argument(
        "--family-similarity",
        required=True,
        type=Path,
        dest="family_similarity",
        help="family_similarity.tsv from tandemx discover.",
    )
    parser.add_argument("--monomers", required=True, type=Path, help="monomers.fa from tandemx discover.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for figures and summary tables.")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top families to annotate and tabulate.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outputs = render_discover_catalog_report(
        DiscoverCatalogConfig(
            families_tsv=args.families,
            family_similarity_tsv=args.family_similarity,
            monomers_fa=args.monomers,
            outdir=args.outdir,
            top_n=args.top_n,
        )
    )
    for path in outputs:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
