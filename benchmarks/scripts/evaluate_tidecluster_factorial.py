#!/usr/bin/env python3
"""Score one frozen TideCluster factorial-assembly validation run."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.tidecluster.factorial_evaluate import evaluate_factorial


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--genome-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    evaluate_factorial(
        args.run_dir / "tc_tidehunter.gff3",
        args.run_dir / "tc_clustering.gff3_1.gff3",
        args.run_dir / "tc_clustering.gff3",
        args.run_dir / "tc_consensus/consensus_sequences_all.fasta",
        args.genome_dir,
        args.outdir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
