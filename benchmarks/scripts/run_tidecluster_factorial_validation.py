#!/usr/bin/env python3
"""Run the frozen multi-seed TideCluster assembly-accuracy validation."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.tidecluster.factorial_run import run_validation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--docker", default="docker")
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    run_validation(args.config, args.outdir, args.docker, args.interval)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
