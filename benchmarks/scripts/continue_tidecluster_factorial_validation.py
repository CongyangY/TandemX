#!/usr/bin/env python3
"""Continue only unattempted cells after a frozen TideCluster factorial failure."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.tidecluster.factorial_continue import run_continuation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--continuation-config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--docker", default="docker")
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    run_continuation(
        args.base_config,
        args.continuation_config,
        args.outdir,
        args.docker,
        args.interval,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
