#!/usr/bin/env python3
"""Generate the deterministic unitFinder installation/interface smoke input."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.unitfinder.smoke import generate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    generate(args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
