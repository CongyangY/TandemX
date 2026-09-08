#!/usr/bin/env python3
"""Audit a clean checkout of the frozen official unitFinder source."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.unitfinder.source_audit import audit_source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    audit_source(args.source_dir, args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
