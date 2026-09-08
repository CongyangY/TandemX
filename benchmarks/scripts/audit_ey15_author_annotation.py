#!/usr/bin/env python3
"""Audit final-assembly TandemX arrays against the Ey15-2 author GFF3."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.retrospective.annotation_overlap import DEFAULT_REPEAT_CLASSES, run_annotation_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arrays", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--annotation-class",
        action="append",
        dest="annotation_classes",
        help="GFF3 feature type to include; repeat for multiple classes",
    )
    args = parser.parse_args()
    run_annotation_audit(
        args.arrays,
        args.annotation,
        args.outdir,
        tuple(args.annotation_classes) if args.annotation_classes else DEFAULT_REPEAT_CLASSES,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
