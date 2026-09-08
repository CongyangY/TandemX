#!/usr/bin/env python3
"""Independently validate an annotation-overlap audit by byte rasterization."""
from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.retrospective.annotation_overlap import DEFAULT_REPEAT_CLASSES
from benchmarks.retrospective.raster_validate import validate_rasterized_overlap, write_validation_receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arrays", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--annotation-class", action="append", dest="annotation_classes")
    args = parser.parse_args()
    receipt = validate_rasterized_overlap(
        args.arrays,
        args.annotation,
        args.audit_dir,
        tuple(args.annotation_classes) if args.annotation_classes else DEFAULT_REPEAT_CLASSES,
    )
    write_validation_receipt(args.receipt, receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
