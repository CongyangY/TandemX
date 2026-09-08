#!/usr/bin/env python3
"""Evaluate a frozen old/new donor-matched collapse experiment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.retrospective.evaluate import EvaluationConfig, run_evaluation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--copy-number", required=True, type=Path)
    parser.add_argument("--old-arrays", required=True, type=Path)
    parser.add_argument("--new-arrays", required=True, type=Path)
    parser.add_argument("--sensitivity-arrays", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    payload = json.loads(args.config.read_text(encoding="utf-8"))
    evaluation = payload["evaluation"]
    config = EvaluationConfig(
        collapse_threshold=float(evaluation["collapse_threshold"]),
        overexpansion_threshold=float(evaluation["overexpansion_threshold"]),
        primary_min_new_bp=int(evaluation["primary_min_new_bp"]),
        min_new_bp_sensitivity=tuple(int(value) for value in evaluation["min_new_bp_sensitivity"]),
    )
    run_evaluation(
        args.copy_number,
        args.old_arrays,
        args.new_arrays,
        args.outdir,
        config,
        args.sensitivity_arrays,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
