"""Run the frozen two-flank length baseline on an input-only JSONL bundle."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from .flank_length import predict_case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("predictions", type=Path)
    args = parser.parse_args()
    rows = [json.loads(line) for line in args.inputs.read_text().splitlines() if line]
    if len({row["case_id"] for row in rows}) != len(rows):
        raise ValueError("duplicate case_id")
    t0 = time.perf_counter()
    predictions = [predict_case(row) for row in rows]
    elapsed = time.perf_counter() - t0
    args.predictions.parent.mkdir(parents=True, exist_ok=True)
    args.predictions.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in predictions)
    )
    print(json.dumps({"cases": len(rows), "elapsed_seconds": elapsed}, sort_keys=True))


if __name__ == "__main__":
    main()
