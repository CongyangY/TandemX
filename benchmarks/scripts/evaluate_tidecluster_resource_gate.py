#!/usr/bin/env python3
"""Record whether the frozen 1-Gb TideCluster resource gate permits a run."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

from benchmarks.tidecluster.resource_gate import evaluate_resource_gate, write_gate_receipt


def command_output(command: list[str]) -> str:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return completed.stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--docker", default="docker")
    args = parser.parse_args()
    runtime_memory_bytes = int(
        command_output([args.docker, "info", "--format", "{{.MemTotal}}"])
    )
    project_root = Path(__file__).resolve().parents[2]
    project_commit = command_output(["git", "-C", str(project_root), "rev-parse", "HEAD"])
    receipt = evaluate_resource_gate(args.config, runtime_memory_bytes, project_commit)
    write_gate_receipt(args.outdir / "resource_gate_receipt.json", receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
