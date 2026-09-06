#!/usr/bin/env python3
"""Install a wheel in a temporary target and exercise the complete toy workflow.

Run with the tandemx-dev interpreter. Existing environment dependencies are used;
the wheel and Rust extension must import from the isolated target, not checkout.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def verify_wheel(wheel: Path, outdir: Path) -> None:
    wheel, outdir = wheel.resolve(), outdir.resolve()
    if not wheel.is_file():
        raise ValueError(f"Missing wheel: {wheel}")
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Choose an empty verification directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    receipt = {"wheel": str(wheel), "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
               "python": sys.executable, "status": "running", "commands": []}
    with tempfile.TemporaryDirectory(prefix="tandemx-wheel-") as temporary:
        stage = Path(temporary)
        target = stage / "installed"
        environment = {**os.environ, "PYTHONPATH": str(target), "PYTHONNOUSERSITE": "1"}

        def execute(name: str, command: list[str]) -> str:
            result = subprocess.run(command, cwd=stage, env=environment, capture_output=True, text=True, timeout=120)
            (outdir / f"{name}.stdout.log").write_text(result.stdout)
            (outdir / f"{name}.stderr.log").write_text(result.stderr)
            receipt["commands"].append({"name": name, "argv": command, "exit_code": result.returncode})
            if result.returncode:
                receipt["status"] = "failed"
                (outdir / "wheel_validation.json").write_text(json.dumps(receipt, indent=2) + "\n")
                raise ValueError(f"Wheel verification failed at {name}; see {outdir}")
            return result.stdout

        execute("install", [sys.executable, "-m", "pip", "install", "--no-index", "--no-deps", "--target", str(target), str(wheel)])
        origin = execute("import", [sys.executable, "-c",
            "import json,tandemx,tandemx._rust_core; print(json.dumps([tandemx.__file__,tandemx._rust_core.__file__]))"])
        paths = json.loads(origin)
        if not all(Path(p).resolve().is_relative_to(target.resolve()) for p in paths):
            raise ValueError(f"Wheel test accidentally imported from outside the isolated target: {paths}")
        receipt["isolated_imports"] = paths
        cli = [sys.executable, "-m", "tandemx.cli"]
        toy = outdir / "toy"
        execute("simulate", [*cli, "simulate", "toy", "--outdir", str(toy)])
        workflow = outdir / "workflow"
        execute("workflow", [*cli, "run", "--reads", str(toy / "reads.fa"), "--assembly", str(toy / "assembly.fa"),
                              "--genome-size", "7744", "--threads", "1", "--kmer-backend", "rust", "--outdir", str(workflow)])
        with (workflow / "pipeline_summary.tsv").open() as handle:
            steps = list(csv.DictReader(handle, delimiter="\t"))
        if len(steps) != 7 or any(int(r["exit_status"]) != 0 or r["output_validated"].lower() != "true" for r in steps):
            raise ValueError("The wheel did not complete and validate all seven toy workflow steps")
        receipt.update({"status": "passed", "validated_steps": [r["step"] for r in steps],
                        "warning": "toy_scale_smoke_test_not_production_or_cross_platform_validation"})
        (outdir / "wheel_validation.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    verify_wheel(args.wheel, args.outdir)
