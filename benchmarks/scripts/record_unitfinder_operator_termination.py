#!/usr/bin/env python3
"""Record a post-start operator termination without recasting it as tool failure."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def fasta_counts(path: Path) -> dict[str, int]:
    records = bases = 0
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                records += 1
            elif line:
                bases += len(line)
    return {"records": records, "bases": bases}


def record(
    result_dir: Path,
    config_path: Path,
    cap_seconds: float,
    stop_command: str,
) -> dict[str, Any]:
    output = result_dir / "operator_termination.json"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite termination record: {output}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt_path = result_dir / "run_receipt.json"
    profile_path = result_dir / "profile/receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    with (result_dir / "profile/stages.tsv").open(encoding="utf-8", newline="") as handle:
        stages = list(csv.DictReader(handle, delimiter="\t"))
    smoke = next((row for row in stages if row["stage"] == "smoke"), None)
    if (
        receipt.get("complete") is not False
        or receipt.get("fate") != "external_process_failure"
        or receipt.get("accuracy_available") is not False
        or profile.get("complete") is not False
        or smoke is None
        or int(smoke["exit_code"]) != 137
        or float(smoke["wall_seconds"]) < cap_seconds
    ):
        raise ValueError("result is not the expected operator-stopped smoke")
    required: dict[str, Any] = {}
    for name in config["acceptance"]["required_nonempty_outputs"]:
        path = result_dir / "run" / name
        required[name] = {
            "exists": path.is_file(),
            "bytes": path.stat().st_size if path.is_file() else None,
            "sha256": digest(path) if path.is_file() else None,
        }
    partial_fasta: dict[str, Any] = {}
    for name in ("unitfinder_smoke.Units.fa", "unitfinder_smoke.Consensus.fa"):
        path = result_dir / "run" / name
        if path.is_file() and path.stat().st_size:
            partial_fasta[name] = {
                **fasta_counts(path),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
    result = {
        "schema_version": 1,
        "complete": True,
        "operator_stop_requested": True,
        "stop_command": stop_command,
        "operational_cap_seconds": cap_seconds,
        "cap_declared": "after_run_start_before_cap_in_agent_commentary",
        "formal_benchmark_timeout": False,
        "termination_classification": (
            "post_start_host_operational_termination_not_formal_tool_timeout"
        ),
        "smoke_stage_wall_seconds": float(smoke["wall_seconds"]),
        "smoke_stage_exit_code": int(smoke["exit_code"]),
        "run_receipt_sha256": digest(receipt_path),
        "profile_receipt_sha256": digest(profile_path),
        "config_sha256": digest(config_path),
        "required_outputs": required,
        "partial_fasta": partial_fasta,
        "accuracy_available": False,
        "oom_interpretation": (
            "not_OOM_evidence_exit_137_followed_operator_docker_stop"
        ),
        "boundary": (
            "The operational cap was not present in the frozen v3 config. This record "
            "must not be cited as a formal unitFinder timeout, OOM, algorithm failure "
            "or zero-accuracy result."
        ),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--cap-seconds", required=True, type=float)
    parser.add_argument("--stop-command", required=True)
    args = parser.parse_args()
    record(args.result_dir, args.config, args.cap_seconds, args.stop_command)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
