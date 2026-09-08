#!/usr/bin/env python3
"""Independently verify a completed unitFinder interface smoke result."""
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


def fasta_summary(path: Path) -> dict[str, int]:
    records = bases = 0
    current_bases = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise ValueError(f"empty FASTA header: {path}:{line_number}")
                if records and current_bases == 0:
                    raise ValueError(f"empty FASTA record: {path}:{line_number}")
                records += 1
                current_bases = 0
            else:
                if records == 0:
                    raise ValueError(f"sequence before FASTA header: {path}:{line_number}")
                current_bases += len(line)
                bases += len(line)
    if records == 0 or current_bases == 0 or bases == 0:
        raise ValueError(f"empty or incomplete FASTA: {path}")
    return {"records": records, "bases": bases}


def verify(config_path: Path, result_dir: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt = json.loads((result_dir / "run_receipt.json").read_text(encoding="utf-8"))
    environment = json.loads((result_dir / "environment.json").read_text(encoding="utf-8"))
    profile_receipt = json.loads(
        (result_dir / "profile/receipt.json").read_text(encoding="utf-8")
    )
    if (
        receipt.get("complete") is not True
        or receipt.get("fate") != "smoke_passed"
        or receipt.get("accuracy_available") is not False
        or receipt.get("missing_required_outputs") != []
    ):
        raise ValueError("unitFinder receipt is not a successful interface-only smoke")
    if (
        environment.get("complete") is not True
        or environment.get("image", {}).get("image_id") != config["image_id"]
        or environment.get("parent_build_success", {}).get("image_id")
        != config["image_id"]
    ):
        raise ValueError("unitFinder environment/image identity differs")
    if (
        profile_receipt.get("complete") is not True
        or profile_receipt.get("requested_stage_count") != 2
        or profile_receipt.get("completed_stage_count") != 2
    ):
        raise ValueError("unitFinder profile is incomplete")
    with (result_dir / "profile/stages.tsv").open(encoding="utf-8", newline="") as handle:
        stages = list(csv.DictReader(handle, delimiter="\t"))
    if [row["stage"] for row in stages] != ["help", "smoke"] or any(
        row["exit_code"] != "0" for row in stages
    ):
        raise ValueError("unitFinder stage results differ from the frozen contract")

    parent_records = config["parent_build_success"]["artifacts"]
    parent_checks = 0
    for name, expected in parent_records.items():
        path = result_dir / "parent_build_snapshot" / name
        if not path.is_file() or digest(path) != expected:
            raise ValueError(f"unitFinder parent snapshot changed: {name}")
        parent_checks += 1

    fasta: dict[str, dict[str, Any]] = {}
    for name in config["acceptance"]["required_nonempty_outputs"]:
        path = result_dir / "run" / name
        record = receipt.get("outputs", {}).get(name)
        if (
            record is None
            or not path.is_file()
            or path.stat().st_size != record["bytes"]
            or digest(path) != record["sha256"]
        ):
            raise ValueError(f"unitFinder required output changed: {name}")
        fasta[name] = {
            **fasta_summary(path),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }

    result = {
        "schema_version": 1,
        "complete": True,
        "independent_of_runner_implementation": True,
        "accuracy_available": False,
        "image_id": config["image_id"],
        "parent_build_artifact_checks": parent_checks,
        "stage_exit_codes": {row["stage"]: int(row["exit_code"]) for row in stages},
        "required_fasta": fasta,
        "boundary": (
            "This verifier establishes interface completion and artifact integrity only. "
            "It does not compare unitFinder outputs with the planted truth and therefore "
            "provides no accuracy endpoint."
        ),
    }
    (result_dir / "independent_verification.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--result-dir", required=True, type=Path)
    args = parser.parse_args()
    verify(args.config, args.result_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
