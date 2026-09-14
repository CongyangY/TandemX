#!/usr/bin/env python3
"""Reconstruct terminal evidence after the frozen historical controller failed.

This script is intentionally postmortem-only. It reads and hashes preserved
artifacts, writes a separate evidence package, and never launches a native tool.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any


EXPECTED_CONFIG_SHA256 = "c2d2d3888c665ff8a099ac7ac08652e2702afba9d72e177cb8ef2479faafd900"
EXPECTED_FROZEN_RUNNER_SHA256 = "a81a4444ab4e5583c4486cfdd2c3f05be8a7843d5e27e37a65575168bb02f36c"
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
FAILURE = "CellFailure:stage_failed:assemble"
METHODS = ("srf_k151", "srf_k101", "competitive_mapping")
FIELDS = (
    "family_id", "method", "reference_proxy_positive", "old_assembly_bp", "new_assembly_bp",
    "state", "native_retained_read_bp", "read_estimated_bp", "old_read_ratio",
    "predicted_proxy_positive", "outcome",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _elapsed_seconds(start: str, end: str) -> float:
    return (datetime.fromisoformat(end) - datetime.fromisoformat(start)).total_seconds()


def verify_stage(receipt: dict[str, Any], name: str, *, success: bool) -> None:
    if receipt.get("stage") != name:
        raise ValueError(f"Wrong stage receipt: expected {name}")
    if success:
        if receipt.get("exit_code") != 0 or receipt.get("timed_out") is not False:
            raise ValueError(f"Expected successful {name} receipt")
    elif receipt.get("exit_code") != -9 or receipt.get("timed_out") is not True:
        raise ValueError("Assemble receipt is not the frozen timeout terminal state")
    for output in receipt.get("outputs", []):
        path = Path(output["path"])
        if not path.is_file() or path.stat().st_size != output["bytes"] or digest(path) != output["sha256"]:
            raise ValueError(f"Stage output binding mismatch: {path}")


def reconstruct(formal: Path, outdir: Path, repository: Path, terminal_stderr: Path) -> None:
    if outdir.exists():
        raise FileExistsError(outdir)
    if (formal / "summary.tsv").exists() or (formal / "completion.json").exists():
        raise ValueError("Formal root already has controller terminal outputs; refuse reconstruction")
    config_path = formal / "run_config.json"
    frozen_runner = formal / "source_snapshot/benchmarks/scripts/run_historical_prioritization_srf_mapping.py"
    if digest(config_path) != EXPECTED_CONFIG_SHA256 or digest(frozen_runner) != EXPECTED_FROZEN_RUNNER_SHA256:
        raise ValueError("Frozen config/source snapshot binding mismatch")
    config = load_json(config_path)
    plan = read_tsv(formal / "plan.tsv")
    expected_plan = [(species["id"], method) for species in config["species"] for method in METHODS]
    if [(row["species"], row["method"]) for row in plan] != expected_plan:
        raise ValueError("Formal six-cell plan differs from the frozen order")

    cell = formal / "ey15_2/srf_k151"
    receipts = {name: load_json(cell / f"{name}.receipt.json") for name in ("count", "dump", "assemble")}
    verify_stage(receipts["count"], "count", success=True)
    verify_stage(receipts["dump"], "dump", success=True)
    verify_stage(receipts["assemble"], "assemble", success=False)
    srf = cell / "srf.fa"
    if srf.stat().st_size != 0 or digest(srf) != EMPTY_SHA256:
        raise ValueError("Timed-out SRF FASTA is not the preserved zero-byte output")
    header_only = cell / "family_prioritization.tsv"
    if len(header_only.read_text().splitlines()) != 1:
        raise ValueError("Expected the serializer-aborted header-only TSV")
    for species, method in expected_plan[1:]:
        if (formal / species / method).exists():
            raise ValueError(f"Later cell unexpectedly exists: {species}/{method}")
    terminal_text = terminal_stderr.read_text()
    required_terminal = "ValueError: dict contains fields not in fieldnames: 'native_retained_read_bp'"
    if required_terminal not in terminal_text:
        raise ValueError("Controller terminal evidence lacks the exact serializer exception")

    source_species = config["species"][0]
    reference_path = repository / source_species["reference_metrics"]["path"]
    if digest(reference_path) != source_species["reference_metrics"]["sha256"]:
        raise ValueError("Frozen reference-metrics binding mismatch")
    eligible = [row for row in read_tsv(reference_path) if row["eligibility"] == "eligible"]
    if len(eligible) != source_species["reference_metrics"]["eligible_records"]:
        raise ValueError("Frozen eligible-family count mismatch")

    outdir.mkdir(parents=True)
    shutil.copyfile(terminal_stderr, outdir / "controller_terminal.stderr.txt")
    family_rows = []
    for row in eligible:
        family_rows.append({
            "family_id": row["family_id"], "method": "srf_k151",
            "reference_proxy_positive": row["reference_state"] == "reference_collapse",
            "old_assembly_bp": float(row["old_assembly_bp"]),
            "new_assembly_bp": float(row["new_assembly_bp"]),
            "state": FAILURE, "native_retained_read_bp": "N/A", "read_estimated_bp": "N/A",
            "old_read_ratio": "N/A", "predicted_proxy_positive": "N/A", "outcome": "N/A",
        })
    write_tsv(outdir / "first_cell_family_terminal.tsv", family_rows, list(FIELDS))

    stage_seconds = sum(float(receipt["runtime_seconds"]) for receipt in receipts.values())
    workflow_seconds = _elapsed_seconds(receipts["count"]["started_utc"], receipts["assemble"]["ended_utc"])
    summary = []
    for index, (species, method) in enumerate(expected_plan, 1):
        if index == 1:
            summary.append({
                "cell_index": index, "species": species, "method": method, "status": "technical_failure",
                "failure": FAILURE, "eligible_families": len(eligible), "available_families": 0,
                "unavailable_families": len(eligible), "zero_supported_eligible_families": 0,
                "TP": 0, "FN": 0, "FP": 0, "TN": 0,
                "native_stages_completed": "count;dump", "native_terminal_stage": "assemble",
                "native_terminal_exit_code": -9, "native_terminal_timed_out": True,
                "native_stage_runtime_seconds": stage_seconds,
                "observed_cell_wall_seconds": workflow_seconds,
                "peak_stage_rss_mib": float(receipts["assemble"]["peak_rss_mib"]),
                "result_origin": "postmortem_reconstruction_from_frozen_receipts",
            })
        else:
            summary.append({
                "cell_index": index, "species": species, "method": method,
                "status": "not_run_prior_failure", "failure": FAILURE,
                "eligible_families": "N/A", "available_families": 0, "unavailable_families": "N/A",
                "zero_supported_eligible_families": "N/A", "TP": 0, "FN": 0, "FP": 0, "TN": 0,
                "native_stages_completed": "", "native_terminal_stage": "not_run",
                "native_terminal_exit_code": "N/A", "native_terminal_timed_out": "N/A",
                "native_stage_runtime_seconds": 0, "observed_cell_wall_seconds": 0,
                "peak_stage_rss_mib": 0,
                "result_origin": "postmortem_reconstruction_from_absent_later_cell_directories",
            })
    summary_fields = list(dict.fromkeys(key for row in summary for key in row))
    write_tsv(outdir / "terminal_summary.tsv", summary, summary_fields)
    write_json(outdir / "terminal_completion.json", {
        "complete": False, "formal_result_state": "technical_failure_not_completed",
        "cells_planned": 6, "cells_completed": 0, "failed_cells": 1, "not_run_cells": 5,
        "stop_reason": FAILURE,
        "native_interpretation": "SRF k151 assembly timed out at the preregistered 7200-second limit",
        "controller_interpretation": "terminal serializer then failed on an omitted field name",
        "rerun_or_tuning_performed": False,
        "proxy_metrics_available": False,
    })

    original_paths = sorted(path for path in formal.rglob("*") if path.is_file())
    inventory = [{"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}
                 for path in original_paths]
    write_json(outdir / "formal_artifact_inventory.json", {"formal_root": str(formal), "files": inventory})
    receipt = {
        "schema_version": 1,
        "scope": "postmortem_terminal_evidence_only_no_native_execution",
        "formal_root": str(formal),
        "frozen_run_config_sha256": digest(config_path),
        "frozen_runner_snapshot_sha256": digest(frozen_runner),
        "count_receipt_sha256": digest(cell / "count.receipt.json"),
        "dump_receipt_sha256": digest(cell / "dump.receipt.json"),
        "assemble_receipt_sha256": digest(cell / "assemble.receipt.json"),
        "controller_terminal_sha256": digest(outdir / "controller_terminal.stderr.txt"),
        "header_only_family_tsv_sha256": digest(header_only),
        "first_cell_state": "technical_failure",
        "later_cell_state": "not_run_prior_failure",
        "native_cells_launched": 1,
        "native_cells_not_run": 5,
        "formal_run_complete": False,
        "native_stage_rerun": False,
        "parameter_tuning": False,
        "preexisting_paf_imported": False,
        "original_formal_files_modified": False,
    }
    write_json(outdir / "reconstruction_receipt.json", receipt)
    artifacts = sorted(path for path in outdir.rglob("*") if path.is_file() and path.name != "archive_manifest.json")
    write_json(outdir / "archive_manifest.json", {"files": [
        {"path": str(path.relative_to(outdir)), "bytes": path.stat().st_size, "sha256": digest(path)}
        for path in artifacts
    ]})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--controller-stderr", type=Path, required=True)
    args = parser.parse_args()
    reconstruct(args.formal_root.resolve(), args.outdir.resolve(), args.repository.resolve(),
                args.controller_stderr.resolve())


if __name__ == "__main__":
    main()
