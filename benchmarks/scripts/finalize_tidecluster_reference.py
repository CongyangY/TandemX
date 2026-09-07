#!/usr/bin/env python3
"""Finalize completed TideCluster stages without rerunning external tools."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.run_tidecluster_docker_reference import parse_gnu_time
from benchmarks.scripts.run_tidecluster_reference import (
    command_output,
    load_sample,
    normalize_real_outputs,
)
from benchmarks.scripts.sample_reference_windows import reference_metadata


EXTERNAL_OUTPUTS = (
    "tc_tidehunter.gff3",
    "tc_clustering.gff3_1.gff3",
    "tc_clustering.gff3",
    "tc_cmd_args.json",
    "tc_consensus/consensus_sequences_all.fasta",
)


def validate_completed_stages(run_dir: Path) -> dict[str, object]:
    profile_dir = run_dir / "profile"
    receipt = json.loads((profile_dir / "receipt.json").read_text())
    with (profile_dir / "stages.tsv").open(newline="", encoding="utf-8") as handle:
        stages = list(csv.DictReader(handle, delimiter="\t"))
    if (
        receipt.get("complete") is not True
        or receipt.get("requested_stage_count") != 2
        or receipt.get("completed_stage_count") != 2
        or {row["stage"] for row in stages} != {"tidehunter", "clustering"}
        or any(int(row["exit_code"]) != 0 for row in stages)
    ):
        raise ValueError("TideCluster external stages are not complete")
    missing = [name for name in EXTERNAL_OUTPUTS if not (run_dir / name).is_file()]
    if missing:
        raise ValueError(f"Completed TideCluster run lacks outputs: {missing}")
    internal = {
        stage: parse_gnu_time(run_dir / f"{stage}.gnu_time.txt")
        for stage in ("tidehunter", "clustering")
    }
    return {"profile": receipt, "stages": stages, "internal_gnu_time": internal}


def finalize(
    run_dir: Path,
    sampling_receipt: Path,
    sample_id: str,
) -> dict[str, object]:
    run_dir = run_dir.resolve()
    if not run_dir.is_dir():
        raise ValueError(f"TideCluster run directory is missing: {run_dir}")
    if (run_dir / "finalization_receipt.json").exists() or (run_dir / "normalized").exists():
        raise ValueError("TideCluster run is already finalized or has ambiguous normalized output")
    stage_evidence = validate_completed_stages(run_dir)
    sampling_receipt = sampling_receipt.resolve()
    sample_path, sample = load_sample(sampling_receipt, sample_id)
    sequence_lengths, observed_sha256 = reference_metadata(sample_path)
    if (
        observed_sha256 != sample["sha256"]
        or sum(sequence_lengths.values()) != sample["total_bases"]
        or len(sequence_lengths) != sample["window_count"]
    ):
        raise ValueError("Sample FASTA differs from its receipt during finalization")

    project_root = Path(__file__).resolve().parents[2]
    external_hashes_before = {
        name: digest_file(run_dir / name) for name in EXTERNAL_OUTPUTS
    }
    summary = normalize_real_outputs(
        run_dir / "tc_tidehunter.gff3",
        run_dir / "tc_clustering.gff3_1.gff3",
        run_dir / "tc_clustering.gff3",
        sequence_lengths,
        run_dir / "normalized",
    )
    external_hashes_after = {
        name: digest_file(run_dir / name) for name in EXTERNAL_OUTPUTS
    }
    if external_hashes_before != external_hashes_after:
        raise ValueError("External TideCluster outputs changed during finalization")

    source_paths = (
        Path(__file__),
        project_root / "benchmarks/tidecluster/normalize.py",
        project_root / "benchmarks/scripts/run_tidecluster_reference.py",
        project_root / "benchmarks/scripts/run_tidecluster_docker_reference.py",
    )
    finalization_environment = {
        "schema_version": 1,
        "complete": True,
        "action": "postprocess_completed_external_stages_without_rerun",
        "sample_id": sample_id,
        "sampling_receipt": str(sampling_receipt),
        "sampling_receipt_sha256": digest_file(sampling_receipt),
        "sample_sha256": observed_sha256,
        "external_output_sha256": external_hashes_before,
        "project_commit": command_output(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"]
        ),
        "project_status": command_output(
            ["git", "-C", str(project_root), "status", "--short"]
        ),
        "source_sha256": {path.name: digest_file(path) for path in source_paths},
        "warning": "normalization_only_external_tidecluster_stages_not_rerun",
    }
    environment_path = run_dir / "finalization_environment.json"
    environment_path.write_text(json.dumps(finalization_environment, indent=2) + "\n")
    snapshot = run_dir / "postprocess_source_snapshot"
    snapshot.mkdir()
    for path in source_paths:
        shutil.copyfile(path, snapshot / path.name)

    result = {
        "complete": True,
        "sample_id": sample_id,
        "summary": summary,
        "internal_gnu_time": stage_evidence["internal_gnu_time"],
        "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
        "warning": "sampled_real_reference_no_whole_genome_context",
    }
    (run_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    tracked = (
        *EXTERNAL_OUTPUTS,
        "environment.json",
        "stage_manifest.json",
        "profile/receipt.json",
        "profile/stages.tsv",
        "profile/samples.tsv",
        "profile/tidehunter.stdout.log",
        "profile/tidehunter.stderr.log",
        "profile/clustering.stdout.log",
        "profile/clustering.stderr.log",
        "tidehunter.gnu_time.txt",
        "clustering.gnu_time.txt",
        "normalized/normalized_arrays.tsv",
        "normalized/summary.json",
        "finalization_environment.json",
        "result.json",
    )
    receipt = {
        "complete": True,
        "external_stages_reused": True,
        "files": [
            {
                "file": name,
                "bytes": (run_dir / name).stat().st_size,
                "sha256": digest_file(run_dir / name),
            }
            for name in tracked
        ],
    }
    (run_dir / "finalization_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--sampling-receipt", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()
    try:
        finalize(args.run_dir, args.sampling_receipt, args.sample_id)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
