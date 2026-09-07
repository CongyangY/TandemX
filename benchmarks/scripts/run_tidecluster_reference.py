#!/usr/bin/env python3
"""Run and normalize TideCluster on a verified reference-window sample."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Iterable

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.profile_stage_resources import profile
from benchmarks.scripts.sample_reference_windows import reference_metadata
from benchmarks.tidecluster.normalize import (
    ResolvedTideClusterRecord,
    TideClusterRecord,
    normalize_resolved_tidecluster,
)


NORMALIZED_FIELDS = (
    "sequence_id",
    "start",
    "end",
    "family_id",
    "period",
    "consensus_sequence",
    "copy_number",
    "representative_tidehunter_id",
    "copy_number_source",
    "source",
    "warning",
)


def command_output(command: list[str]) -> str:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return (completed.stdout or completed.stderr).strip()


def load_sample(receipt_path: Path, sample_id: str) -> tuple[Path, dict[str, object]]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    samples = receipt.get("samples", {})
    if receipt.get("complete") is not True or sample_id not in samples:
        raise ValueError(f"Unknown or incomplete reference sample: {sample_id}")
    sample = samples[sample_id]
    path = receipt_path.parent / str(sample["path"])
    if (
        not path.is_file()
        or path.stat().st_size != sample["bytes"]
        or digest_file(path) != sample["sha256"]
    ):
        raise ValueError(f"Reference sample differs from its receipt: {path}")
    return path.resolve(), sample


def union_bases(records: Iterable[TideClusterRecord | ResolvedTideClusterRecord]) -> int:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for record in records:
        grouped[record.sequence_id].append((record.start, record.end))
    total = 0
    for intervals in grouped.values():
        current_end = -1
        for start, end in sorted(intervals):
            total += max(0, end - max(start, current_end))
            current_end = max(current_end, end)
    return total


def quantile(values: list[int], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def summarize_records(
    records: list[TideClusterRecord | ResolvedTideClusterRecord],
    sequence_lengths: dict[str, int],
) -> dict[str, object]:
    for record in records:
        length = sequence_lengths.get(record.sequence_id)
        if length is None or not 0 <= record.start < record.end <= length:
            raise ValueError(
                f"TideCluster interval is outside the sampled reference: "
                f"{record.sequence_id}:{record.start}-{record.end}"
            )
    covered = union_bases(records)
    periods = [record.period for record in records]
    total_bases = sum(sequence_lengths.values())
    return {
        "predicted_array_count": len(records),
        "predicted_family_count": len({record.family_id for record in records}),
        "predicted_positive_sequence_count": len({record.sequence_id for record in records}),
        "predicted_union_bp": covered,
        "predicted_union_base_fraction": covered / total_bases,
        "period_bp_min": min(periods) if periods else None,
        "period_bp_median": quantile(periods, 0.5),
        "period_bp_p95": quantile(periods, 0.95),
        "period_bp_max": max(periods) if periods else None,
        "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
    }


def normalize_real_outputs(
    tidehunter_gff: Path,
    intermediate_clustering_gff: Path,
    clustering_gff: Path,
    sequence_lengths: dict[str, int],
    outdir: Path,
) -> dict[str, object]:
    records = normalize_resolved_tidecluster(
        tidehunter_gff, intermediate_clustering_gff, clustering_gff
    )
    summary = summarize_records(records, sequence_lengths)
    rows = [
        {
            "sequence_id": record.sequence_id,
            "start": record.start,
            "end": record.end,
            "family_id": record.family_id,
            "period": record.period,
            "consensus_sequence": record.consensus_sequence,
            "copy_number": record.copy_number,
            "representative_tidehunter_id": record.representative_tidehunter_id,
            "copy_number_source": record.copy_number_source,
            "source": "TideCluster_clustering_joined_to_TideHunter",
            "warning": "family_sequence_from_cluster_representative;copy_number_unavailable_for_merged_intervals;descriptive_real_reference_output_without_accuracy_truth",
        }
        for record in records
    ]
    outdir.mkdir(parents=True, exist_ok=False)
    write_table(outdir / "normalized_arrays.tsv", rows, list(NORMALIZED_FIELDS))
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def tool_metadata(source_dir: Path, env_bin: Path) -> dict[str, object]:
    python = env_bin / "python"
    tidecluster = source_dir / "TideCluster.py"
    tools = {
        "python": python,
        "tidecluster": tidecluster,
        "tidehunter": env_bin / "TideHunter",
        "mmseqs": env_bin / "mmseqs",
        "blastn": env_bin / "blastn",
    }
    missing = [str(path) for path in tools.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Missing TideCluster runtime files: {missing}")
    source_files = ("TideCluster.py", "tc_utils.py", "version.py")
    versions = {
        "tidecluster": command_output([str(python), str(tidecluster), "--version"]),
        "tidehunter": command_output([str(tools["tidehunter"]), "-v"]).splitlines()[-1],
        "mmseqs": command_output([str(tools["mmseqs"]), "version"]).splitlines()[-1],
        "blastn": command_output([str(tools["blastn"]), "-version"]).splitlines()[0],
    }
    if versions["tidecluster"] == "1.21.2" and versions["tidehunter"] != "1.4.3":
        raise ValueError(
            "TideCluster 1.21.2 requires TideHunter 1.4.3; use the pinned "
            "container instead of bypassing its version check"
        )
    return {
        "versions": versions,
        "paths": {name: str(path.resolve()) for name, path in tools.items()},
        "sha256": {name: digest_file(path) for name, path in tools.items()},
        "source_sha256": {
            name: digest_file(source_dir / name) for name in source_files
        },
        "tidecluster_source_commit": command_output(
            ["git", "-C", str(source_dir), "rev-parse", "HEAD"]
        ),
        "tidecluster_source_status": command_output(
            ["git", "-C", str(source_dir), "status", "--short"]
        ),
    }


def run(
    sampling_receipt: Path,
    sample_id: str,
    outdir: Path,
    source_dir: Path,
    env_bin: Path,
    cpus: int,
    minimum_length: int,
    interval: float,
) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    if not 1 <= cpus <= 64 or minimum_length <= 0:
        raise ValueError("cpus must be in [1,64] and minimum length must be positive")
    sampling_receipt = sampling_receipt.resolve()
    source_dir = source_dir.resolve()
    env_bin = env_bin.resolve()
    sample_path, sample = load_sample(sampling_receipt, sample_id)
    sequence_lengths, observed_sample_sha256 = reference_metadata(sample_path)
    if (
        observed_sample_sha256 != sample["sha256"]
        or sum(sequence_lengths.values()) != sample["total_bases"]
        or len(sequence_lengths) != sample["window_count"]
    ):
        raise ValueError("Sample FASTA content differs from declared sampling metrics")

    outdir.mkdir(parents=True)
    profile_dir = outdir / "profile"
    prefix = outdir / "tc"
    metadata = tool_metadata(source_dir, env_bin)
    python = str(env_bin / "python")
    tidecluster = str(source_dir / "TideCluster.py")
    base = ["-c", str(cpus), "-pr", str(prefix), "-f", str(sample_path)]
    stages = [
        {
            "name": "tidehunter",
            "command": [python, tidecluster, "tidehunter", *base],
            "cwd": str(outdir),
            "scratch_dir": str(outdir),
        },
        {
            "name": "clustering",
            "command": [
                python,
                tidecluster,
                "clustering",
                "-m",
                str(minimum_length),
                *base,
            ],
            "cwd": str(outdir),
            "scratch_dir": str(outdir),
        },
    ]
    stage_manifest = outdir / "stage_manifest.json"
    stage_manifest.write_text(json.dumps({"stages": stages}, indent=2) + "\n")
    project_root = Path(__file__).resolve().parents[2]
    project_commit = command_output(["git", "-C", str(project_root), "rev-parse", "HEAD"])
    environment = {
        "schema_version": 1,
        "complete": False,
        "scope": "descriptive_real_reference_window_comparator_without_accuracy_truth",
        "sample_id": sample_id,
        "sample": sample,
        "sampling_receipt": str(sampling_receipt),
        "sampling_receipt_sha256": digest_file(sampling_receipt),
        "sample_sequence_lengths": sequence_lengths,
        "tool": metadata,
        "project_commit": project_commit,
        "runner_sha256": digest_file(Path(__file__)),
        "normalizer_sha256": digest_file(
            project_root / "benchmarks/tidecluster/normalize.py"
        ),
        "profiler_sha256": digest_file(
            project_root / "benchmarks/scripts/profile_stage_resources.py"
        ),
        "cpus": cpus,
        "minimum_length": minimum_length,
        "resource_scope": "aggregate_live_process_tree_plus_scratch_high_water",
        "warning": "sampled_real_reference_no_independent_accuracy_truth_no_whole_genome_context",
    }
    environment_path = outdir / "environment.json"
    environment_path.write_text(json.dumps(environment, indent=2) + "\n")
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    for path in (
        Path(__file__),
        project_root / "benchmarks/tidecluster/normalize.py",
        project_root / "benchmarks/scripts/profile_stage_resources.py",
        project_root / "benchmarks/scripts/sample_reference_windows.py",
    ):
        shutil.copyfile(path, snapshot / path.name)

    old_path = os.environ.get("PATH")
    os.environ["PATH"] = str(env_bin) + os.pathsep + (old_path or "")
    try:
        profile_exit = profile(stage_manifest, profile_dir, interval)
    finally:
        if old_path is None:
            os.environ.pop("PATH", None)
        else:
            os.environ["PATH"] = old_path
    profile_receipt = json.loads((profile_dir / "receipt.json").read_text())
    if profile_exit != 0 or profile_receipt.get("complete") is not True:
        (outdir / "run_receipt.json").write_text(
            json.dumps({"complete": False, "profile": profile_receipt}, indent=2) + "\n"
        )
        raise RuntimeError("TideCluster stage failed; inspect retained profile logs")

    normalized_dir = outdir / "normalized"
    summary = normalize_real_outputs(
        Path(str(prefix) + "_tidehunter.gff3"),
        Path(str(prefix) + "_clustering.gff3_1.gff3"),
        Path(str(prefix) + "_clustering.gff3"),
        sequence_lengths,
        normalized_dir,
    )
    stage_rows = (profile_dir / "stages.tsv").read_text().splitlines()
    result = {
        "complete": True,
        "sample_id": sample_id,
        "summary": summary,
        "profile_stage_count": len(stage_rows) - 1,
        "accuracy": "not_assessed_without_independent_real_array_and_family_truth",
        "warning": "sampled_real_reference_no_whole_genome_context",
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    environment["complete"] = True
    environment_path.write_text(json.dumps(environment, indent=2) + "\n")
    tracked_outputs = (
        "environment.json",
        "stage_manifest.json",
        "profile/receipt.json",
        "profile/stages.tsv",
        "profile/samples.tsv",
        "profile/tidehunter.stdout.log",
        "profile/tidehunter.stderr.log",
        "profile/clustering.stdout.log",
        "profile/clustering.stderr.log",
        "tc_tidehunter.gff3",
        "tc_clustering.gff3",
        "tc_clustering.gff3_1.gff3",
        "tc_cmd_args.json",
        "tc_consensus/consensus_sequences_all.fasta",
        "normalized/normalized_arrays.tsv",
        "normalized/summary.json",
        "result.json",
    )
    missing = [name for name in tracked_outputs if not (outdir / name).is_file()]
    if missing:
        raise ValueError(f"TideCluster completed but required outputs are missing: {missing}")
    receipt = {
        "complete": True,
        "files": [
            {
                "file": name,
                "bytes": (outdir / name).stat().st_size,
                "sha256": digest_file(outdir / name),
            }
            for name in tracked_outputs
        ],
    }
    (outdir / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampling-receipt", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--tidecluster-source", required=True, type=Path)
    parser.add_argument("--env-bin", required=True, type=Path)
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--minimum-length", type=int, default=100)
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    try:
        run(
            args.sampling_receipt,
            args.sample_id,
            args.outdir,
            args.tidecluster_source,
            args.env_bin,
            args.cpus,
            args.minimum_length,
            args.interval,
        )
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
