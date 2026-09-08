#!/usr/bin/env python3
"""Run pinned TideCluster in Docker on a verified reference-window sample."""
from __future__ import annotations

import argparse
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.profile_stage_resources import profile
from benchmarks.scripts.run_tidecluster_reference import (
    command_output,
    load_sample,
    normalize_real_outputs,
)
from benchmarks.scripts.sample_reference_windows import reference_metadata


EXPECTED_VERSIONS = {
    "tidecluster": "1.21.2",
    "tidehunter": "1.4.3",
    "kitehor": "0.13.2",
    "mmseqs": "747c64cc8db3b4803a0f1194a3f75b3ba9f81bcb",
    "blastn": "2.16.0+",
}


def container_path(path: Path, host_root: Path, container_root: str) -> str:
    try:
        relative = path.resolve().relative_to(host_root.resolve())
    except ValueError as error:
        raise ValueError(f"Path is outside the mounted benchmark root: {path}") from error
    return str(PurePosixPath(container_root) / PurePosixPath(relative.as_posix()))


def parse_gnu_time(
    path: Path, *, require_success: bool = True
) -> dict[str, float | int]:
    text = path.read_text(encoding="utf-8")
    patterns = {
        "user_seconds": r"User time \(seconds\): ([0-9.]+)",
        "system_seconds": r"System time \(seconds\): ([0-9.]+)",
        "maximum_rss_kb": r"Maximum resident set size \(kbytes\): ([0-9]+)",
        "exit_status": r"Exit status: ([0-9]+)",
    }
    result: dict[str, float | int] = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        if match is None:
            raise ValueError(f"GNU time field {name} is missing: {path}")
        result[name] = (
            int(match.group(1))
            if name in {"maximum_rss_kb", "exit_status"}
            else float(match.group(1))
        )
    elapsed = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([0-9:.]+)", text
    )
    if elapsed is None:
        raise ValueError(f"GNU time elapsed field is missing: {path}")
    parts = [float(value) for value in elapsed.group(1).split(":")]
    result["wall_seconds"] = sum(
        value * 60**index for index, value in enumerate(reversed(parts))
    )
    if require_success and (
        result["exit_status"] != 0 or result["maximum_rss_kb"] <= 0
    ):
        raise ValueError(f"GNU time does not describe a successful stage: {path}")
    return result


def docker_metadata(docker: str, image: str) -> dict[str, object]:
    image_id = command_output([docker, "image", "inspect", image, "--format", "{{.Id}}"])
    if re.fullmatch(r"sha256:[0-9a-f]{64}", image_id) is None:
        raise ValueError(f"Docker image lacks a stable sha256 ID: {image_id}")

    def run_entrypoint(entrypoint: str, *arguments: str) -> str:
        return command_output(
            [
                docker,
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                "--entrypoint",
                entrypoint,
                image,
                *arguments,
            ]
        )

    versions = {
        "tidecluster": run_entrypoint("TideCluster.py", "--version").splitlines()[-1],
        "tidehunter": run_entrypoint("TideHunter", "-v").splitlines()[-1],
        "kitehor": run_entrypoint("kitehor", "--version").splitlines()[-1],
        "mmseqs": run_entrypoint("mmseqs", "version").splitlines()[-1],
        "blastn": run_entrypoint("blastn", "-version").splitlines()[0].split()[-1],
    }
    for name, expected in EXPECTED_VERSIONS.items():
        observed = versions[name]
        if expected not in observed:
            raise ValueError(f"Unexpected {name} version: {observed}; expected {expected}")
    return {
        "image": image,
        "image_id": image_id,
        "platform": "linux/amd64",
        "versions": versions,
    }


def run(
    sampling_receipt: Path,
    sample_id: str,
    outdir: Path,
    image: str,
    docker: str,
    cpus: int,
    minimum_length: int,
    interval: float,
) -> dict[str, object]:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    if not 1 <= cpus <= 64 or minimum_length <= 0:
        raise ValueError("cpus must be in [1,64] and minimum length must be positive")
    sampling_receipt = sampling_receipt.resolve()
    input_root = sampling_receipt.parent
    sample_path, sample = load_sample(sampling_receipt, sample_id)
    sequence_lengths, observed_sha256 = reference_metadata(sample_path)
    if (
        observed_sha256 != sample["sha256"]
        or sum(sequence_lengths.values()) != sample["total_bases"]
        or len(sequence_lengths) != sample["window_count"]
    ):
        raise ValueError("Sample FASTA content differs from declared sampling metrics")
    metadata = docker_metadata(docker, image)

    outdir = outdir.resolve()
    outdir.mkdir(parents=True)
    sample_container = container_path(sample_path, input_root, "/input")
    output_container = "/output"
    prefix_container = "/output/tc"
    profile_dir = outdir / "profile"
    base_docker = [
        docker,
        "run",
        "--rm",
        "--platform",
        "linux/amd64",
        "-v",
        f"{input_root}:/input:ro",
        "-v",
        f"{outdir}:/output",
        "-w",
        "/output",
        "--entrypoint",
        "/opt/conda/envs/tidecluster/bin/time",
        image,
        "-v",
    ]
    stage_specs = (
        ("tidehunter", ["tidehunter", "-c", str(cpus), "-pr", prefix_container, "-f", sample_container]),
        (
            "clustering",
            [
                "clustering",
                "-c",
                str(cpus),
                "-m",
                str(minimum_length),
                "-pr",
                prefix_container,
                "-f",
                sample_container,
            ],
        ),
    )
    stages = []
    for name, tidecluster_args in stage_specs:
        internal_time = f"{output_container}/{name}.gnu_time.txt"
        stages.append(
            {
                "name": name,
                "command": [
                    *base_docker,
                    "-o",
                    internal_time,
                    "TideCluster.py",
                    *tidecluster_args,
                ],
                "cwd": str(outdir),
                "scratch_dir": str(outdir),
            }
        )
    stage_manifest = outdir / "stage_manifest.json"
    stage_manifest.write_text(json.dumps({"stages": stages}, indent=2) + "\n")
    project_root = Path(__file__).resolve().parents[2]
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
        "project_commit": command_output(["git", "-C", str(project_root), "rev-parse", "HEAD"]),
        "project_status": command_output(["git", "-C", str(project_root), "status", "--short"]),
        "runner_sha256": digest_file(Path(__file__)),
        "normalizer_sha256": digest_file(project_root / "benchmarks/tidecluster/normalize.py"),
        "profiler_sha256": digest_file(project_root / "benchmarks/scripts/profile_stage_resources.py"),
        "cpus": cpus,
        "minimum_length": minimum_length,
        "resource_scope": {
            "host": "docker_client_process_tree_plus_host_output_scratch_high_water",
            "container": "gnu_time_maximum_rss_for_each_tidecluster_stage",
        },
        "warning": "sampled_real_reference_no_independent_accuracy_truth_no_whole_genome_context",
    }
    environment_path = outdir / "environment.json"
    environment_path.write_text(json.dumps(environment, indent=2) + "\n")
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    for path in (
        Path(__file__),
        project_root / "benchmarks/scripts/run_tidecluster_reference.py",
        project_root / "benchmarks/tidecluster/normalize.py",
        project_root / "benchmarks/scripts/profile_stage_resources.py",
        project_root / "benchmarks/scripts/sample_reference_windows.py",
    ):
        shutil.copyfile(path, snapshot / path.name)

    profile_exit = profile(stage_manifest, profile_dir, interval)
    profile_receipt = json.loads((profile_dir / "receipt.json").read_text())
    if profile_exit != 0 or profile_receipt.get("complete") is not True:
        (outdir / "run_receipt.json").write_text(
            json.dumps({"complete": False, "profile": profile_receipt}, indent=2) + "\n"
        )
        raise RuntimeError("TideCluster stage failed; inspect retained profile logs")

    internal_resources = {
        name: parse_gnu_time(outdir / f"{name}.gnu_time.txt")
        for name, _arguments in stage_specs
    }
    summary = normalize_real_outputs(
        outdir / "tc_tidehunter.gff3",
        outdir / "tc_clustering.gff3_1.gff3",
        outdir / "tc_clustering.gff3",
        sequence_lengths,
        outdir / "normalized",
        family_consensus_fasta=outdir / "tc_consensus/consensus_sequences_all.fasta",
    )
    result = {
        "complete": True,
        "sample_id": sample_id,
        "summary": summary,
        "internal_gnu_time": internal_resources,
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
        "tidehunter.gnu_time.txt",
        "clustering.gnu_time.txt",
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
    parser.add_argument("--image", default="tandemx/tidecluster:1.21.2")
    parser.add_argument("--docker", default="docker")
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--minimum-length", type=int, default=100)
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    try:
        run(
            args.sampling_receipt,
            args.sample_id,
            args.outdir,
            args.image,
            args.docker,
            args.cpus,
            args.minimum_length,
            args.interval,
        )
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
