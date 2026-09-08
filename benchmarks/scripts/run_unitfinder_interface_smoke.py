#!/usr/bin/env python3
"""Run the frozen unitFinder single-chromosome installation/interface smoke."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any

from benchmarks.scripts.profile_stage_resources import profile


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def docker_image_metadata(docker: str, image: str) -> dict[str, Any]:
    completed = subprocess.run(
        [docker, "image", "inspect", image],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = json.loads(completed.stdout)
    if len(rows) != 1:
        raise ValueError(f"expected one Docker image: {image}")
    row = rows[0]
    return {
        "tag": image,
        "image_id": row["Id"],
        "architecture": row["Architecture"],
        "os": row["Os"],
        "repo_digests": row.get("RepoDigests") or [],
        "labels": row.get("Config", {}).get("Labels") or {},
    }


def validate_inputs(config: dict[str, Any], config_path: Path) -> tuple[Path, Path]:
    root = Path(__file__).resolve().parents[2]
    dockerfile = root / "benchmarks/containers/unitfinder/Dockerfile"
    if digest(dockerfile) != config["dockerfile_sha256"]:
        raise ValueError("unitFinder Dockerfile differs from the frozen config")
    input_dir = Path(config["input"]["directory"]).resolve()
    fasta = input_dir / config["input"]["fasta"]
    truth = input_dir / config["input"]["truth"]
    manifest = input_dir / "manifest.json"
    checks = (
        (fasta, config["input"]["fasta_sha256"]),
        (truth, config["input"]["truth_sha256"]),
        (manifest, config["input"]["manifest_sha256"]),
    )
    for path, expected in checks:
        if not path.is_file() or digest(path) != expected:
            raise ValueError(f"unitFinder frozen smoke input differs: {path}")
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    if manifest_payload.get("sequence_length_bp") != config["input"]["sequence_length_bp"]:
        raise ValueError("unitFinder smoke sequence length differs from frozen config")
    if digest(config_path) == "":
        raise AssertionError("unreachable empty config digest")
    return input_dir, dockerfile


def build_stage_manifest(
    config: dict[str, Any], input_dir: Path, run_dir: Path, docker: str
) -> dict[str, Any]:
    image = config["image_tag"]
    base = [
        docker,
        "run",
        "--rm",
        "--platform",
        config["docker"]["platform"],
        "--network",
        config["docker"]["network"],
        "-v",
        f"{input_dir}:/input:ro",
        "-v",
        f"{run_dir}:/output",
        "-w",
        "/output",
    ]
    return {
        "stages": [
            {
                "name": "help",
                "command": [*base, image, "-h"],
                "cwd": str(run_dir),
                "scratch_dir": str(run_dir),
            },
            {
                "name": "smoke",
                "command": [
                    *base,
                    "--entrypoint",
                    "/opt/conda/envs/unitfinder/bin/time",
                    image,
                    "-v",
                    "-o",
                    "/output/unitfinder.gnu_time.txt",
                    "python",
                    "/opt/unitFinder/bin/unitFinder.py",
                    *config["command_arguments"],
                ],
                "cwd": str(run_dir),
                "scratch_dir": str(run_dir),
            },
        ]
    }


def run(
    config_path: Path,
    outdir: Path,
    docker: str = "docker",
    interval: float = 0.2,
) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"choose a new output directory: {outdir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    input_dir, dockerfile = validate_inputs(config, config_path)
    metadata = docker_image_metadata(docker, config["image_tag"])
    labels = metadata["labels"]
    if labels.get("org.opencontainers.image.revision") != config["source_commit"]:
        raise ValueError("unitFinder image source revision differs from frozen config")
    if metadata["architecture"] != "amd64" or metadata["os"] != "linux":
        raise ValueError("unitFinder image platform differs from frozen config")

    outdir.mkdir(parents=True)
    run_dir = outdir / "run"
    run_dir.mkdir()
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    for path in (config_path, dockerfile, Path(__file__)):
        shutil.copyfile(path, snapshot / path.name)
    stage_manifest = outdir / "stage_manifest.json"
    stage_manifest.write_text(
        json.dumps(build_stage_manifest(config, input_dir, run_dir, docker), indent=2)
        + "\n",
        encoding="utf-8",
    )
    environment = {
        "schema_version": 1,
        "complete": False,
        "experiment_id": config["experiment_id"],
        "config_sha256": digest(config_path),
        "dockerfile_sha256": digest(dockerfile),
        "image": metadata,
        "input_manifest_sha256": digest(input_dir / "manifest.json"),
        "boundary": config["boundary"],
    }
    environment_path = outdir / "environment.json"
    environment_path.write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    exit_code = profile(stage_manifest, outdir / "profile", interval)
    profile_receipt = json.loads(
        (outdir / "profile/receipt.json").read_text(encoding="utf-8")
    )
    if exit_code != 0 or profile_receipt.get("complete") is not True:
        receipt = {
            "schema_version": 1,
            "complete": False,
            "fate": "external_process_failure",
            "profile": profile_receipt,
            "warning": "retained_failure_is_not_a_zero_accuracy_measurement",
        }
        (outdir / "run_receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return receipt

    required = [run_dir / name for name in config["acceptance"]["required_nonempty_outputs"]]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    output_files = sorted(path for path in run_dir.rglob("*") if path.is_file())
    if missing:
        fate, complete = "interface_contract_failure", False
    else:
        fate, complete = "smoke_passed", True
    receipt = {
        "schema_version": 1,
        "complete": complete,
        "fate": fate,
        "missing_required_outputs": missing,
        "output_file_count": len(output_files),
        "outputs": {
            path.relative_to(run_dir).as_posix(): {
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
            for path in output_files
        },
        "boundary": config["boundary"],
    }
    (outdir / "run_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    environment["complete"] = complete
    environment_path.write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--docker", default="docker")
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    receipt = run(args.config, args.outdir, args.docker, args.interval)
    return 0 if receipt["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
