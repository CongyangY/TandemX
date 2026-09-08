#!/usr/bin/env python3
"""Export and classify a failed frozen unitFinder Docker build."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any


BUILD_REF = re.compile(r"[a-z0-9]+")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate_failure_log(text: str, source_commit: str) -> dict[str, Any]:
    required = {
        "source_checkout_reached": f'HEAD is now at {source_commit[:7]}',
        "relative_hash_path_failure": "sha256sum: LICENSE: No such file or directory",
        "build_exit_123": "exit code: 123",
    }
    missing = [name for name, marker in required.items() if marker not in text]
    if missing:
        raise ValueError(f"unitFinder build log lacks expected failure evidence: {missing}")
    missing_paths = sorted(
        set(re.findall(r"sha256sum: ([^:\n]+): No such file or directory", text))
    )
    if not missing_paths:
        raise ValueError("unitFinder build log contains no failed tracked paths")
    return {
        "classification": "container_provenance_hash_working_directory_error",
        "missing_relative_path_count": len(missing_paths),
        "first_missing_relative_path": missing_paths[0],
        "source_checkout_reached": True,
        "dependency_install_reached": "Transaction finished" in text,
        "upstream_interface_started": False,
    }


def command_output(command: list[str]) -> str:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return (completed.stdout or completed.stderr).strip()


def record(
    config_path: Path,
    dockerfile: Path,
    build_ref: str,
    outdir: Path,
    docker: str = "docker",
) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"choose a new output directory: {outdir}")
    if BUILD_REF.fullmatch(build_ref) is None:
        raise ValueError("unsafe Docker BuildKit build reference")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if digest(dockerfile) != config["dockerfile_sha256"]:
        raise ValueError("failed-build Dockerfile differs from frozen config")
    completed = subprocess.run(
        [docker, "buildx", "history", "logs", build_ref],
        check=True,
        capture_output=True,
        text=True,
    )
    log_text = completed.stdout or completed.stderr
    classification = validate_failure_log(log_text, config["source_commit"])
    image = subprocess.run(
        [docker, "image", "inspect", config["image_tag"]],
        check=False,
        capture_output=True,
        text=True,
    )
    if image.returncode == 0:
        raise ValueError("failed unitFinder build unexpectedly produced the target image")

    outdir.mkdir(parents=True)
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    shutil.copyfile(config_path, snapshot / config_path.name)
    shutil.copyfile(dockerfile, snapshot / dockerfile.name)
    shutil.copyfile(Path(__file__), snapshot / Path(__file__).name)
    log_path = outdir / "docker_build.log"
    log_path.write_text(log_text, encoding="utf-8")
    project_root = Path(__file__).resolve().parents[2]
    receipt = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "complete": False,
        "fate": "container_build_failure",
        "build_ref": build_ref,
        "attempted_command": [
            docker,
            "build",
            "--platform",
            config["docker"]["platform"],
            "-t",
            config["image_tag"],
            str(dockerfile.parent),
        ],
        "project_commit": command_output(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"]
        ),
        "config_sha256": digest(config_path),
        "dockerfile_sha256": digest(dockerfile),
        "docker_build_log": {
            "file": log_path.name,
            "bytes": log_path.stat().st_size,
            "sha256": digest(log_path),
        },
        "classification": classification,
        "target_image_present_after_failure": False,
        "smoke_execution_started": False,
        "accuracy_available": False,
        "warning": (
            "retained_container_definition_failure_not_unitFinder_algorithm_failure;"
            "no_accuracy_or_runtime_claim"
        ),
    }
    (outdir / "run_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dockerfile", required=True, type=Path)
    parser.add_argument("--build-ref", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--docker", default="docker")
    args = parser.parse_args()
    record(args.config, args.dockerfile, args.build_ref, args.outdir, args.docker)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
