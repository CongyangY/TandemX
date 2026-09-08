#!/usr/bin/env python3
"""Export a successful frozen unitFinder image build and its provenance."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

from benchmarks.scripts.record_unitfinder_build_failure import BUILD_REF, command_output


PROVENANCE_FILES = (
    "unitfinder_tracked_files.sha256",
    "conda_explicit.txt",
    "help.txt",
    "trf_version.txt",
    "nucmer_version.txt",
    "clustalo_version.txt",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate_success_log(text: str, image_id: str) -> None:
    markers = (
        "python /opt/unitFinder/bin/unitFinder.py -h",
        f"writing image {image_id}",
        "naming to docker.io/tandemx/unitfinder:e80bff38-v2 done",
    )
    missing = [marker for marker in markers if marker not in text]
    if missing or "ERROR: failed to solve" in text:
        raise ValueError(f"unitFinder success log is inconsistent: {missing}")


def image_metadata(docker: str, image: str) -> dict[str, Any]:
    rows = json.loads(command_output([docker, "image", "inspect", image]))
    if len(rows) != 1:
        raise ValueError("expected exactly one unitFinder image")
    row = rows[0]
    return {
        "tag": image,
        "image_id": row["Id"],
        "architecture": row["Architecture"],
        "os": row["Os"],
        "size_bytes": row["Size"],
        "repo_digests": row.get("RepoDigests") or [],
        "labels": row.get("Config", {}).get("Labels") or {},
    }


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
        raise ValueError("successful-build Dockerfile differs from frozen config")
    metadata = image_metadata(docker, config["image_tag"])
    if metadata["architecture"] != "amd64" or metadata["os"] != "linux":
        raise ValueError("unitFinder image platform differs from frozen config")
    if metadata["labels"].get("org.opencontainers.image.revision") != config[
        "source_commit"
    ]:
        raise ValueError("unitFinder image source revision differs from frozen config")
    log_completed = subprocess.run(
        [docker, "buildx", "history", "logs", build_ref],
        check=True,
        capture_output=True,
    )
    log_bytes = log_completed.stdout or log_completed.stderr
    log_text = log_bytes.decode("utf-8")
    validate_success_log(log_text, metadata["image_id"])

    outdir.mkdir(parents=True)
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    shutil.copyfile(config_path, snapshot / config_path.name)
    shutil.copyfile(dockerfile, snapshot / dockerfile.name)
    shutil.copyfile(Path(__file__), snapshot / Path(__file__).name)
    log_path = outdir / "docker_build.log"
    log_path.write_bytes(log_bytes)
    provenance_dir = outdir / "image_provenance"
    provenance_dir.mkdir()
    provenance: dict[str, Any] = {}
    for name in PROVENANCE_FILES:
        completed = subprocess.run(
            [
                docker,
                "run",
                "--rm",
                "--platform",
                config["docker"]["platform"],
                "--entrypoint",
                "cat",
                config["image_tag"],
                f"/opt/provenance/{name}",
            ],
            check=True,
            capture_output=True,
        )
        path = provenance_dir / name
        path.write_bytes(completed.stdout)
        if path.stat().st_size == 0:
            raise ValueError(f"empty unitFinder image provenance: {name}")
        provenance[name] = {
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }
    tracked_lines = (provenance_dir / "unitfinder_tracked_files.sha256").read_text(
        encoding="utf-8"
    ).splitlines()
    if not tracked_lines or any(
        re.fullmatch(r"[0-9a-f]{64}  .+", line) is None for line in tracked_lines
    ):
        raise ValueError("invalid unitFinder tracked-source hash manifest")
    project_root = Path(__file__).resolve().parents[2]
    receipt = {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "complete": True,
        "fate": "container_build_passed",
        "build_ref": build_ref,
        "project_commit": command_output(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"]
        ),
        "config_sha256": digest(config_path),
        "dockerfile_sha256": digest(dockerfile),
        "image": metadata,
        "tracked_source_file_count": len(tracked_lines),
        "docker_build_log": {
            "file": log_path.name,
            "bytes": log_path.stat().st_size,
            "sha256": digest(log_path),
        },
        "image_provenance": provenance,
        "smoke_execution_started": False,
        "accuracy_available": False,
        "warning": "container_build_success_only;unitFinder_smoke_and_accuracy_unrun",
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
