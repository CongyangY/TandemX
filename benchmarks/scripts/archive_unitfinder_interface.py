#!/usr/bin/env python3
"""Archive compact unitFinder build or interface-smoke evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

from benchmarks.scripts.run_tidecluster_docker_reference import parse_gnu_time


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate(source: Path, config_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt_path = source / "run_receipt.json"
    if not receipt_path.is_file():
        raise FileNotFoundError(f"unitFinder receipt is missing: {receipt_path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("accuracy_available") not in (None, False):
        raise ValueError("interface evidence cannot claim unitFinder accuracy")
    if receipt.get("complete") is False:
        fate = receipt.get("fate")
        if fate not in {
            "container_build_failure",
            "external_process_failure",
            "interface_contract_failure",
        }:
            raise ValueError("unexpected incomplete unitFinder fate")
        if fate == "container_build_failure":
            log = source / receipt["docker_build_log"]["file"]
            record = receipt["docker_build_log"]
            if (
                not log.is_file()
                or log.stat().st_size != record["bytes"]
                or digest(log) != record["sha256"]
            ):
                raise ValueError("unitFinder failed-build log changed")
            if receipt.get("smoke_execution_started") is not False:
                raise ValueError("failed-build receipt changed smoke execution state")
        else:
            if receipt.get("smoke_execution_started") is not True:
                raise ValueError("failed smoke receipt changed execution state")
            environment = json.loads(
                (source / "environment.json").read_text(encoding="utf-8")
            )
            profile = json.loads(
                (source / "profile/receipt.json").read_text(encoding="utf-8")
            )
            if environment.get("complete") is not False:
                raise ValueError("failed smoke environment changed completion state")
            if fate == "external_process_failure" and (
                profile.get("complete") is not False or receipt.get("profile") != profile
            ):
                raise ValueError("failed smoke profile changed")
            if fate == "interface_contract_failure" and profile.get("complete") is not True:
                raise ValueError("interface-contract failure lacks complete stages")
    elif receipt.get("complete") is True:
        if receipt.get("fate") == "container_build_passed":
            if receipt.get("smoke_execution_started") is not False:
                raise ValueError("unitFinder build-only receipt changed smoke state")
            log = source / receipt["docker_build_log"]["file"]
            record = receipt["docker_build_log"]
            if (
                not log.is_file()
                or log.stat().st_size != record["bytes"]
                or digest(log) != record["sha256"]
            ):
                raise ValueError("unitFinder successful-build log changed")
            for name, record in receipt.get("image_provenance", {}).items():
                path = source / "image_provenance" / name
                if (
                    not path.is_file()
                    or path.stat().st_size != record["bytes"]
                    or digest(path) != record["sha256"]
                ):
                    raise ValueError(f"unitFinder image provenance changed: {name}")
        elif receipt.get("fate") == "smoke_passed":
            if receipt.get("missing_required_outputs") != []:
                raise ValueError("unitFinder success receipt is incomplete")
            for name, record in receipt.get("outputs", {}).items():
                path = source / "run" / name
                if (
                    not path.is_file()
                    or path.stat().st_size != record["bytes"]
                    or digest(path) != record["sha256"]
                ):
                    raise ValueError(f"unitFinder smoke output changed: {name}")
            environment = json.loads(
                (source / "environment.json").read_text(encoding="utf-8")
            )
            profile = json.loads(
                (source / "profile/receipt.json").read_text(encoding="utf-8")
            )
            if environment.get("complete") is not True or profile.get(
                "complete"
            ) is not True:
                raise ValueError("unitFinder success environment/profile is incomplete")
        else:
            raise ValueError("unexpected complete unitFinder fate")
    else:
        raise ValueError("unitFinder receipt lacks an explicit completion state")
    return config, receipt


def archive(source: Path, config_path: Path, outdir: Path) -> dict[str, Any]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite archive: {outdir}")
    source = source.resolve()
    config, receipt = validate(source, config_path)
    outdir.mkdir(parents=True)
    copied: list[dict[str, Any]] = []

    def copy(source_path: Path, target: Path, source_label: str) -> None:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, target)
        if target.stat().st_size != source_path.stat().st_size or digest(target) != digest(
            source_path
        ):
            raise OSError(f"unitFinder archive copy differs: {source_path}")
        copied.append(
            {
                "file": target.relative_to(outdir).as_posix(),
                "source": source_label,
                "bytes": target.stat().st_size,
                "sha256": digest(target),
            }
        )

    copy(config_path, outdir / "preregistration.json", str(config_path.resolve()))
    for path in sorted(
        (path for path in source.rglob("*") if path.is_file()),
        key=lambda path: path.relative_to(source).as_posix(),
    ):
        copy(path, outdir / "results" / path.relative_to(source), str(path))

    resources = None
    time_path = source / "run/unitfinder.gnu_time.txt"
    if time_path.is_file():
        resources = parse_gnu_time(time_path, require_success=receipt["complete"] is True)
    headline = {
        "schema_version": 1,
        "archive_complete": True,
        "experiment_complete": receipt["complete"],
        "experiment_fate": receipt["fate"],
        "evidence_scope": (
            "container_build_only"
            if receipt["fate"] == "container_build_passed"
            else (
                "interface_smoke_only"
                if receipt["fate"] == "smoke_passed"
                else "retained_failure_only"
            )
        ),
        "accuracy_available": False,
        "source_commit": config["source_commit"],
        "resources": resources,
        "boundary": config["boundary"],
    }
    headline_path = outdir / "headline_summary.json"
    headline_path.write_text(
        json.dumps(headline, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    readme = outdir / "README.md"
    if receipt["fate"] == "container_build_passed":
        description = (
            "The corrected pinned unitFinder image built successfully, and its "
            "BuildKit log, image ID, dependency lock and upstream tracked-file "
            "hashes are retained. The interface smoke had not started, so this is "
            "container-build evidence only and not an accuracy result.\n"
        )
    elif receipt["complete"]:
        description = (
            "The pinned unitFinder source completed its frozen single-chromosome "
            "installation/interface smoke. This establishes executable packaging "
            "and required native outputs only; it is not an accuracy result.\n"
        )
    elif receipt["fate"] == "container_build_failure":
        description = (
            "The first frozen unitFinder container build failed before the upstream "
            "interface started. The complete BuildKit log and source snapshot are "
            "retained; this is a container-definition failure, not an algorithm-"
            "accuracy measurement.\n"
        )
    else:
        description = (
            "The frozen unitFinder interface smoke started but did not satisfy its "
            "predeclared execution/output contract. All partial outputs, logs and "
            "resource records are retained. This is an interface execution fate, "
            "not a zero-accuracy measurement.\n"
        )
    readme.write_text(f"# {config['experiment_id']}\n\n{description}", encoding="utf-8")
    for path in (headline_path, readme):
        copied.append(
            {
                "file": path.name,
                "source": "generated_by_archive_unitfinder_interface",
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        )
    manifest = {
        "schema_version": 1,
        "archive_complete": True,
        "experiment_complete": receipt["complete"],
        "file_count": len(copied),
        "files": copied,
    }
    (outdir / "archive_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    archive(args.source, args.config, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
