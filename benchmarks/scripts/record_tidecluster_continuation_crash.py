#!/usr/bin/env python3
"""Record the observed v2 TideCluster continuation orchestration failure."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.run_tidecluster_docker_reference import parse_gnu_time


STAGE_PATTERN = re.compile(
    r"s(?P<seed>[0-9]+)_(?P<setting>[A-Za-z0-9_]+)_(?P<component>tidehunter|clustering)"
)


def record_failure(source: Path, stage: str) -> dict[str, Any]:
    receipt = source / "run_receipt.json"
    if receipt.exists():
        raise FileExistsError(f"refusing to overwrite failure receipt: {receipt}")
    match = STAGE_PATTERN.fullmatch(stage)
    if match is None or match.group("component") != "tidehunter":
        raise ValueError("expected one TideHunter continuation stage name")
    environment_path = source / "environment.json"
    environment = json.loads(environment_path.read_text(encoding="utf-8"))
    if environment.get("complete") is not False:
        raise ValueError("continuation environment does not retain an incomplete run")
    seed = int(match.group("seed"))
    setting = match.group("setting")
    run_dir = source / f"seed{seed}" / setting
    time_path = run_dir / "tidehunter.gnu_time.txt"
    internal = parse_gnu_time(time_path, require_success=False)
    if internal.get("exit_status") != 0:
        raise ValueError("the externally launched TideHunter stage did not complete")
    outputs = {}
    for name in ("tc_cmd_args.json", "tc_chunks.bed", "tc_tidehunter.gff3"):
        path = run_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"observed external-stage output is missing: {path}")
        outputs[path.relative_to(source).as_posix()] = {
            "bytes": path.stat().st_size,
            "sha256": digest_file(path),
        }
    profile = source / "profile"
    log_paths = {
        "stdout": profile / f"{stage}.stdout.log",
        "stderr": profile / f"{stage}.stderr.log",
    }
    for path in log_paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"continuation log is missing: {path}")
    payload = {
        "schema_version": 1,
        "experiment_id": environment["experiment_id"],
        "complete": False,
        "fate": "orchestration_failure_after_external_stage_start",
        "orchestration_error": {
            "type": "AttributeError",
            "message": "'str' object has no attribute 'exists'",
            "location": "profile_stage_resources.directory_size(stage_scratch_dir)",
        },
        "affected_stage": stage,
        "profile_stage_row_written": False,
        "external_stage_observed_status": "completed_after_profiler_exception",
        "internal_gnu_time": internal,
        "retained_outputs": outputs,
        "logs": {
            name: {
                "file": path.relative_to(source).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": digest_file(path),
            }
            for name, path in log_paths.items()
        },
        "warning": (
            "technical_orchestration_failure_not_TideCluster_accuracy_failure;"
            "completed_external_stage_must_not_be_rerun;dependent_clustering_unstarted"
        ),
    }
    receipt.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--stage", required=True)
    args = parser.parse_args()
    record_failure(args.source, args.stage)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
