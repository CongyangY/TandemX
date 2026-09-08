#!/usr/bin/env python3
"""Profile sequential benchmark stages, including descendant processes.

The manifest is JSON with a top-level ``stages`` list.  Every stage contains a
unique ``name`` and a shell-free ``command`` array; ``cwd`` and ``scratch_dir``
are optional.  Outputs are written incrementally so a failed stage still leaves
its logs and resource trace.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from typing import Any


SAMPLE_FIELDS = (
    "stage",
    "elapsed_seconds",
    "rss_mib",
    "process_count",
    "scratch_bytes",
)
STAGE_FIELDS = (
    "stage",
    "exit_code",
    "wall_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "peak_process_tree_rss_mib",
    "peak_process_count",
    "peak_scratch_bytes",
    "sample_count",
    "stdout_log",
    "stderr_log",
)


def load_manifest(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    stages = payload.get("stages")
    if not isinstance(stages, list) or not stages:
        raise ValueError("manifest must contain a nonempty stages list")
    names: set[str] = set()
    normalized: list[dict[str, Any]] = []
    for index, stage in enumerate(stages, 1):
        if not isinstance(stage, dict):
            raise ValueError(f"stage {index} must be an object")
        name = stage.get("name")
        command = stage.get("command")
        if not isinstance(name, str) or not name or any(char in name for char in "/\\"):
            raise ValueError(f"stage {index} has an invalid name")
        if name in names:
            raise ValueError(f"duplicate stage name: {name}")
        if not isinstance(command, list) or not command or not all(
            isinstance(item, str) and item for item in command
        ):
            raise ValueError(f"stage {name} command must be a nonempty string array")
        cwd = Path(stage["cwd"]).resolve() if stage.get("cwd") else None
        scratch = Path(stage["scratch_dir"]).resolve() if stage.get("scratch_dir") else None
        if cwd is not None and not cwd.is_dir():
            raise ValueError(f"stage {name} cwd is not a directory: {cwd}")
        names.add(name)
        normalized.append({"name": name, "command": command, "cwd": cwd, "scratch_dir": scratch})
    return normalized


def directory_size(path: Path | None) -> int:
    if path is None or not path.exists():
        return 0
    total = 0
    for root, _directories, files in os.walk(path):
        for filename in files:
            try:
                total += (Path(root) / filename).stat().st_size
            except FileNotFoundError:
                continue
    return total


def _linux_process_rows() -> list[tuple[int, int, int]]:
    rows: list[tuple[int, int, int]] = []
    proc = Path("/proc")
    if not proc.is_dir():
        return rows
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            status = (entry / "status").read_text(encoding="utf-8")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        parent = rss_kib = None
        for line in status.splitlines():
            if line.startswith("PPid:"):
                parent = int(line.split()[1])
            elif line.startswith("VmRSS:"):
                rss_kib = int(line.split()[1])
        if parent is not None:
            rows.append((int(entry.name), parent, rss_kib or 0))
    return rows


def _ps_process_rows() -> list[tuple[int, int, int]]:
    try:
        completed = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            check=True,
            text=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return []
    rows = []
    for line in completed.stdout.splitlines():
        fields = line.split()
        if len(fields) == 3 and all(field.isdigit() for field in fields):
            rows.append(tuple(map(int, fields)))
    return rows


def process_tree_rss(root_pid: int) -> tuple[int, int]:
    """Return aggregate resident KiB and process count for a live process tree."""
    rows = _linux_process_rows() or _ps_process_rows()
    children: dict[int, list[int]] = {}
    rss: dict[int, int] = {}
    for pid, parent, rss_kib in rows:
        children.setdefault(parent, []).append(pid)
        rss[pid] = rss_kib
    pending = [root_pid]
    members: set[int] = set()
    while pending:
        pid = pending.pop()
        if pid in members:
            continue
        members.add(pid)
        pending.extend(children.get(pid, ()))
    live = [pid for pid in members if pid in rss]
    return sum(rss[pid] for pid in live), len(live)


def write_table(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run_stage(stage: dict[str, Any], outdir: Path, interval: float) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    name = stage["name"]
    stdout_path = outdir / f"{name}.stdout.log"
    stderr_path = outdir / f"{name}.stderr.log"
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    samples: list[dict[str, Any]] = []
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(stage["command"], cwd=stage["cwd"], stdout=stdout, stderr=stderr)
        while True:
            rss_kib, process_count = process_tree_rss(process.pid)
            elapsed = time.monotonic() - started
            samples.append(
                {
                    "stage": name,
                    "elapsed_seconds": f"{elapsed:.6f}",
                    "rss_mib": f"{rss_kib / 1024:.6f}",
                    "process_count": process_count,
                    "scratch_bytes": directory_size(stage["scratch_dir"]),
                }
            )
            exit_code = process.poll()
            if exit_code is not None:
                break
            time.sleep(interval)
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    wall = time.monotonic() - started
    result = {
        "stage": name,
        "exit_code": exit_code,
        "wall_seconds": f"{wall:.6f}",
        "user_cpu_seconds": f"{usage_after.ru_utime - usage_before.ru_utime:.6f}",
        "system_cpu_seconds": f"{usage_after.ru_stime - usage_before.ru_stime:.6f}",
        "peak_process_tree_rss_mib": f"{max(float(row['rss_mib']) for row in samples):.6f}",
        "peak_process_count": max(int(row["process_count"]) for row in samples),
        "peak_scratch_bytes": max(int(row["scratch_bytes"]) for row in samples),
        "sample_count": len(samples),
        "stdout_log": stdout_path.name,
        "stderr_log": stderr_path.name,
    }
    return result, samples


def profile(manifest: Path, outdir: Path, interval: float) -> int:
    if not 0.05 <= interval <= 60:
        raise ValueError("sample interval must be between 0.05 and 60 seconds")
    stages = load_manifest(manifest)
    outdir.mkdir(parents=True, exist_ok=True)
    stage_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    for stage in stages:
        result, samples = run_stage(stage, outdir, interval)
        stage_rows.append(result)
        sample_rows.extend(samples)
        write_table(outdir / "samples.tsv", SAMPLE_FIELDS, sample_rows)
        write_table(outdir / "stages.tsv", STAGE_FIELDS, stage_rows)
        if result["exit_code"] != 0:
            break
    complete = len(stage_rows) == len(stages) and all(row["exit_code"] == 0 for row in stage_rows)
    receipt = {
        "schema_version": 1,
        "complete": complete,
        "requested_stage_count": len(stages),
        "completed_stage_count": len(stage_rows),
        "sample_interval_seconds": interval,
        "rss_scope": "aggregate_live_process_tree",
        "commands_use_shell": False,
    }
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return 0 if complete else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--interval", type=float, default=0.2)
    args = parser.parse_args()
    try:
        return profile(args.manifest, args.outdir, args.interval)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
