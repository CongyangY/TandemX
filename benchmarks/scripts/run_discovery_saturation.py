"""Run the pre-specified TandemX discovery-saturation validation."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import shutil
import sys

from benchmarks.challenge.adapters import build_command
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, write_table


def validate_dataset(dataset: Path, config_path: Path, seed: int) -> list[dict]:
    receipt = json.loads((dataset / "generation_receipt.json").read_text())
    if (
        not receipt.get("complete")
        or receipt.get("split") != "validation"
        or not receipt.get("validation_used")
        or receipt.get("heldout_used")
        or int(receipt.get("seed", -1)) != seed
        or receipt.get("config_sha256") != digest_file(config_path)
    ):
        raise ValueError(f"Invalid validation dataset: {dataset}")
    copied = dataset / "run_config.json"
    if digest_file(copied) != digest_file(config_path):
        raise ValueError(f"Copied config differs: {dataset}")
    conditions = receipt.get("conditions_completed", [])
    if len(conditions) != 7 or any(row.get("label") != "hifi_like_iid" for row in conditions):
        raise ValueError(f"Unexpected condition set: {dataset}")
    for row in conditions:
        folder = dataset / "reads" / row["condition_id"]
        if digest_file(folder / "manifest.json") != row["manifest_sha256"]:
            raise ValueError(f"Condition manifest differs: {folder}")
    return conditions


def run_one(job: dict, snapshot: Path, timeout: float) -> dict:
    folder = job["run_path"]
    folder.mkdir(parents=True, exist_ok=False)
    command, _ = build_command(
        "tandemx",
        sys.executable,
        job["reads"],
        folder,
        30,
        1000,
        100,
        threads=1,
    )
    command[:1] = [sys.executable, "-m", "tandemx.cli"]
    command += [
        "--discovery-method",
        "elastic",
        "--clustering-method",
        "sequence",
        "--cluster-identity",
        ".95",
        "--family-audit",
        "related",
    ]
    (folder / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    measured = run_process(
        command,
        folder / "stdout.log",
        folder / "stderr.log",
        timeout,
        {**os.environ, "PYTHONPATH": str(snapshot)},
        snapshot,
    )
    status = "ok" if measured["exit_code"] == 0 and not measured["timed_out"] else "failed"
    row = {
        "seed": job["seed"],
        "condition_id": job["condition_id"],
        "coverage": job["coverage"],
        "read_count": job["read_count"],
        "total_bases": job["total_bases"],
        "run_path": str(folder),
        "status": status,
        **measured,
    }
    (folder / "execution.json").write_text(json.dumps(row, indent=2) + "\n")
    return row


def run(
    config_path: Path,
    datasets: list[Path],
    outdir: Path,
    timeout: float,
    workers: int,
) -> None:
    config = json.loads(config_path.read_text())
    seeds = [int(value) for value in config["seeds"]["validation"]]
    coverages = [float(value) for value in config["coverages"]]
    if (
        len(datasets) != len(seeds)
        or len({path.resolve() for path in datasets}) != len(datasets)
        or coverages != [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
        or not math.isfinite(timeout)
        or timeout <= 0
        or not 1 <= workers <= 3
    ):
        raise ValueError("Dataset, depth, timeout or worker contract differs from protocol")
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    snapshot = outdir / "source_snapshot"
    provenance = source_manifest(root, snapshot)
    extra_paths = [
        Path(__file__),
        root / "benchmarks/discovery/saturation.py",
        root / "benchmarks/configs/discovery_saturation_validation_v1.json",
        root / "docs/discovery_saturation_protocol_20260910.md",
    ]
    extra_hashes = {}
    for path in extra_paths:
        relative = path.relative_to(root)
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        extra_hashes[str(relative)] = digest_file(target)
    jobs = []
    dataset_rows = []
    for seed, dataset in zip(seeds, datasets, strict=True):
        dataset = dataset.resolve()
        conditions = validate_dataset(dataset, config_path, seed)
        observed = [float(row["coverage"]) for row in conditions]
        if observed != coverages:
            raise ValueError(f"Coverage order differs for seed {seed}: {observed}")
        dataset_rows.append(
            {
                "seed": seed,
                "path": str(dataset),
                "generation_receipt_sha256": digest_file(dataset / "generation_receipt.json"),
            }
        )
        for row in conditions:
            jobs.append(
                {
                    "seed": seed,
                    "condition_id": row["condition_id"],
                    "coverage": float(row["coverage"]),
                    "read_count": int(row["read_count"]),
                    "total_bases": int(row["total_bases"]),
                    "reads": dataset / "reads" / row["condition_id"] / "reads.fa",
                    "run_path": outdir
                    / "runs"
                    / f"s{seed}"
                    / f"depth_{float(row['coverage']):g}x",
                }
            )
    provenance.update(
        config=str(config_path.resolve()),
        config_sha256=digest_file(config_path),
        protocol=str((root / "docs/discovery_saturation_protocol_20260910.md").resolve()),
        datasets=dataset_rows,
        extra_source_sha256=extra_hashes,
        workers=workers,
        threads_per_run=1,
        timing_use="descriptive_only_concurrent_across_seeds",
    )
    (outdir / "environment.json").write_text(json.dumps(provenance, indent=2) + "\n")
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_one, job, snapshot, timeout): job for job in jobs}
        for future in as_completed(futures):
            rows.append(future.result())
            rows.sort(key=lambda row: (row["seed"], row["coverage"]))
            write_table(outdir / "run_summary.tsv", rows, list(rows[0]))
    complete = len(rows) == 21 and all(row["status"] == "ok" for row in rows)
    receipt = {
        "complete": complete,
        "expected_runs": 21,
        "completed_runs": sum(row["status"] == "ok" for row in rows),
        "failed_runs": sum(row["status"] != "ok" for row in rows),
    }
    (outdir / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if not complete:
        raise RuntimeError("Discovery saturation run incomplete; outputs retained")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--datasets", type=Path, nargs=3, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    run(args.config, args.datasets, args.outdir, args.timeout, args.workers)
