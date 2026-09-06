"""Run controlled per-array comparisons with strict failure reporting.

Example: python -m benchmarks.challenge.run --config
benchmarks/configs/challenge_v1.yaml --outdir /tmp/tandemx-challenge
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import importlib.metadata
import os
import platform
import random
import shutil
import signal
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import yaml

from .adapters import build_command, parse_arrays, read_fasta
from .evaluate import score_arrays, score_families
from .sequence_metrics import score_cyclic_recovery
from .schema import ArrayRecord, digest_file, read_table, write_table
from .simulate import Scenario, generate_dataset

LOGGER = logging.getLogger(__name__)


def json_safe(value):
    """Standard JSON uses null for unavailable metrics, never NaN tokens."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return None if isinstance(value, float) and not math.isfinite(value) else value


def run_process(command: list[str], stdout: Path, stderr: Path, timeout: float,
                environment: dict[str, str] | None = None, working_directory: Path | None = None) -> dict:
    if timeout <= 0:
        raise ValueError("Timeout must be positive")
    start = time.perf_counter()
    timed_out = False
    with stdout.open("w") as out, stderr.open("w") as err:
        process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True,
                                   env=environment, cwd=working_directory)
        if hasattr(os, "wait4"):
            while True:
                pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                if pid:
                    process.returncode = os.waitstatus_to_exitcode(status)
                    break
                if time.perf_counter() - start > timeout:
                    timed_out = True
                    os.killpg(process.pid, signal.SIGKILL)
                    _, status, usage = os.wait4(process.pid, 0)
                    process.returncode = os.waitstatus_to_exitcode(status)
                    break
                time.sleep(0.02)
            peak = usage.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024)
            user_cpu, system_cpu = usage.ru_utime, usage.ru_stime
        else:  # pragma: no cover - publication benchmarks require a Unix host
            try:
                process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                process.kill()
                process.wait()
            peak = math.nan
            user_cpu = system_cpu = math.nan
    return {"exit_code": process.returncode, "runtime_seconds": time.perf_counter() - start,
            "peak_rss_mib": peak, "cpu_user_seconds": user_cpu, "cpu_system_seconds": system_cpu,
            "timed_out": timed_out}


def source_manifest(root: Path, snapshot: Path | None = None) -> dict:
    paths = [p for folder in (root / "tandemx", root / "benchmarks" / "challenge", root / "rust-core" / "src")
             for p in folder.rglob("*") if p.is_file() and p.suffix in {".py", ".rs", ".so", ".pyd"}]
    paths.extend(p for p in (root / "pyproject.toml", root / "environment.yml", root / "LICENSE", root / "README.md",
                            root / "rust-core" / "Cargo.toml", root / "rust-core" / "Cargo.lock") if p.is_file())
    hashes = {str(p.relative_to(root)): digest_file(p) for p in sorted(paths)}
    if snapshot is not None:
        snapshot.mkdir(parents=True, exist_ok=False)
        for path in paths:
            target = snapshot / path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            if digest_file(target) != hashes[str(path.relative_to(root))]:
                raise ValueError(f"Source changed while snapshotting: {path}")
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=False)
    commit = revision.stdout.strip() if revision.returncode == 0 else None
    return {"git_head": commit, "revision_warning": None if commit else "not_a_git_checkout_use_source_digest", "file_hashes": hashes,
            "source_snapshot": str(snapshot) if snapshot else None,
            "source_digest": hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()}


def run_suite(config_path: Path, outdir: Path, split: str, selected: list[str] | None = None) -> int:
    config = yaml.safe_load(config_path.read_text())
    root = Path(__file__).resolve().parents[2]
    outdir = outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Choose a new empty output directory: {outdir}")
    if split not in config["seeds"]:
        raise ValueError(f"Unknown split: {split}")
    if set(config["tools"]) - {"tandemx", "trf", "tidehunter", "ultra"}:
        raise ValueError("Unsupported comparator in this array-level benchmark")
    ultra_options = config.get("ultra_options", {})
    if set(ultra_options) - {"window_size", "windows", "tune", "tune_indel"}:
        raise ValueError("Unsupported ULTRA option")
    for key in ("window_size", "windows"):
        if key in ultra_options and (type(ultra_options[key]) is not int or ultra_options[key] < 1):
            raise ValueError(f"ULTRA {key} must be a positive integer")
    for key in ("tune", "tune_indel"):
        if key in ultra_options and type(ultra_options[key]) is not bool:
            raise ValueError(f"ULTRA {key} must be boolean")
    executables = {}
    for tool, value in config["tools"].items():
        resolved = shutil.which(value) or str((root / value).resolve())
        if not Path(resolved).is_file() or not os.access(resolved, os.X_OK):
            raise ValueError(f"Missing executable {tool}: {resolved}")
        executables[tool] = resolved
    scenarios = [Scenario(**row) for row in config["scenarios"] if selected is None or row["name"] in selected]
    if not scenarios or (selected and set(selected) - {s.name for s in scenarios}):
        raise ValueError("No scenarios or an unknown scenario was requested")
    for scenario in scenarios:
        scenario.validate()
    if int(config["repetitions"]) < 1 or float(config["timeout_seconds"]) <= 0:
        raise ValueError("Repetitions and timeout must be positive")
    # Independent random families/seeds: no truth is passed through build_command.
    if len({s for values in config["seeds"].values() for s in values}) != sum(map(len, config["seeds"].values())):
        raise ValueError("Seed leakage: development and test seeds must be disjoint")
    outdir.mkdir(parents=True)
    handler = logging.FileHandler(outdir / "run.log")
    LOGGER.addHandler(handler)
    LOGGER.setLevel(logging.INFO)
    (outdir / "run_config.yaml").write_text(yaml.safe_dump({**config, "selected_split": split,
                                                          "selected_scenarios": selected}, sort_keys=False))
    manifest = source_manifest(root, outdir / "source_snapshot")
    versions = {}
    for tool, executable in executables.items():
        version_arg = "--version" if tool == "tandemx" else "-h" if tool == "ultra" else "-v"
        result = subprocess.run([executable, version_arg],
                                capture_output=True, text=True, timeout=10, check=False)
        versions[tool] = {"exit_code": result.returncode, "self_report": (result.stdout + result.stderr)[:4000]}
    manifest.update({"config_sha256": digest_file(config_path), "python": sys.version, "platform": platform.platform(),
                     "independent_evaluator": {"edlib": importlib.metadata.version("edlib")},
                     "executables": {t: {"path": p, "sha256": digest_file(Path(p))} for t, p in executables.items()},
                     "tool_versions": versions, "memory_method": "wait4 direct child ru_maxrss", "threads": 1, "split": split,
                     "family_recovery_source": "TandemX final catalog; other tools per-array consensuses",
                     "tandemx_execution_source": "source_snapshot via PYTHONPATH and working directory"})
    (outdir / "environment.json").write_text(json.dumps(manifest, indent=2) + "\n")
    raw: list[dict] = []
    min_period, max_period = config["period_range"]
    min_span = config["minimum_span_bp"]
    for scenario in scenarios:
        for seed in config["seeds"][split]:
            dataset_id = f"{scenario.name}_s{seed}"
            dataset = outdir / "datasets" / dataset_id
            generated = generate_dataset(scenario, seed, dataset)
            truth = [ArrayRecord(r["read_id"], int(r["start"]), int(r["end"]), int(r["period"]), r["sequence"], r["family_id"])
                     for r in read_table(dataset / "truth_arrays.tsv")]
            lengths = {r["read_id"]: int(r["length_bp"]) for r in read_table(dataset / "truth_reads.tsv")}
            observed_families = {r.family_id: r.sequence for r in truth}
            for repetition in range(1, int(config["repetitions"]) + 1):
                order = sorted(executables)
                random.Random(seed + repetition).shuffle(order)
                for tool in order:
                    run = outdir / "runs" / dataset_id / tool / f"rep{repetition}"
                    run.mkdir(parents=True)
                    command, output = build_command(tool, executables[tool], dataset / "reads.fa", run,
                                                    min_period, max_period, min_span)
                    if tool == "tandemx" and "discovery_method" in config:
                        command.extend(["--discovery-method", str(config["discovery_method"])])
                    if tool == "tandemx":
                        for key in ("clustering_method", "cluster_identity"):
                            if key in config:
                                command.extend(["--" + key.replace("_", "-"), str(config[key])])
                    if tool == "ultra":
                        for key, flag in (("window_size", "--win_size"), ("windows", "--windows")):
                            if key in ultra_options:
                                command.extend([flag, str(ultra_options[key])])
                        for key in ("tune", "tune_indel"):
                            if ultra_options.get(key):
                                command.append("--" + key)
                    (run / "command.json").write_text(json.dumps(command, indent=2) + "\n")
                    row = {"scenario": scenario.name, "dataset_id": dataset_id, "seed": seed, "split": split,
                           "tool": tool, "repetition": repetition, "total_bases": generated["total_bases"],
                           **run_process(command, output if tool == "trf" else run / "stdout.log", run / "stderr.log",
                                         config["timeout_seconds"],
                                         {**os.environ, "PYTHONPATH": str(outdir / "source_snapshot")} if tool == "tandemx" else None,
                                         outdir / "source_snapshot" if tool == "tandemx" else None)}
                    row["status"] = "failed" if row["exit_code"] or row["timed_out"] else "ok"
                    row["error"] = ""
                    if row["status"] == "ok":
                        try:
                            predictions = parse_arrays(tool, output, min_period, max_period, min_span)
                            metrics, matches = score_arrays(predictions, truth, lengths, config["minimum_iou"])
                            sequences = (list(read_fasta(output.parent / "monomers.fa").values()) if tool == "tandemx"
                                         else [r.sequence for r in predictions])
                            families, family_rows = score_families(sequences, observed_families)
                            cyclic_metrics, cyclic_rows = score_cyclic_recovery(sequences, observed_families)
                            row.update(metrics)
                            row.update(families)
                            row.update(cyclic_metrics)
                            norm = [asdict(r) for r in sorted(predictions, key=lambda r: (r.read_id, r.start, r.end, r.period, r.sequence))]
                            row["prediction_sha256"] = hashlib.sha256(json.dumps(norm, sort_keys=True).encode()).hexdigest()
                            row["catalog_sha256"] = hashlib.sha256(json.dumps(sorted(set(sequences))).encode()).hexdigest()
                            write_table(run / "predictions.tsv", norm, list(ArrayRecord.__dataclass_fields__))
                            write_table(run / "matches.tsv", matches, list(matches[0]) if matches else ["prediction_index", "status"])
                            write_table(run / "family_recovery.tsv", family_rows,
                                        list(family_rows[0]) if family_rows else ["family_id", "recovered", "criterion"])
                            write_table(run / "cyclic_monomer_recovery.tsv", cyclic_rows,
                                        list(cyclic_rows[0]) if cyclic_rows else ["truth_id", "recovered", "criterion"])
                        except (OSError, ValueError, KeyError, IndexError) as error:
                            row["status"], row["error"] = "invalid_output", str(error)
                    raw.append(row)
                    fields = list(dict.fromkeys(key for r in raw for key in r))
                    write_table(outdir / "raw_runs.tsv", ({k: r.get(k, "NA") for k in fields} for r in raw), fields)
                    (run / "receipt.json").write_text(json.dumps(json_safe(row), indent=2, allow_nan=False) + "\n")
                    LOGGER.info("%s %s rep%d status=%s time=%.3fs array_recall=%s", dataset_id, tool, repetition,
                                row["status"], row["runtime_seconds"], row.get("array_recall", "NA"))
    summary = summarize(raw)
    write_table(outdir / "summary.tsv", summary, list(summary[0]))
    receipt = {"complete": True, "successful_runs": sum(r["status"] == "ok" for r in raw),
               "failed_runs": sum(r["status"] != "ok" for r in raw), "total_runs": len(raw)}
    (outdir / "validation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    LOGGER.removeHandler(handler)
    handler.close()
    return int(receipt["failed_runs"] > 0)


def summarize(raw: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for row in raw:
        groups.setdefault((row["scenario"], row["dataset_id"], row["tool"]), []).append(row)
    result = []
    metric_names = ["array_recall", "array_precision", "array_f1", "read_detection_recall", "read_detection_precision",
                    "negative_read_call_rate", "sequence_family_recall", "matched_period_mae_bp", "matched_boundary_mae_bp",
                    "cyclic_monomer_recall", "homologous_consensus_fraction", "mean_best_cyclic_edit_similarity",
                    "base_union_recall", "base_union_precision", "base_union_f1", "duplicate_bp_fraction"]
    for (scenario, dataset, tool), rows in sorted(groups.items()):
        good = [r for r in rows if r["status"] == "ok"]
        valid = len(good) == len(rows)
        deterministic = valid and len({(r["prediction_sha256"], r["catalog_sha256"]) for r in good}) == 1
        row = {"scenario": scenario, "dataset_id": dataset, "tool": tool, "seed": rows[0]["seed"],
               "split": rows[0]["split"], "successful_runs": len(good), "attempted_runs": len(rows),
               "deterministic": deterministic if len(rows) > 1 else "not_tested_single_run",
               "median_runtime_seconds": statistics.median(r["runtime_seconds"] for r in good) if valid else "NA",
               "median_peak_rss_mib": statistics.median(r["peak_rss_mib"] for r in good) if valid else "NA"}
        row.update({name: good[0].get(name, "NA") if valid and deterministic else "NA" for name in metric_names})
        result.append(row)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--split", default="development")
    parser.add_argument("--scenarios", nargs="+")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    try:
        return run_suite(args.config, args.outdir, args.split, args.scenarios)
    except (ValueError, OSError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
