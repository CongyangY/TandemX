"""Archive a frozen cascade evaluation, including any failed gates and runs."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any

import yaml

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.evaluate_cascade_heldout import validate_matrix


ROOT_FILES = (
    "environment.json",
    "run_config.yaml",
    "validation.json",
    "raw_runs.tsv",
    "summary.tsv",
    "run.log",
)
EVALUATION_FILES = ("gate_results.json", "paired_tidehunter.tsv")
FAILURE_FILES = ("command.json", "receipt.json", "stderr.log", "stdout.log")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _copy(source: Path, destination: Path, relative: Path) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    source_digest = digest_file(source)
    if destination.stat().st_size != source.stat().st_size or digest_file(destination) != source_digest:
        raise OSError(f"Evaluation archive copy differs: {source}")
    return {
        "file": relative.as_posix(),
        "source": str(source),
        "sha256": source_digest,
        "bytes": destination.stat().st_size,
    }


def validate_gate_receipt(
    config: dict[str, Any],
    config_path: Path,
    run: Path,
    evaluation: Path,
) -> tuple[dict[str, object], list[dict[str, str]], dict[str, int]]:
    """Validate matrix and ensure the saved gate receipt describes that matrix."""
    missing = [name for name in ROOT_FILES if not (run / name).is_file()]
    missing.extend(name for name in EVALUATION_FILES if not (evaluation / name).is_file())
    if missing:
        raise ValueError(f"Evaluation evidence lacks required files: {', '.join(missing)}")
    raw = _read_tsv(run / "raw_runs.tsv")
    summary = _read_tsv(run / "summary.tsv")
    matrix = validate_matrix(config, config_path, run, raw, summary)
    gates = json.loads((evaluation / "gate_results.json").read_text())
    environment = json.loads((run / "environment.json").read_text())
    if gates.get("benchmark_id") != config.get("benchmark_id"):
        raise ValueError("Gate receipt benchmark ID differs from the frozen config")
    if gates.get("config_sha256") != digest_file(config_path):
        raise ValueError("Gate receipt config hash differs from the frozen config")
    if gates.get("run_source_digest") != environment.get("source_digest"):
        raise ValueError("Gate receipt source digest differs from the run")
    if gates.get("matrix") != matrix:
        raise ValueError("Gate receipt matrix differs from the validated run matrix")
    gate_rows = gates.get("gates")
    if not isinstance(gate_rows, list) or not gate_rows:
        raise ValueError("Gate receipt contains no gates")
    failed_names = [str(gate["name"]) for gate in gate_rows if gate.get("passed") is not True]
    expected_status = "passed" if not failed_names else "failed"
    if gates.get("status") != expected_status or gates.get("failed_gate_names") != failed_names:
        raise ValueError("Gate receipt status or failed-gate list is internally inconsistent")
    return gates, raw, matrix


def archive(config_path: Path, run: Path, evaluation: Path, outdir: Path) -> dict[str, object]:
    """Create a compact, hash-checked archive without discarding failed runs."""
    config_path = config_path.resolve()
    run = run.resolve()
    evaluation = evaluation.resolve()
    outdir = outdir.resolve()
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    config = yaml.safe_load(config_path.read_text())
    gates, raw, matrix = validate_gate_receipt(config, config_path, run, evaluation)
    failed = [row for row in raw if row["status"] != "ok"]
    evaluation_split = str(config.get("evaluation_split", "heldout"))
    evaluation_seeds = [str(seed) for seed in config["seeds"][evaluation_split]]

    outdir.mkdir(parents=True)
    manifest: list[dict[str, object]] = []
    for name in ROOT_FILES:
        relative = Path("run") / name
        manifest.append(_copy(run / name, outdir / relative, relative))
    for name in EVALUATION_FILES:
        relative = Path("evaluation") / name
        manifest.append(_copy(evaluation / name, outdir / relative, relative))
    config_relative = Path("frozen_config.yaml")
    manifest.append(_copy(config_path, outdir / config_relative, config_relative))

    failure_artifacts: list[dict[str, object]] = []
    for row in failed:
        source_dir = (
            run
            / "runs"
            / row["dataset_id"]
            / row["tool"]
            / f"rep{row['repetition']}"
        )
        for name in FAILURE_FILES:
            source = source_dir / name
            if not source.is_file():
                continue
            relative = (
                Path("failed_runs")
                / row["dataset_id"]
                / row["tool"]
                / f"rep{row['repetition']}"
                / name
            )
            manifest.append(_copy(source, outdir / relative, relative))
            failure_artifacts.append(
                {
                    "dataset_id": row["dataset_id"],
                    "tool": row["tool"],
                    "repetition": row["repetition"],
                    "artifact": relative.as_posix(),
                    "source_sha256": digest_file(source),
                    "source_bytes": source.stat().st_size,
                }
            )

    write_table(outdir / "failed_runs.tsv", failed, list(raw[0]))
    artifact_fields = [
        "dataset_id",
        "tool",
        "repetition",
        "artifact",
        "source_sha256",
        "source_bytes",
    ]
    write_table(outdir / "failure_artifacts.tsv", failure_artifacts, artifact_fields)
    summary = {
        "benchmark_id": config["benchmark_id"],
        "gate_status": gates["status"],
        "failed_gate_names": gates["failed_gate_names"],
        "matrix": matrix,
        "failed_run_count": len(failed),
        "failed_run_tools": sorted({row["tool"] for row in failed}),
        "evaluation_split": evaluation_split,
        "evaluation_seeds": evaluation_seeds,
        "warning": (
            f"{evaluation_split}_seeds_consumed_once;"
            "failed_runs_and_failed_gates_retained;do_not_rerun_or_relabel"
        ),
    }
    (outdir / "archive_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (outdir / "README.md").write_text(
        f"# Cascade {evaluation_split} evidence\n\n"
        "This compact archive preserves the complete frozen matrix, gate receipt, "
        "and receipts/logs for every failed execution. A failed process remains a "
        "process failure and is never converted to zero accuracy. Evaluation seed"
        f"{'s' if len(evaluation_seeds) != 1 else ''} "
        f"{', '.join(evaluation_seeds)} "
        f"{'were' if len(evaluation_seeds) != 1 else 'was'} consumed once and must "
        "not be rerun for model selection.\n"
    )
    for name in ("failed_runs.tsv", "failure_artifacts.tsv", "archive_summary.json", "README.md"):
        path = outdir / name
        manifest.append(
            {
                "file": name,
                "source": "derived_after_validation",
                "sha256": digest_file(path),
                "bytes": path.stat().st_size,
            }
        )
    (outdir / "archive_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--evaluation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        archive(args.config, args.run, args.evaluation, args.outdir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
