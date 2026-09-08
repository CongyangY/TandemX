"""Fail-soft continuation of a frozen TideCluster factorial validation."""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
import json
from pathlib import Path
import shutil
from typing import Any

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.profile_stage_resources import (
    SAMPLE_FIELDS,
    STAGE_FIELDS,
    run_stage,
    write_table as write_profile_table,
)
from benchmarks.scripts.run_tidecluster_docker_reference import (
    docker_metadata,
    parse_gnu_time,
)
from benchmarks.tidecluster.factorial_evaluate import evaluate_factorial
from benchmarks.tidecluster.factorial_run import (
    SUMMARY_FIELDS,
    build_stage_specs,
    command_output,
    validate_datasets,
)


CELL_FATE_FIELDS = (
    "seed",
    "setting",
    "status",
    "tidehunter_status",
    "clustering_status",
    "evaluation_status",
    "source",
    "reason",
)
STAGE_FATE_FIELDS = (
    "stage",
    "seed",
    "setting",
    "component",
    "status",
    "execution_source",
    "exit_code",
    "reason",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_parent_failure(config: dict[str, Any]) -> dict[str, Any]:
    parent = Path(config["parent_failure_dir"]).resolve()
    if not parent.is_dir():
        raise FileNotFoundError(f"parent failure directory is missing: {parent}")
    artifacts = config.get("parent_artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise ValueError("continuation config requires parent artifact hashes")
    for name, expected in artifacts.items():
        path = parent / name
        if not path.is_file() or digest_file(path) != expected:
            raise ValueError(f"parent failure artifact changed: {name}")
    receipt = json.loads((parent / "run_receipt.json").read_text(encoding="utf-8"))
    if receipt.get("complete") is not False or receipt.get("fate") != "external_process_failure":
        raise ValueError("parent does not retain the expected external-process failure")
    stage_rows = read_tsv(parent / "profile/stages.tsv")
    failed = {row["stage"] for row in stage_rows if int(row["exit_code"]) != 0}
    frozen_skip = set(config.get("never_rerun_stages", []))
    if failed != frozen_skip or len(failed) != 1:
        raise ValueError("never-rerun stages differ from the parent failure")
    failed_stage = next(iter(failed))
    if not failed_stage.endswith("_tidehunter"):
        raise ValueError("continuation supports a parent TideHunter-stage failure only")
    return {
        "directory": parent,
        "receipt": receipt,
        "stage_rows": stage_rows,
        "failed_stage": failed_stage,
    }


def blank_summary(
    seed: int,
    setting: str,
    status: str,
    warning: str,
    tidehunter: dict[str, float | int] | None = None,
    clustering: dict[str, float | int] | None = None,
) -> dict[str, Any]:
    row: dict[str, Any] = {field: "" for field in SUMMARY_FIELDS}
    row.update(
        {
            "seed": seed,
            "setting": setting,
            "status": status,
            "warning": warning,
        }
    )
    for name, resource in (("tidehunter", tidehunter), ("clustering", clustering)):
        if resource is not None:
            row[f"{name}_wall_seconds"] = resource["wall_seconds"]
            row[f"{name}_maximum_rss_kb"] = resource["maximum_rss_kb"]
    return row


def retained_gnu_time(path: Path) -> dict[str, float | int] | None:
    try:
        return parse_gnu_time(path, require_success=False)
    except (OSError, ValueError):
        return None


def success_summary(
    seed: int,
    setting: str,
    metrics: dict[str, Any],
    tidehunter: dict[str, float | int],
    clustering: dict[str, float | int],
) -> dict[str, Any]:
    row = blank_summary(seed, setting, "ok", metrics.get("warning", ""), tidehunter, clustering)
    for field in SUMMARY_FIELDS[7:]:
        if field != "warning":
            row[field] = metrics.get(field, "")
    return row


def _stage_fate(
    stage: str,
    seed: int,
    setting: str,
    component: str,
    status: str,
    source: str,
    exit_code: int | str,
    reason: str,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "seed": seed,
        "setting": setting,
        "component": component,
        "status": status,
        "execution_source": source,
        "exit_code": exit_code,
        "reason": reason,
    }


@dataclass
class ContinuationRecords:
    profile_dir: Path
    interval: float
    profile_rows: list[dict[str, Any]] = field(default_factory=list)
    sample_rows: list[dict[str, Any]] = field(default_factory=list)
    stage_fates: list[dict[str, Any]] = field(default_factory=list)
    cell_fates: list[dict[str, Any]] = field(default_factory=list)
    summary_rows: list[dict[str, Any]] = field(default_factory=list)

    def execute(self, stage: dict[str, Any]) -> dict[str, Any]:
        result, samples = run_stage(stage, self.profile_dir, self.interval)
        self.profile_rows.append(result)
        self.sample_rows.extend(samples)
        write_profile_table(
            self.profile_dir / "stages.tsv", STAGE_FIELDS, self.profile_rows
        )
        write_profile_table(
            self.profile_dir / "samples.tsv", SAMPLE_FIELDS, self.sample_rows
        )
        return result


def record_tidehunter_failure(
    records: ContinuationRecords,
    seed: int,
    setting: str,
    tidehunter_name: str,
    clustering_name: str,
    tidehunter_time: dict[str, float | int] | None,
    exit_code: int,
    source: str,
    parent: bool = False,
) -> None:
    status = (
        "external_resource_failure_parent_v1" if parent else "external_resource_failure"
    )
    records.summary_rows.append(
        blank_summary(
            seed,
            setting,
            status,
            "parent_inner_exit_137_accuracy_unavailable_not_zero"
            if parent
            else "tidehunter_failed_accuracy_unavailable_not_zero",
            tidehunter=tidehunter_time,
        )
    )
    inner_exit = (
        tidehunter_time["exit_status"] if tidehunter_time is not None else exit_code
    )
    records.cell_fates.append(
        {
            "seed": seed,
            "setting": setting,
            "status": status,
            "tidehunter_status": (
                "imported_failed_outer_exit_1_inner_exit_137"
                if parent
                else f"failed_exit_{inner_exit}"
                + ("_gnu_time_unavailable" if tidehunter_time is None else "")
            ),
            "clustering_status": "not_started_dependency_failure",
            "evaluation_status": "not_started_dependency_failure",
            "source": source,
            "reason": "tidehunter_container_resource_failure"
            if parent
            else "tidehunter_external_process_failure",
        }
    )
    records.stage_fates.extend(
        [
            _stage_fate(
                tidehunter_name,
                seed,
                setting,
                "tidehunter",
                "failed",
                source,
                exit_code,
                "container_inner_exit_137" if parent else "external_process_failure",
            ),
            _stage_fate(
                clustering_name,
                seed,
                setting,
                "clustering",
                "not_started_dependency_failure",
                "continuation_v2",
                "",
                "parent_tidehunter_failed" if parent else "tidehunter_failed",
            ),
        ]
    )


def record_clustering_failure(
    records: ContinuationRecords,
    seed: int,
    setting: str,
    clustering_name: str,
    tidehunter_time: dict[str, float | int],
    clustering_time: dict[str, float | int] | None,
    exit_code: int,
) -> None:
    records.summary_rows.append(
        blank_summary(
            seed,
            setting,
            "external_resource_failure",
            "clustering_failed_accuracy_unavailable_not_zero",
            tidehunter_time,
            clustering_time,
        )
    )
    inner_exit = (
        clustering_time["exit_status"] if clustering_time is not None else exit_code
    )
    records.cell_fates.append(
        {
            "seed": seed,
            "setting": setting,
            "status": "external_resource_failure",
            "tidehunter_status": "ok",
            "clustering_status": f"failed_exit_{inner_exit}"
            + ("_gnu_time_unavailable" if clustering_time is None else ""),
            "evaluation_status": "not_started_dependency_failure",
            "source": "continuation_v2",
            "reason": "clustering_external_process_failure",
        }
    )
    records.stage_fates.append(
        _stage_fate(
            clustering_name,
            seed,
            setting,
            "clustering",
            "failed",
            "continuation_v2",
            exit_code,
            "external_process_failure",
        )
    )


def record_evaluation_failure(
    records: ContinuationRecords,
    seed: int,
    setting: str,
    error: Exception,
    tidehunter_time: dict[str, float | int],
    clustering_time: dict[str, float | int],
) -> None:
    reason = f"{type(error).__name__}:{error}"
    records.summary_rows.append(
        blank_summary(
            seed,
            setting,
            "normalization_or_evaluation_failure",
            reason,
            tidehunter_time,
            clustering_time,
        )
    )
    records.cell_fates.append(
        {
            "seed": seed,
            "setting": setting,
            "status": "normalization_or_evaluation_failure",
            "tidehunter_status": "ok",
            "clustering_status": "ok",
            "evaluation_status": "failed",
            "source": "continuation_v2",
            "reason": reason,
        }
    )


def record_success(
    records: ContinuationRecords,
    seed: int,
    setting: str,
    metrics: dict[str, Any],
    tidehunter_time: dict[str, float | int],
    clustering_time: dict[str, float | int],
) -> None:
    records.summary_rows.append(
        success_summary(seed, setting, metrics, tidehunter_time, clustering_time)
    )
    records.cell_fates.append(
        {
            "seed": seed,
            "setting": setting,
            "status": "ok",
            "tidehunter_status": "ok",
            "clustering_status": "ok",
            "evaluation_status": "ok",
            "source": "continuation_v2",
            "reason": "",
        }
    )


def execute_cells(
    runs: list[dict[str, object]],
    stage_by_name: dict[str, dict[str, object]],
    frozen_skip: set[str],
    parent_time: dict[str, float | int],
    profile_dir: Path,
    interval: float,
) -> ContinuationRecords:
    records = ContinuationRecords(profile_dir, interval)
    for run in runs:
        seed = int(run["seed"])
        setting = str(run["setting"])
        run_dir = Path(run["run_dir"])
        prefix = f"s{seed}_{setting}"
        tidehunter_name = f"{prefix}_tidehunter"
        clustering_name = f"{prefix}_clustering"
        if tidehunter_name in frozen_skip:
            record_tidehunter_failure(
                records,
                seed,
                setting,
                tidehunter_name,
                clustering_name,
                parent_time,
                1,
                "parent_v1_hash_verified",
                parent=True,
            )
            continue

        tidehunter_result = records.execute(stage_by_name[tidehunter_name])
        tidehunter_time = retained_gnu_time(run_dir / "tidehunter.gnu_time.txt")
        if int(tidehunter_result["exit_code"]) != 0:
            record_tidehunter_failure(
                records,
                seed,
                setting,
                tidehunter_name,
                clustering_name,
                tidehunter_time,
                int(tidehunter_result["exit_code"]),
                "continuation_v2",
            )
            continue
        if tidehunter_time is None:
            raise ValueError("successful TideHunter stage lacks a GNU-time record")
        records.stage_fates.append(
            _stage_fate(
                tidehunter_name,
                seed,
                setting,
                "tidehunter",
                "ok",
                "continuation_v2",
                0,
                "",
            )
        )

        clustering_result = records.execute(stage_by_name[clustering_name])
        clustering_time = retained_gnu_time(run_dir / "clustering.gnu_time.txt")
        if int(clustering_result["exit_code"]) != 0:
            record_clustering_failure(
                records,
                seed,
                setting,
                clustering_name,
                tidehunter_time,
                clustering_time,
                int(clustering_result["exit_code"]),
            )
            continue
        if clustering_time is None:
            raise ValueError("successful clustering stage lacks a GNU-time record")
        records.stage_fates.append(
            _stage_fate(
                clustering_name,
                seed,
                setting,
                "clustering",
                "ok",
                "continuation_v2",
                0,
                "",
            )
        )

        try:
            metrics = evaluate_factorial(
                run_dir / "tc_tidehunter.gff3",
                run_dir / "tc_clustering.gff3_1.gff3",
                run_dir / "tc_clustering.gff3",
                run_dir / "tc_consensus/consensus_sequences_all.fasta",
                Path(run["genome_dir"]),
                run_dir / "evaluation",
            )
        except Exception as error:
            record_evaluation_failure(
                records,
                seed,
                setting,
                error,
                tidehunter_time,
                clustering_time,
            )
            continue
        record_success(
            records,
            seed,
            setting,
            metrics,
            tidehunter_time,
            clustering_time,
        )
    return records


def run_continuation(
    base_config_path: Path,
    continuation_config_path: Path,
    outdir: Path,
    docker: str = "docker",
    interval: float = 0.2,
) -> list[dict[str, Any]]:
    if outdir.exists():
        raise FileExistsError(f"Choose a new output directory: {outdir}")
    if not 0.05 <= interval <= 60:
        raise ValueError("sample interval must be between 0.05 and 60 seconds")
    base = json.loads(base_config_path.read_text(encoding="utf-8"))
    continuation = json.loads(continuation_config_path.read_text(encoding="utf-8"))
    if digest_file(base_config_path) != continuation.get("base_config_sha256"):
        raise ValueError("base TideCluster factorial config changed")
    if continuation.get("continuation_rule", {}).get(
        "all_original_datasets_settings_commands_and_accuracy_rules_unchanged"
    ) is not True:
        raise ValueError("continuation did not freeze the original experiment")
    parent = validate_parent_failure(continuation)
    datasets = validate_datasets(base)
    metadata = docker_metadata(docker, base["image"])
    if metadata["image_id"] != base["image_id"]:
        raise ValueError("Docker image ID differs from frozen config")
    project_root = Path(__file__).resolve().parents[2]
    project_commit = command_output(["git", "-C", str(project_root), "rev-parse", "HEAD"])
    project_status = command_output(
        ["git", "-C", str(project_root), "status", "--short", "--untracked-files=no"]
    )
    if project_status:
        raise ValueError("tracked worktree must be clean before continuation execution")

    outdir.mkdir(parents=True)
    stages, runs = build_stage_specs(base, datasets, outdir, docker)
    stage_by_name = {str(stage["name"]): stage for stage in stages}
    if parent["failed_stage"] not in stage_by_name:
        raise ValueError("parent failed stage is absent from the frozen stage plan")
    (outdir / "stage_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "stages": stages,
                "never_rerun_stages": continuation["never_rerun_stages"],
                "failure_policy": continuation["continuation_rule"],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    source_paths = (
        base_config_path,
        continuation_config_path,
        Path(__file__),
        project_root / "benchmarks/tidecluster/factorial_run.py",
        project_root / "benchmarks/tidecluster/factorial_evaluate.py",
        project_root / "benchmarks/tidecluster/normalize.py",
        project_root / "benchmarks/scripts/profile_stage_resources.py",
        project_root / "benchmarks/scripts/run_tidecluster_docker_reference.py",
    )
    for source in source_paths:
        shutil.copyfile(source, snapshot / source.name)
    parent_snapshot = outdir / "parent_failure_snapshot"
    for name in continuation["parent_artifacts"]:
        source = parent["directory"] / name
        target = parent_snapshot / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)

    environment = {
        "schema_version": 1,
        "complete": False,
        "experiment_id": continuation["experiment_id"],
        "project_commit": project_commit,
        "project_status": project_status,
        "base_config_sha256": digest_file(base_config_path),
        "continuation_config_sha256": digest_file(continuation_config_path),
        "tool": metadata,
        "parent_failure": {
            "directory": str(parent["directory"]),
            "failed_stage": parent["failed_stage"],
            "artifacts": continuation["parent_artifacts"],
        },
        "boundary": continuation["boundary"],
    }
    environment_path = outdir / "environment.json"
    environment_path.write_text(json.dumps(environment, indent=2) + "\n", encoding="utf-8")

    profile_dir = outdir / "profile"
    profile_dir.mkdir()
    frozen_skip = set(continuation["never_rerun_stages"])
    parent_time_names = [
        name
        for name in continuation["parent_artifacts"]
        if name.endswith("/tidehunter.gnu_time.txt")
    ]
    if len(parent_time_names) != 1:
        raise ValueError("continuation requires one parent TideHunter GNU-time record")
    parent_time = retained_gnu_time(parent["directory"] / parent_time_names[0])
    if parent_time is None:
        raise ValueError("parent failed-stage GNU-time record is unreadable")
    records = execute_cells(
        runs,
        stage_by_name,
        frozen_skip,
        parent_time,
        profile_dir,
        interval,
    )

    if len(records.summary_rows) != len(runs) or len(records.cell_fates) != len(runs):
        raise RuntimeError("continuation did not record every frozen cell fate")
    if len(records.stage_fates) != len(stages):
        raise RuntimeError("continuation did not record every frozen stage fate")
    write_table(outdir / "summary.tsv", records.summary_rows, SUMMARY_FIELDS)
    write_table(outdir / "cell_fates.tsv", records.cell_fates, CELL_FATE_FIELDS)
    write_table(outdir / "stage_fates.tsv", records.stage_fates, STAGE_FATE_FIELDS)
    profile_receipt = {
        "schema_version": 1,
        "complete": True,
        "planned_stage_count": len(stages),
        "parent_imported_failed_stage_count": 1,
        "continuation_attempted_stage_count": len(records.profile_rows),
        "successful_continuation_stage_count": sum(
            int(row["exit_code"]) == 0 for row in records.profile_rows
        ),
        "failed_continuation_stage_count": sum(
            int(row["exit_code"]) != 0 for row in records.profile_rows
        ),
        "dependency_not_started_stage_count": sum(
            row["status"] == "not_started_dependency_failure"
            for row in records.stage_fates
        ),
        "all_stage_fates_recorded": len(records.stage_fates) == len(stages),
        "sample_interval_seconds": interval,
        "rss_scope": "aggregate_live_host_process_tree;container_gnu_time_in_summary",
    }
    (profile_dir / "receipt.json").write_text(
        json.dumps(profile_receipt, indent=2) + "\n", encoding="utf-8"
    )
    for row in datasets:
        genome_dir = Path(row["genome_dir"])
        for name, expected in (
            ("genome.fa", row["manifest"]["files"]["genome.fa"]),
            ("catalogue.fa", row["manifest"]["files"]["catalogue.fa"]),
            (
                "truth_copy_number.tsv",
                row["manifest"]["files"]["truth_copy_number.tsv"],
            ),
        ):
            if digest_file(genome_dir / name) != expected:
                raise ValueError("a frozen validation input changed during continuation")
    environment["complete"] = True
    environment_path.write_text(json.dumps(environment, indent=2) + "\n", encoding="utf-8")
    successful = sum(row["status"] == "ok" for row in records.summary_rows)
    receipt = {
        "schema_version": 1,
        "experiment_id": continuation["experiment_id"],
        "complete": True,
        "fate": (
            "completed_all_cells_successful"
            if successful == len(records.summary_rows)
            else "completed_with_retained_external_or_technical_failures"
        ),
        "accuracy_complete": successful == len(records.summary_rows),
        "run_count": len(records.summary_rows),
        "successful_run_count": successful,
        "failed_or_unavailable_run_count": len(records.summary_rows) - successful,
        "summary_sha256": digest_file(outdir / "summary.tsv"),
        "cell_fates_sha256": digest_file(outdir / "cell_fates.tsv"),
        "stage_fates_sha256": digest_file(outdir / "stage_fates.tsv"),
        "parent_failure_never_rerun": sorted(frozen_skip),
        "warning": (
            "same_process_simulated_validation_genomes_not_independent_plants;"
            "failed_accuracy_measurements_are_NA_not_zero"
        ),
    }
    (outdir / "run_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return records.summary_rows
