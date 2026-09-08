import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.tidecluster.factorial_continue import (
    ContinuationRecords,
    blank_summary,
    execute_cells,
    validate_imported_successes,
    validate_parent_failure,
)


def _parent_config(tmp_path: Path) -> dict:
    parent = tmp_path / "parent"
    (parent / "profile").mkdir(parents=True)
    (parent / "run_receipt.json").write_text(
        json.dumps({"complete": False, "fate": "external_process_failure"})
    )
    (parent / "profile/stages.tsv").write_text(
        "stage\texit_code\n"
        "s6401_default_primary_tidehunter\t1\n"
    )
    artifacts = {
        name: digest_file(parent / name)
        for name in ("run_receipt.json", "profile/stages.tsv")
    }
    return {
        "parent_failure_dir": str(parent),
        "parent_artifacts": artifacts,
        "never_rerun_stages": ["s6401_default_primary_tidehunter"],
    }


def test_validate_parent_failure_requires_exact_never_rerun_stage(
    tmp_path: Path,
) -> None:
    config = _parent_config(tmp_path)
    result = validate_parent_failure(config)
    assert result["failed_stage"] == "s6401_default_primary_tidehunter"
    config["never_rerun_stages"] = []
    with pytest.raises(ValueError, match="never-rerun"):
        validate_parent_failure(config)


def test_validate_parent_failure_rehashes_parent_artifacts(tmp_path: Path) -> None:
    config = _parent_config(tmp_path)
    parent = Path(config["parent_failure_dir"])
    (parent / "run_receipt.json").write_text("changed")
    with pytest.raises(ValueError, match="artifact changed"):
        validate_parent_failure(config)


def test_blank_summary_keeps_failed_resource_and_missing_accuracy() -> None:
    row = blank_summary(
        6401,
        "default_primary",
        "external_resource_failure_parent_v1",
        "accuracy_unavailable_not_zero",
        tidehunter={"wall_seconds": 10.89, "maximum_rss_kb": 7776184},
    )
    assert row["tidehunter_maximum_rss_kb"] == 7776184
    assert row["array_recall"] == ""
    assert row["status"] == "external_resource_failure_parent_v1"


def test_execute_cells_never_reruns_parent_failure_and_continues(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    runs = []
    stages = {}
    for setting in ("default_primary", "matched_period_sensitivity"):
        run_dir = tmp_path / setting
        run_dir.mkdir()
        runs.append(
            {"seed": 6401, "setting": setting, "run_dir": run_dir, "genome_dir": tmp_path}
        )
        for component in ("tidehunter", "clustering"):
            name = f"s6401_{setting}_{component}"
            stages[name] = {"name": name, "cwd": str(run_dir)}

    attempted = []

    def fake_run_stage(stage: dict, _profile_dir: Path, _interval: float):
        attempted.append(stage["name"])
        component = stage["name"].rsplit("_", 1)[-1]
        run_dir = Path(stage["cwd"])
        (run_dir / f"{component}.gnu_time.txt").write_text(
            "User time (seconds): 1.0\n"
            "System time (seconds): 0.2\n"
            "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00\n"
            "Maximum resident set size (kbytes): 1000\n"
            "Exit status: 0\n"
        )
        return (
            {
                "stage": stage["name"],
                "exit_code": 0,
                "wall_seconds": "2.0",
                "user_cpu_seconds": "1.0",
                "system_cpu_seconds": "0.2",
                "peak_process_tree_rss_mib": "1.0",
                "peak_process_count": 1,
                "peak_scratch_bytes": 0,
                "sample_count": 1,
                "stdout_log": "stdout.log",
                "stderr_log": "stderr.log",
            },
            [],
        )

    monkeypatch.setattr(
        "benchmarks.tidecluster.factorial_continue.run_stage", fake_run_stage
    )
    monkeypatch.setattr(
        "benchmarks.tidecluster.factorial_continue.evaluate_factorial",
        lambda *_args: {"array_recall": 0.5, "warning": "test"},
    )
    (tmp_path / "profile").mkdir()
    records = execute_cells(
        runs,
        stages,
        {"s6401_default_primary_tidehunter"},
        {
            "wall_seconds": 10.89,
            "maximum_rss_kb": 7776184,
            "exit_status": 1,
            "user_seconds": 1.0,
            "system_seconds": 1.0,
        },
        {},
        "test_continuation",
        tmp_path / "profile",
        0.2,
    )
    assert attempted == [
        "s6401_matched_period_sensitivity_tidehunter",
        "s6401_matched_period_sensitivity_clustering",
    ]
    assert [row["status"] for row in records.summary_rows] == [
        "external_resource_failure_parent_v1",
        "ok",
    ]
    assert len(records.stage_fates) == 4


def test_continuation_records_normalize_manifest_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    observed = {}

    def fake_run_stage(stage: dict, _profile_dir: Path, _interval: float):
        observed.update(stage)
        return (
            {
                "stage": stage["name"],
                "exit_code": 0,
                "wall_seconds": "0",
                "user_cpu_seconds": "0",
                "system_cpu_seconds": "0",
                "peak_process_tree_rss_mib": "0",
                "peak_process_count": 1,
                "peak_scratch_bytes": 0,
                "sample_count": 0,
                "stdout_log": "stdout",
                "stderr_log": "stderr",
            },
            [],
        )

    monkeypatch.setattr(
        "benchmarks.tidecluster.factorial_continue.run_stage", fake_run_stage
    )
    profile = tmp_path / "profile"
    profile.mkdir()
    ContinuationRecords(profile, 0.2).execute(
        {"name": "stage", "cwd": str(tmp_path), "scratch_dir": str(tmp_path)}
    )
    assert observed["cwd"] == tmp_path
    assert observed["scratch_dir"] == tmp_path


def test_validate_imported_success_rehashes_and_parses_resources(tmp_path: Path) -> None:
    source = tmp_path / "prior"
    run_dir = source / "seed6401/matched_period_sensitivity"
    run_dir.mkdir(parents=True)
    receipt = source / "run_receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "fate": "orchestration_failure_after_external_stage_start",
                "affected_stage": "s6401_matched_period_sensitivity_tidehunter",
            }
        )
    )
    time_path = run_dir / "tidehunter.gnu_time.txt"
    time_path.write_text(
        "User time (seconds): 1.0\n"
        "System time (seconds): 0.2\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00\n"
        "Maximum resident set size (kbytes): 3250036\n"
        "Exit status: 0\n"
    )
    output = run_dir / "tc_tidehunter.gff3"
    output.write_text("##gff-version 3\n")
    names = [
        "run_receipt.json",
        "seed6401/matched_period_sensitivity/tidehunter.gnu_time.txt",
        "seed6401/matched_period_sensitivity/tc_tidehunter.gff3",
    ]
    config = {
        "imported_successful_stages": {
            "s6401_matched_period_sensitivity_tidehunter": {
                "source_dir": str(source),
                "artifacts": {name: digest_file(source / name) for name in names},
                "copy_outputs": names[1:],
            }
        }
    }
    result = validate_imported_successes(config)
    imported = result["s6401_matched_period_sensitivity_tidehunter"]
    assert imported["resources"]["maximum_rss_kb"] == 3250036
    output.write_text("changed\n")
    with pytest.raises(ValueError, match="artifact changed"):
        validate_imported_successes(config)


def test_execute_cells_imports_success_without_rerunning_tidehunter(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "matched_period_sensitivity"
    run_dir.mkdir()
    stages = {
        "s6401_matched_period_sensitivity_tidehunter": {
            "name": "s6401_matched_period_sensitivity_tidehunter",
            "cwd": str(run_dir),
        },
        "s6401_matched_period_sensitivity_clustering": {
            "name": "s6401_matched_period_sensitivity_clustering",
            "cwd": str(run_dir),
        },
    }
    attempted = []

    def fake_run_stage(stage: dict, _profile_dir: Path, _interval: float):
        attempted.append(stage["name"])
        (run_dir / "clustering.gnu_time.txt").write_text(
            "User time (seconds): 1.0\n"
            "System time (seconds): 0.2\n"
            "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00\n"
            "Maximum resident set size (kbytes): 1000\n"
            "Exit status: 0\n"
        )
        return (
            {
                "stage": stage["name"],
                "exit_code": 0,
                "wall_seconds": "2.0",
                "user_cpu_seconds": "1.0",
                "system_cpu_seconds": "0.2",
                "peak_process_tree_rss_mib": "1.0",
                "peak_process_count": 1,
                "peak_scratch_bytes": 0,
                "sample_count": 1,
                "stdout_log": "stdout.log",
                "stderr_log": "stderr.log",
            },
            [],
        )

    monkeypatch.setattr(
        "benchmarks.tidecluster.factorial_continue.run_stage", fake_run_stage
    )
    monkeypatch.setattr(
        "benchmarks.tidecluster.factorial_continue.evaluate_factorial",
        lambda *_args: {"array_recall": 0.5, "warning": "test"},
    )
    profile = tmp_path / "profile"
    profile.mkdir()
    imported = {
        "s6401_matched_period_sensitivity_tidehunter": {
            "resources": {
                "wall_seconds": 11.08,
                "maximum_rss_kb": 3250036,
                "exit_status": 0,
                "user_seconds": 32.99,
                "system_seconds": 1.13,
            }
        }
    }
    records = execute_cells(
        [
            {
                "seed": 6401,
                "setting": "matched_period_sensitivity",
                "run_dir": run_dir,
                "genome_dir": tmp_path,
            }
        ],
        stages,
        {"s6401_matched_period_sensitivity_tidehunter"},
        imported["s6401_matched_period_sensitivity_tidehunter"]["resources"],
        imported,
        "test_continuation",
        profile,
        0.2,
    )
    assert attempted == ["s6401_matched_period_sensitivity_clustering"]
    assert records.summary_rows[0]["status"] == "ok"
    assert records.stage_fates[0]["execution_source"] == (
        "prior_continuation_hash_verified"
    )
