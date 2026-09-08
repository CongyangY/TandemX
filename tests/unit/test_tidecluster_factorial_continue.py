import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.tidecluster.factorial_continue import (
    blank_summary,
    execute_cells,
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
