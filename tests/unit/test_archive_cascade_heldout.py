from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
import yaml

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_cascade_heldout import archive, validate_gate_receipt
from benchmarks.scripts.evaluate_cascade_heldout import evaluate


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _failed_comparator_run(tmp_path: Path) -> tuple[Path, Path, Path]:
    config = {
        "benchmark_id": "heldout_archive_test",
        "repetitions": 1,
        "seeds": {"development": [1], "heldout": [3]},
        "tools": {"tandemx": "x", "tidehunter": "y", "trf": "z"},
        "scenarios": [
            {"name": "positive"},
            {"name": "control", "positive_fraction": 0.0},
            {"name": "related_families"},
        ],
        "acceptance_gates": {
            "tandemx_failed_runs_max": 0,
            "all_comparator_failed_runs_max": 0,
            "tandemx_deterministic_groups_fraction_min": 1.0,
            "paired_tidehunter_groups_fraction_min": 1.0,
            "positive_array_recall_min": 0.95,
            "positive_array_precision_min": 0.95,
            "negative_read_call_rate_max": 0.0,
            "related_family_cyclic_monomer_recall_min": 1.0,
            "tidehunter_array_recall_noninferiority_margin": 0.02,
            "tidehunter_array_precision_noninferiority_margin": 0.02,
            "tidehunter_runtime_geometric_mean_ratio_max": 2.0,
            "tidehunter_peak_rss_geometric_mean_ratio_max": 1.25,
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    run = tmp_path / "run"
    run.mkdir()
    raw: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    for scenario in ("positive", "control", "related_families"):
        for tool in ("tandemx", "tidehunter", "trf"):
            failed = scenario == "control" and tool == "trf"
            raw.append(
                {
                    "scenario": scenario,
                    "dataset_id": f"{scenario}_s3",
                    "seed": 3,
                    "tool": tool,
                    "repetition": 1,
                    "status": "failed" if failed else "ok",
                }
            )
            positive = scenario != "control"
            summary.append(
                {
                    "scenario": scenario,
                    "seed": 3,
                    "tool": tool,
                    "successful_runs": 0 if failed else 1,
                    "attempted_runs": 1,
                    "deterministic": not failed,
                    "array_recall": 1.0 if positive else "NA",
                    "array_precision": 1.0 if positive else "NA",
                    "negative_read_call_rate": 0.0 if not failed else "NA",
                    "cyclic_monomer_recall": 1.0 if positive else "NA",
                    "median_runtime_seconds": 1.5 if tool == "tandemx" else 1.0,
                    "median_peak_rss_mib": 10.0,
                }
            )
    _write(run / "raw_runs.tsv", raw)
    _write(run / "summary.tsv", summary)
    (run / "environment.json").write_text(
        json.dumps(
            {
                "config_sha256": digest_file(config_path),
                "split": "heldout",
                "source_digest": "abc",
            }
        )
    )
    (run / "run_config.yaml").write_text(
        yaml.safe_dump(
            {**config, "selected_split": "heldout", "selected_scenarios": None}
        )
    )
    (run / "validation.json").write_text(
        json.dumps({"complete": True, "total_runs": len(raw)})
    )
    (run / "run.log").write_text("complete with retained comparator failure\n")
    failed_dir = run / "runs" / "control_s3" / "trf" / "rep1"
    failed_dir.mkdir(parents=True)
    (failed_dir / "command.json").write_text('{"command": ["trf"]}\n')
    (failed_dir / "receipt.json").write_text('{"timed_out": true}\n')
    (failed_dir / "stderr.log").write_text("timeout\n")
    evaluation = tmp_path / "evaluation"
    result = evaluate(config_path, run, evaluation)
    assert result["status"] == "failed"
    return config_path, run, evaluation


def test_archive_retains_failed_comparator_and_gate(tmp_path: Path) -> None:
    config, run, evaluation = _failed_comparator_run(tmp_path)
    destination = tmp_path / "archive"
    result = archive(config, run, evaluation, destination)
    assert result["gate_status"] == "failed"
    assert result["failed_run_count"] == 1
    assert result["failed_run_tools"] == ["trf"]
    assert "all_comparator_failed_runs" in result["failed_gate_names"]
    assert (destination / "failed_runs/control_s3/trf/rep1/receipt.json").is_file()
    assert (destination / "archive_manifest.json").is_file()


def test_archive_rejects_inconsistent_gate_receipt(tmp_path: Path) -> None:
    config_path, run, evaluation = _failed_comparator_run(tmp_path)
    receipt_path = evaluation / "gate_results.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["failed_gate_names"] = []
    receipt_path.write_text(json.dumps(receipt))
    config = yaml.safe_load(config_path.read_text())
    with pytest.raises(ValueError, match="internally inconsistent"):
        validate_gate_receipt(config, config_path, run, evaluation)


def test_archive_supports_successful_validation_without_failures(
    tmp_path: Path,
) -> None:
    evidence = Path("paper/evidence/cascade_gap_free_validation_v1")
    destination = tmp_path / "archive"
    result = archive(
        evidence / "frozen_config.yaml",
        evidence / "run",
        evidence / "evaluation",
        destination,
    )
    assert result["gate_status"] == "passed"
    assert result["failed_run_count"] == 0
    assert result["evaluation_split"] == "validation"
    assert result["evaluation_seeds"] == ["2201"]
    with (destination / "failed_runs.tsv").open(newline="", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []
    assert "2201 was consumed once" in (destination / "README.md").read_text()
