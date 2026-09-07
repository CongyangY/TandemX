from __future__ import annotations

import csv
import json
from pathlib import Path

import yaml

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.evaluate_cascade_heldout import evaluate


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _run(
    tmp_path: Path,
    tx_time: float = 1.5,
    evaluation_split: str = "heldout",
    extra_metric_gates: bool = False,
) -> tuple[Path, Path]:
    evaluation_seed = 3 if evaluation_split == "heldout" else 4
    config = {
        "benchmark_id": "test", "repetitions": 2,
        "seeds": {"development": [1], "validation": [4], "heldout": [3]},
        "tools": {"tandemx": "x", "tidehunter": "y", "trf": "z"},
        "scenarios": [{"name": "positive"}, {"name": "control", "positive_fraction": 0.0}],
        "acceptance_gates": {
            "tandemx_failed_runs_max": 0, "all_comparator_failed_runs_max": 0,
            "tandemx_deterministic_groups_fraction_min": 1.0,
            "paired_tidehunter_groups_fraction_min": 1.0,
            "positive_array_recall_min": 0.95, "positive_array_precision_min": 0.95,
            "negative_read_call_rate_max": 0.0, "related_family_cyclic_monomer_recall_min": 1.0,
            "tidehunter_array_recall_noninferiority_margin": 0.02,
            "tidehunter_array_precision_noninferiority_margin": 0.02,
            "tidehunter_runtime_geometric_mean_ratio_max": 2.0,
            "tidehunter_peak_rss_geometric_mean_ratio_max": 1.25,
        },
    }
    if evaluation_split != "heldout":
        config["evaluation_split"] = evaluation_split
    if extra_metric_gates:
        config["acceptance_gates"].update(
            {
                "positive_base_union_f1_min": 0.99,
                "positive_matched_boundary_mae_bp_max": 4.0,
            }
        )
    # The evaluator requires the named related-family scenario for its dedicated gate.
    config["scenarios"].append({"name": "related_families"})
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    run = tmp_path / "run"
    run.mkdir()
    raw = []
    summary = []
    for scenario in ("positive", "control", "related_families"):
        for tool in ("tandemx", "tidehunter", "trf"):
            for repetition in (1, 2):
                raw.append({"scenario": scenario, "seed": evaluation_seed, "tool": tool, "repetition": repetition, "status": "ok"})
            positive = scenario != "control"
            summary.append({
                "scenario": scenario, "seed": evaluation_seed, "tool": tool,
                "successful_runs": 2, "attempted_runs": 2, "deterministic": True,
                "array_recall": 1.0 if positive else "NA", "array_precision": 1.0 if positive else "NA",
                "negative_read_call_rate": 0.0, "cyclic_monomer_recall": 1.0 if positive else "NA",
                "base_union_f1": 0.995 if positive else "NA",
                "matched_boundary_mae_bp": 2.0 if positive else "NA",
                "median_runtime_seconds": tx_time if tool == "tandemx" else 1.0,
                "median_peak_rss_mib": 10.0,
            })
    _write(run / "raw_runs.tsv", raw)
    _write(run / "summary.tsv", summary)
    (run / "environment.json").write_text(json.dumps({"config_sha256": digest_file(config_path), "split": evaluation_split, "source_digest": "abc"}))
    (run / "run_config.yaml").write_text(yaml.safe_dump({**config, "selected_split": evaluation_split, "selected_scenarios": None}))
    (run / "validation.json").write_text(json.dumps({"complete": True, "total_runs": len(raw)}))
    return config_path, run


def test_heldout_gates_pass_complete_noninferior_matrix(tmp_path: Path) -> None:
    config, run = _run(tmp_path)
    result = evaluate(config, run, tmp_path / "evaluation")
    assert result["status"] == "passed"
    assert result["failed_gate_names"] == []


def test_heldout_gates_retain_performance_failure(tmp_path: Path) -> None:
    config, run = _run(tmp_path, tx_time=3.0)
    result = evaluate(config, run, tmp_path / "evaluation")
    assert result["status"] == "failed"
    assert "tidehunter_runtime_geometric_mean_ratio" in result["failed_gate_names"]


def test_validation_split_and_optional_accuracy_gates(tmp_path: Path) -> None:
    config, run = _run(
        tmp_path, evaluation_split="validation", extra_metric_gates=True
    )
    result = evaluate(config, run, tmp_path / "evaluation")
    assert result["status"] == "passed"
    observed = {row["name"]: row["observed"] for row in result["gates"]}
    assert observed["minimum_positive_base_union_f1"] == 0.995
    assert observed["maximum_positive_matched_boundary_mae_bp"] == 2.0
    assert result["warning"].startswith("validation_seeds_consumed_once")
