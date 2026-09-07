from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.scripts.archive_cascade_evidence import compare_to_baseline, validate_run


def _table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _run(tmp_path: Path) -> Path:
    run = tmp_path / "run"
    run.mkdir()
    (run / "environment.json").write_text(json.dumps({"source_digest": "a" * 64}))
    (run / "validation.json").write_text(json.dumps({"complete": True, "failed_runs": 0, "total_runs": 3}))
    (run / "run_config.yaml").write_text("discovery_method: cascade\n")
    (run / "run.log").write_text("complete\n")
    raw = [{"scenario": "clean", "repetition": i, "status": "ok", "exit_code": 0,
            "timed_out": False, "prediction_sha256": "p", "catalog_sha256": "c"}
           for i in range(1, 4)]
    _table(run / "raw_runs.tsv", raw)
    _table(run / "summary.tsv", [{
        "scenario": "clean", "successful_runs": 3, "attempted_runs": 3,
        "deterministic": True, "median_runtime_seconds": 1, "median_peak_rss_mib": 8,
        "array_recall": 1, "array_precision": 1, "array_f1": 1,
        "read_detection_recall": 1, "read_detection_precision": 1,
        "negative_read_call_rate": 0, "sequence_family_recall": 1,
        "matched_period_mae_bp": 0, "matched_boundary_mae_bp": 0,
    }])
    return run


def test_cascade_archive_validation_and_comparison(tmp_path: Path) -> None:
    run = validate_run(_run(tmp_path))
    baseline = tmp_path / "baseline.tsv"
    old = dict(run["summary"][0])
    old.update(median_runtime_seconds=2, median_peak_rss_mib=10)
    _table(baseline, [old])
    comparison, aggregate = compare_to_baseline(baseline, [run])
    assert comparison[0]["runtime_ratio"] == 0.5
    assert aggregate["runtime_faster_scenarios"] == 1


def test_cascade_archive_rejects_nondeterminism_and_regression(tmp_path: Path) -> None:
    path = _run(tmp_path)
    rows = list(csv.DictReader((path / "raw_runs.tsv").open(), delimiter="\t"))
    rows[1]["prediction_sha256"] = "different"
    _table(path / "raw_runs.tsv", rows)
    with pytest.raises(ValueError, match="Non-deterministic"):
        validate_run(path)
