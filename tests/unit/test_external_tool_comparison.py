from pathlib import Path

import sys

from benchmarks.scripts.run_external_tool_comparison import (
    TruthRecord,
    aggregate,
    command_for_tool,
    filter_period_range,
    parse_tidehunter,
    parse_trf,
    prediction_digest,
    run_logged,
    score_predictions,
)


def test_parse_tidehunter(tmp_path: Path) -> None:
    path = tmp_path / "tidehunter.tsv"
    path.write_text("repeat_1 rep0 5 1000 1 1000 60 99.0 0 1,61 ACGT\n", encoding="utf-8")
    assert parse_tidehunter(path) == {"repeat_1": [60]}


def test_parse_trf_ngs(tmp_path: Path) -> None:
    path = tmp_path / "trf.tsv"
    path.write_text("@repeat_1\n1 1000 60 16.7 60 99 0 100 25 25 25 25 2.0\n", encoding="utf-8")
    assert parse_trf(path) == {"repeat_1": [60]}


def test_score_predictions() -> None:
    truth = {"repeat_1": TruthRecord("repeat_1", True, 60), "repeat_2": TruthRecord("repeat_2", True, 60), "background_1": TruthRecord("background_1", False, 0)}
    scores = score_predictions({"repeat_1": [61], "background_1": [30]}, truth)
    assert scores["recall"] == 0.5
    assert scores["precision"] == 0.5
    assert scores["false_positive_rate"] == 1.0
    assert scores["monomer_length_mae_bp"] == 1.0


def test_filter_period_range_and_digest_are_order_stable() -> None:
    first = filter_period_range({"read_b": [60, 10], "read_a": [61, 60]}, 30, 500)
    second = {"read_a": [60, 61], "read_b": [60]}
    assert first == {"read_b": [60], "read_a": [61, 60]}
    assert prediction_digest(first) == prediction_digest(second)


def test_tandemx_comparison_command_is_explicitly_single_threaded(tmp_path: Path) -> None:
    command, _output = command_for_tool(
        "tandemx", Path("tandemx"), Path("reads.fa"), tmp_path, 30, 500, 2000
    )
    assert command[command.index("--threads") + 1] == "1"
    assert "--no-progress" in command


def test_run_logged_records_peak_memory(tmp_path: Path) -> None:
    status, runtime, peak_memory = run_logged(
        [sys.executable, "-c", "print('ok')"], tmp_path / "stdout.log", tmp_path / "stderr.log"
    )
    assert status == 0
    assert runtime > 0
    assert peak_memory is None or peak_memory > 0


def test_aggregate_reports_memory_variability_and_determinism() -> None:
    rows = []
    for repetition, runtime, memory in ((1, 1.0, 10.0), (2, 1.2, 12.0), (3, 1.1, 11.0)):
        rows.append(
            {
                "dataset_id": "data", "tool": "tool", "repetition": repetition,
                "runtime_seconds": runtime, "peak_memory_mb": memory, "prediction_sha256": "same",
                "read_count": 1, "total_bases": 1_000_000, "monomer_length": 60,
                "error_rate": 0.0, "recall": 1.0, "precision": 1.0,
                "false_positive_rate": 0.0, "monomer_length_mae_bp": 0.0, "exit_status": 0,
            }
        )
    summary = aggregate(rows)[0]
    assert summary["median_peak_memory_mb"] == 11.0
    assert summary["deterministic_predictions"] == "true"
    assert summary["unique_prediction_digests"] == 1
    assert float(summary["runtime_cv"]) > 0
