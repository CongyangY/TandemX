from pathlib import Path

from benchmarks.scripts.run_external_tool_comparison import TruthRecord, parse_tidehunter, parse_trf, score_predictions


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
