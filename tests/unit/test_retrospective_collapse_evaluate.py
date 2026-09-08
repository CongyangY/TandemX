import json
from pathlib import Path

import pytest

from benchmarks.retrospective.evaluate import (
    EvaluationConfig,
    assembly_bp_by_family,
    evaluate_families,
    run_evaluation,
    summarize_rows,
)


def test_assembly_bp_unions_overlaps_per_chromosome(tmp_path: Path) -> None:
    bed = tmp_path / "arrays.bed"
    bed.write_text(
        "chr1\t0\t100\tf1\n"
        "chr1\t50\t150\tf1\n"
        "chr2\t0\t25\tf1\n"
        "chr1\t20\t40\tf2\n"
    )
    assert assembly_bp_by_family(bed) == {"f1": 175.0, "f2": 20.0}


def test_evaluation_keeps_eligible_ineligible_and_failure_rows() -> None:
    config = EvaluationConfig(primary_min_new_bp=100, min_new_bp_sensitivity=(50, 100, 200))
    rows = evaluate_families(
        read_bp={"tp": 200, "tn": 100, "small": 50},
        old_bp={"tp": 20, "tn": 90, "small": 0, "failure": 10},
        new_bp={"tp": 190, "tn": 100, "small": 80, "failure": 200},
        sensitivity_bp=None,
        config=config,
    )
    by_family = {row["family_id"]: row for row in rows}
    assert by_family["tp"]["outcome"] == "TP"
    assert by_family["tn"]["outcome"] == "TN"
    assert by_family["small"]["eligibility"] == "not_source_eligible"
    assert by_family["failure"]["eligibility"] == "technical_failure"
    assert len(rows) == 4
    summary = summarize_rows(rows, config)
    assert summary["confusion"] == {"TP": 1, "FN": 0, "FP": 0, "TN": 1}
    assert summary["eligible_family_rows"] == 2
    assert summary["not_source_eligible_rows"] == 1
    assert summary["technical_failure_rows"] == 1
    assert summary["min_new_bp_sensitivity"][0] == {
        "min_new_bp": 50, "families": 3, "TP": 2, "FN": 0, "FP": 0, "TN": 1
    }


def test_source_eligibility_does_not_depend_on_new_read_agreement() -> None:
    config = EvaluationConfig(primary_min_new_bp=100, min_new_bp_sensitivity=(100,))
    rows = evaluate_families(
        read_bp={"discordant": 10_000},
        old_bp={"discordant": 10},
        new_bp={"discordant": 100},
        sensitivity_bp=None,
        config=config,
    )
    assert rows[0]["new_read_ratio"] == 0.01
    assert rows[0]["eligibility"] == "eligible"
    assert rows[0]["outcome"] == "TP"


def test_run_writes_complete_outputs_and_refuses_overwrite(tmp_path: Path) -> None:
    copy_number = tmp_path / "copy_number.tsv"
    copy_number.write_text("family_id\testimated_bp\nf1\t200\n")
    old = tmp_path / "old.bed"
    old.write_text("chr1\t0\t20\tf1\n")
    new = tmp_path / "new.bed"
    new.write_text("chr1\t0\t190\tf1\n")
    outdir = tmp_path / "result"
    config = EvaluationConfig(primary_min_new_bp=100, min_new_bp_sensitivity=(100,))
    rows, summary = run_evaluation(copy_number, old, new, outdir, config)
    assert rows[0]["outcome"] == "TP"
    assert json.loads((outdir / "summary.json").read_text())["confusion"]["TP"] == 1
    assert len((outdir / "family_metrics.tsv").read_text().splitlines()) == 2
    with pytest.raises(FileExistsError):
        run_evaluation(copy_number, old, new, outdir, config)


def test_config_requires_primary_threshold_in_sensitivity_set() -> None:
    with pytest.raises(ValueError, match="included"):
        EvaluationConfig(primary_min_new_bp=100, min_new_bp_sensitivity=(50,)).validate()
