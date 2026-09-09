import json
from pathlib import Path

from benchmarks.retrospective.evaluate import EvaluationConfig, run_evaluation
from benchmarks.scripts.independent_verify_donor_matched_collapse import verify


def make_fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"evaluation": {
        "collapse_threshold": 0.6,
        "overexpansion_threshold": 1.5,
        "primary_min_new_bp": 100,
        "min_new_bp_sensitivity": [50, 100, 200],
    }}))
    copy_number = tmp_path / "copy_number.tsv"
    copy_number.write_text("family_id\testimated_bp\ncollapse\t200\nretained\t100\n")
    old = tmp_path / "old.bed"
    old.write_text("chr1\t0\t20\tcollapse\nchr1\t0\t90\tretained\n")
    new = tmp_path / "new.bed"
    new.write_text("chr1\t0\t190\tcollapse\nchr1\t0\t100\tretained\n")
    sensitivity = tmp_path / "sensitivity.bed"
    sensitivity.write_text("chr1\t0\t195\tcollapse\nchr1\t0\t100\tretained\n")
    result_dir = tmp_path / "result"
    run_evaluation(
        copy_number,
        old,
        new,
        result_dir,
        EvaluationConfig(primary_min_new_bp=100, min_new_bp_sensitivity=(50, 100, 200)),
        sensitivity,
    )
    return config, copy_number, old, new, sensitivity, result_dir / "family_metrics.tsv", result_dir / "summary.json"


def test_independent_verifier_recomputes_complete_evaluation(tmp_path: Path) -> None:
    config, copy_number, old, new, sensitivity, metrics, summary = make_fixture(tmp_path)
    receipt = verify(
        config_path=config,
        copy_number_path=copy_number,
        old_arrays_path=old,
        new_arrays_path=new,
        sensitivity_arrays_path=sensitivity,
        family_metrics_path=metrics,
        summary_path=summary,
    )
    assert receipt["verification_passed"] is True
    assert receipt["family_rows_recomputed"] == 2
    assert all(receipt["checks"].values())


def test_independent_verifier_detects_tampered_family_metric(tmp_path: Path) -> None:
    config, copy_number, old, new, sensitivity, metrics, summary = make_fixture(tmp_path)
    metrics.write_text(metrics.read_text().replace("\t20.0\t", "\t21.0\t", 1))
    receipt = verify(
        config_path=config,
        copy_number_path=copy_number,
        old_arrays_path=old,
        new_arrays_path=new,
        sensitivity_arrays_path=sensitivity,
        family_metrics_path=metrics,
        summary_path=summary,
    )
    assert receipt["verification_passed"] is False
    assert "family_numeric_mismatch:collapse:old_assembly_bp" in receipt["failures"]


def test_independent_verifier_accepts_absent_sensitivity_assembly(tmp_path: Path) -> None:
    config, copy_number, old, new, _sensitivity, _metrics, _summary = make_fixture(tmp_path)
    result_dir = tmp_path / "result_without_sensitivity"
    run_evaluation(
        copy_number,
        old,
        new,
        result_dir,
        EvaluationConfig(primary_min_new_bp=100, min_new_bp_sensitivity=(50, 100, 200)),
    )
    receipt = verify(
        config_path=config,
        copy_number_path=copy_number,
        old_arrays_path=old,
        new_arrays_path=new,
        sensitivity_arrays_path=None,
        family_metrics_path=result_dir / "family_metrics.tsv",
        summary_path=result_dir / "summary.json",
    )

    assert receipt["verification_passed"] is True
    assert "sensitivity_arrays" not in receipt["inputs"]
