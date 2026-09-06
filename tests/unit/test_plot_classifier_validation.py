import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.plot_classifier_validation import load_inputs


def _metric_row(seed: int, family: str, method: str, outcome: str) -> dict:
    return {
        "seed": seed,
        "unit_substitution_rate": .03,
        "array_fragments": 1,
        "coverage": 5,
        "substitution_rate": .001,
        "assembly_fraction": .5,
        "family_id": family,
        "method": method,
        "outcome": outcome,
    }


def _confusion(outcomes: tuple[str, ...]) -> dict:
    tp = outcomes.count("TP")
    fn = outcomes.count("FN")
    fp = outcomes.count("FP")
    tn = outcomes.count("TN")
    return {
        "available": len(outcomes), "unavailable": 0,
        "TP": tp, "FN": fn, "FP": fp, "TN": tn,
        "sensitivity": tp / (tp + fn),
        "false_positive_rate": fp / (fp + tn),
        "precision": tp / (tp + fp),
    }


def test_classifier_plot_loads_failed_heldout_with_localization(tmp_path: Path) -> None:
    development = tmp_path / "development"
    baseline = tmp_path / "baseline"
    heldout = tmp_path / "heldout"
    for path in (development, baseline, heldout):
        path.mkdir()
    (development / "validation.json").write_text(json.dumps({
        "complete": True,
        "selected_candidate": "blend_a0.5_t0.5",
        "acceptance": {"passed": True},
    }))
    localization = [{
        "seed": 41, "unit_substitution_rate": .03, "array_fragments": 1,
        "assembly_fraction": 1, "family_id": "f1", "true_assembly_bp": 100,
        "predicted_assembly_bp": 98, "overlap_bp": 98, "base_recall": .98,
        "base_precision": 1, "fragments": 1, "truth_fragments": 1,
    }]
    write_table(baseline / "localization_metrics.tsv", localization, list(localization[0]))
    (baseline / "validation.json").write_text(json.dumps({
        "complete": True, "localization_family_rows": 1,
    }))
    baseline_outcomes = ("TP", "FN", "FP", "TN")
    selected_outcomes = ("TP", "TP", "FP", "TN")
    rows = [
        _metric_row(41, f"f{index}", method, outcome)
        for method, outcomes in (
            ("single_k21", baseline_outcomes),
            ("blend_a0.5_t0.5", selected_outcomes),
        )
        for index, outcome in enumerate(outcomes, 1)
    ]
    write_table(heldout / "comparison_metrics.tsv", rows, list(rows[0]))
    (heldout / "environment.json").write_text(json.dumps({
        "classifier_model": "blend_a0.5_t0.5", "heldout_seeds": [41],
    }))
    (heldout / "validation.json").write_text(json.dumps({
        "complete": True,
        "paired_family_conditions": 4,
        "selected_candidate": "blend_a0.5_t0.5",
        "acceptance": {"passed": False},
        "baseline_confusion": _confusion(baseline_outcomes),
        "selected_confusion": _confusion(selected_outcomes),
        "comparison_metrics_sha256": digest_file(heldout / "comparison_metrics.tsv"),
    }))
    loaded, localizations, _, validation = load_inputs(development, baseline, heldout)
    assert len(loaded) == 8 and len(localizations) == 1
    assert validation["acceptance"]["passed"] is False

    rows[-1]["family_id"] = "changed"
    write_table(heldout / "comparison_metrics.tsv", rows, list(rows[0]))
    validation = json.loads((heldout / "validation.json").read_text())
    validation["comparison_metrics_sha256"] = digest_file(heldout / "comparison_metrics.tsv")
    (heldout / "validation.json").write_text(json.dumps(validation))
    try:
        load_inputs(development, baseline, heldout)
    except ValueError as error:
        assert "unpaired" in str(error)
    else:
        raise AssertionError("Unpaired classifier figure input was accepted")
