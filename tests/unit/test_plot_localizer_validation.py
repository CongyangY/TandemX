import csv
import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.plot_localizer_validation import load_inputs


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _localization(root: Path, seed: int, split: str) -> None:
    root.mkdir()
    rows = [{
        "seed": seed, "unit_substitution_rate": .05, "array_fragments": 3,
        "assembly_fraction": 1, "family_id": "f1", "true_assembly_bp": 100,
        "predicted_assembly_bp": 98, "base_recall": .98, "base_precision": 1,
        "fragments": 3, "truth_fragments": 3,
    }]
    _write_tsv(root / "localization_metrics.tsv", rows)
    (root / "run_config.json").write_text(json.dumps({"seeds": {split: [seed]}}))
    (root / "validation.json").write_text(json.dumps({
        "complete": True, "executions": 1, "successful": 1,
        "localization_family_rows": 1, "comparison_family_rows": 1,
    }))
    (root / "environment.json").write_text(json.dumps({"split": split}))


def test_localizer_plot_inputs_require_paired_development_and_fresh_heldout(tmp_path: Path) -> None:
    development_v1 = tmp_path / "development_v1"
    development_v2 = tmp_path / "development_v2"
    heldout = tmp_path / "heldout"
    multik = tmp_path / "multik"
    _localization(development_v1, 31, "development")
    _localization(development_v2, 31, "development")
    _localization(heldout, 41, "heldout")
    multik.mkdir()
    rows = []
    confusion = {}
    for method, outcome in (("single_k21", "FN"), ("multik_fallback_depth_rule", "TP")):
        rows.append({
            "seed": 41, "unit_substitution_rate": .05, "array_fragments": 3,
            "coverage": 5, "substitution_rate": .001, "assembly_fraction": .5,
            "family_id": "f1", "method": method, "outcome": outcome,
        })
        confusion[method] = {
            "TP": int(outcome == "TP"), "FN": int(outcome == "FN"), "FP": 0, "TN": 0,
        }
    _write_tsv(multik / "comparison_metrics.tsv", rows)
    (multik / "validation.json").write_text(json.dumps({
        "complete": True, "split": "heldout", "method_confusion": confusion,
    }))
    (multik / "run_config.json").write_text("{}")
    (multik / "environment.json").write_text(json.dumps({
        "previous_validation_sha256": digest_file(heldout / "validation.json"),
        "previous_environment_sha256": digest_file(heldout / "environment.json"),
        "previous_config_sha256": digest_file(heldout / "run_config.json"),
    }))
    localizations, _, selected, _ = load_inputs(
        development_v1, development_v2, heldout, multik
    )
    assert len(localizations) == 3 and len(selected) == 2

    with (development_v2 / "localization_metrics.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[0]["true_assembly_bp"] = "99"
    _write_tsv(development_v2 / "localization_metrics.tsv", rows)
    try:
        load_inputs(development_v1, development_v2, heldout, multik)
    except ValueError as error:
        assert "paired" in str(error)
    else:
        raise AssertionError("Unpaired development truth was accepted")
