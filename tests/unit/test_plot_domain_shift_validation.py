import csv
import json
from pathlib import Path

from benchmarks.scripts.plot_domain_shift_validation import load_inputs


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_domain_shift_plot_requires_exact_scenario_pairing(tmp_path: Path) -> None:
    metrics = tmp_path/"metrics.tsv"
    validation = tmp_path/"validation.json"
    localization = tmp_path/"localization.tsv"
    baseline_validation = tmp_path/"baseline_validation.json"
    rows = []
    counts = {}
    for method, outcome in (("single_k21", "FN"), ("multik_fallback_depth_rule", "TP")):
        rows.append({
            "seed": 1, "unit_substitution_rate": .03, "array_fragments": 3,
            "coverage": 5, "substitution_rate": .001, "assembly_fraction": .5,
            "family_id": "f1", "method": method, "outcome": outcome,
        })
        counts[method] = {
            "TP": int(outcome == "TP"), "FN": int(outcome == "FN"), "FP": 0, "TN": 0,
        }
    write_tsv(metrics, rows)
    validation.write_text(json.dumps({
        "complete": True, "split": "heldout", "method_confusion": counts,
    }))
    write_tsv(localization, [{
        "seed": 1, "unit_substitution_rate": .03, "array_fragments": 3,
        "assembly_fraction": .5, "family_id": "f1", "base_recall": .8,
        "base_precision": 1, "fragments": 2, "truth_fragments": 3,
    }])
    baseline_validation.write_text(json.dumps({
        "complete": True, "localization_family_rows": 1,
    }))
    selected, locations, _ = load_inputs(metrics, validation, localization, baseline_validation)
    assert len(selected) == 2 and len(locations) == 1

    rows[1]["array_fragments"] = 1
    write_tsv(metrics, rows)
    try:
        load_inputs(metrics, validation, localization, baseline_validation)
    except ValueError as error:
        assert "paired" in str(error)
    else:
        raise AssertionError("Unpaired domain-shift rows were accepted")
