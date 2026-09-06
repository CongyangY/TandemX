import csv
import json
from pathlib import Path

from benchmarks.scripts.plot_collapse_validation import confusion, load_pairs


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_confusion_uses_endpoint_denominators() -> None:
    result = confusion([{"outcome": value} for value in ("TP", "TP", "FN", "FP", "TN", "TN")])
    assert result["sensitivity"] == 2 / 3
    assert result["false_positive_rate"] == 1 / 3
    assert result["precision"] == 2 / 3


def test_load_pairs_rejects_unpaired_or_changed_counts(tmp_path: Path) -> None:
    metrics = tmp_path / "metrics.tsv"
    validation = tmp_path / "validation.json"
    rows = []
    for method, outcome in (("single_k21", "FN"), ("multik_fallback_depth_rule", "TP")):
        rows.append({"seed": 1, "coverage": 1, "substitution_rate": 0,
                     "assembly_fraction": 0.5, "family_id": "f1",
                     "method": method, "outcome": outcome})
    write_tsv(metrics, rows)
    validation.write_text(json.dumps({
        "complete": True, "split": "heldout", "method_confusion": {
            "single_k21": {"TP": 0, "FN": 1, "FP": 0, "TN": 0},
            "multik_fallback_depth_rule": {"TP": 1, "FN": 0, "FP": 0, "TN": 0},
        },
    }))
    selected, _ = load_pairs(metrics, validation)
    assert len(selected) == 2
    rows[1]["family_id"] = "f2"
    write_tsv(metrics, rows)
    try:
        load_pairs(metrics, validation)
    except ValueError as error:
        assert "paired" in str(error)
    else:
        raise AssertionError("Unpaired figure inputs were accepted")
