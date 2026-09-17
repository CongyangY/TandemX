"""Score TandemX locate with supplied truth-derived monomers on one control."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


BASE = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    truth_path = BASE / "control/truth_monomers.tsv"
    prediction_path = BASE / "native/tandemx_catalogue_supplied/arrays.bed"
    truth = [(int(r["start0"]), int(r["end0"]), r["label"])
             for r in csv.DictReader(truth_path.open(), delimiter="\t")]
    with prediction_path.open() as source:
        calls = [(int(r[1]), int(r[2]), r[3]) for r in csv.reader(source, delimiter="\t")]
    truth_bp = set().union(*(set(range(start, end)) for start, end, _ in truth))
    called_bp = set().union(*(set(range(start, end)) for start, end, _ in calls))
    overlap = len(truth_bp & called_bp)
    near_boundary_label_matches = sum(
        any(label == predicted_label and abs(start - predicted_start) <= 10
            and abs(end - predicted_end) <= 10
            for predicted_start, predicted_end, predicted_label in calls)
        for start, end, label in truth
    )
    result = {
        "condition": "supplied_truth_derived_candidate_monomer_catalogue",
        "comparison_status": "not_rankable_against_de_novo_TideCluster_or_CENdetectHOR",
        "task": "assembly_only_monomer_localization_not_HOR_order",
        "truth_calls": len(truth),
        "predicted_calls": len(calls),
        "truth_array_bp": len(truth_bp),
        "called_array_bp": len(called_bp),
        "overlap_bp": overlap,
        "array_bp_recall": overlap / len(truth_bp),
        "array_bp_precision": overlap / len(called_bp),
        "truth_monomer_calls_with_matching_label_and_both_boundaries_within_10bp": near_boundary_label_matches,
        "truth_monomer_call_denominator": len(truth),
        "truth_sha256": sha(truth_path),
        "prediction_sha256": sha(prediction_path),
    }
    (BASE / "score_tandemx_catalogue_supplied.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
