"""Score frozen B1 v2 injected edits without claiming biological accuracy."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from benchmarks.scripts.build_b1_structure_cases import sha256_file, write_jsonl


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def keyed(rows: list[dict], label: str) -> dict[str, dict]:
    result = {}
    for row in rows:
        case_id = row.get("case_id")
        if not isinstance(case_id, str) or case_id in result:
            raise ValueError(f"{label}: missing or duplicate case_id")
        result[case_id] = row
    return result


def edit_distance(a: list[str], b: list[str]) -> int:
    previous = list(range(len(b) + 1))
    for i, left in enumerate(a, 1):
        current = [i]
        for j, right in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (left != right)))
        previous = current
    return previous[-1]


def interval_iou(a: list[int], b: list[int]) -> float:
    if not (isinstance(a, list) and isinstance(b, list) and len(a) == len(b) == 2
            and all(isinstance(x, int) and not isinstance(x, bool) for x in a + b)
            and 0 <= a[0] <= a[1] and 0 <= b[0] <= b[1]):
        raise ValueError("invalid half-open interval")
    union = max(a[1], b[1]) - min(a[0], b[0])
    return 1.0 if union == 0 else max(0, min(a[1], b[1]) - max(a[0], b[0])) / union


def finite_number(value: object) -> bool:
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)


def mean_or_block(values: list[float | None], expected_count: int) -> float | None:
    return sum(values) / expected_count if len(values) == expected_count and all(
        value is not None for value in values) else None


def score(bundle: Path, metrics_path: Path, predictions_path: Path, outdir: Path) -> dict:
    receipt = json.loads((bundle / "receipt.json").read_text())
    metrics = json.loads(metrics_path.read_text())
    if metrics["schema_version"] != 2 or metrics["case_count"] != receipt["case_count"]:
        raise ValueError("frozen metrics and bundle disagree")
    for field, filename in (("inputs_sha256", "inputs.jsonl"), ("truth_sha256", "truth.jsonl")):
        if receipt[field] != sha256_file(bundle / filename):
            raise ValueError(f"bundle {filename} hash mismatch")
    inputs = keyed(read_jsonl(bundle / "inputs.jsonl"), "inputs")
    truth = keyed(read_jsonl(bundle / "truth.jsonl"), "truth")
    predictions = keyed(read_jsonl(predictions_path), "predictions")
    ids = receipt["case_ids"]
    if (len(ids) != metrics["case_count"] or set(ids) != set(inputs) or set(ids) != set(truth)
            or set(predictions) - set(ids)):
        raise ValueError("case ID set mismatch or unknown prediction")
    if any(truth[c]["edited_sequence_sha256"] != inputs[c]["assembly_sequence_sha256"] for c in ids):
        raise ValueError("truth/input edited sequence mismatch")
    equivalence = {}
    for case_id in ids:
        equivalence.setdefault(inputs[case_id]["assembly_sequence_sha256"], []).append(case_id)
    per_case = []
    for case_id in ids:
        target = truth[case_id]
        pred = predictions.get(case_id)
        status = "not_reported" if pred is None else pred.get("status")
        if status not in {"ok", "abstain", "unmapped", "failed", "not_run", "not_reported"}:
            raise ValueError(f"invalid status: {case_id}")
        value = None if pred is None else pred.get("event_score")
        if status == "ok":
            if not finite_number(value) or not 0 <= value <= 1:
                raise ValueError(f"ok requires finite event_score in [0,1]: {case_id}")
        elif value is not None:
            raise ValueError(f"non-ok event_score must be null: {case_id}")
        positive = target["event_status"] == "positive"
        equivalent_ids = equivalence[inputs[case_id]["assembly_sequence_sha256"]]
        equivalent_types = sorted({truth[other]["event_type"] for other in equivalent_ids})
        identifiable_type = len(equivalent_types) == 1
        row = dict(case_id=case_id, status=status, truth_event_status=target["event_status"],
                   truth_event_type=target["event_type"], event_score=value,
                   observable_equivalence_class_sha256=inputs[case_id]["assembly_sequence_sha256"],
                   observable_equivalence_class_size=len(equivalent_ids),
                   observationally_equivalent_case_ids=equivalent_ids,
                   observationally_equivalent_event_types=equivalent_types,
                   event_type_identifiability_status="identifiable_within_bundle" if identifiable_type else "not_identifiable_from_available_sequences",
                   decision=None, confusion_cell=None, event_type_match=None,
                   edited_label_path_distance=None, orientation_path_accuracy=None,
                   signed_bp_delta_abs_error=None, edited_array_interval_iou=None)
        if status == "ok":
            called = value >= 0.5
            row["decision"] = "positive" if called else "intact_negative"
            row["confusion_cell"] = ("tp" if positive else "fp") if called else ("fn" if positive else "tn")
            if positive and identifiable_type and pred.get("event_type") is not None:
                row["event_type_match"] = pred["event_type"] == target["event_type"]
            labels = pred.get("predicted_edited_label_path")
            if isinstance(labels, list) and all(isinstance(x, str) for x in labels):
                row["edited_label_path_distance"] = edit_distance(target["edited_label_path"], labels)
                orientations = pred.get("predicted_edited_orientation_path")
                if (isinstance(orientations, list) and len(orientations) == len(labels)
                        == len(target["edited_orientation_path"])
                        and all(x in {"+", "-"} for x in orientations)):
                    row["orientation_path_accuracy"] = sum(
                        a == b for a, b in zip(orientations, target["edited_orientation_path"])) / len(orientations) if orientations else 1.0
            bp = pred.get("predicted_signed_bp_delta")
            if finite_number(bp):
                row["signed_bp_delta_abs_error"] = abs(bp - target["signed_bp_delta"])
            interval = pred.get("predicted_edited_interval_bp")
            if interval is not None:
                row["edited_array_interval_iou"] = interval_iou(target["edited_array_interval_bp"], interval)
        per_case.append(row)
    all_ok = all(row["status"] == "ok" for row in per_case)
    counts = {key: sum(row["confusion_cell"] == key for row in per_case) for key in ("tp", "fn", "fp", "tn")}
    positive_count = sum(row["truth_event_status"] == "positive" for row in per_case)
    negative_count = len(per_case) - positive_count
    primary = {key: counts[key] if all_ok else None for key in counts}
    primary.update(sensitivity=counts["tp"] / positive_count if all_ok else None,
                   fpr=counts["fp"] / negative_count if all_ok else None,
                   ppv=counts["tp"] / (counts["tp"] + counts["fp"]) if all_ok and counts["tp"] + counts["fp"] else None)
    type_values = [row["event_type_match"] for row in per_case if row["truth_event_status"] == "positive"]
    secondary = dict(
        event_type_accuracy=mean_or_block(type_values, positive_count) if all_ok else None,
        edited_label_path_mean_distance=mean_or_block([row["edited_label_path_distance"] for row in per_case], len(per_case)) if all_ok else None,
        edited_orientation_path_mean_accuracy=mean_or_block([row["orientation_path_accuracy"] for row in per_case], len(per_case)) if all_ok else None,
        signed_bp_delta_mae=mean_or_block([row["signed_bp_delta_abs_error"] for row in per_case], len(per_case)) if all_ok else None,
        edited_array_interval_mean_iou=mean_or_block([row["edited_array_interval_iou"] for row in per_case], len(per_case)) if all_ok else None)
    summary = dict(schema_version=2, benchmark="B1_structure_development_v2",
                   split="development", material_count=receipt["material_count"],
                   denominator=len(per_case), positive_count=positive_count, intact_negative_count=negative_count,
                   observable_equivalence_class_count=len(equivalence),
                   event_type_nonidentifiable_case_count=sum(row["event_type_identifiability_status"] == "not_identifiable_from_available_sequences" for row in per_case),
                   status_counts={status: sum(row["status"] == status for row in per_case)
                                  for status in ("ok", "abstain", "unmapped", "failed", "not_run", "not_reported")},
                   primary_status="complete" if all_ok else "blocked_incomplete_predictions",
                   primary=primary, secondary=secondary,
                   secondary_status="complete" if all_ok and all(value is not None for value in secondary.values()) else "blocked_incomplete_predictions_or_fields",
                   auprc_status=metrics["auprc_status"],
                   biological_accuracy_status=metrics["physical_or_read_assembly_accuracy_status"],
                   confidence_interval_status=metrics["confidence_interval_status"],
                   prediction_sha256=sha256_file(predictions_path), metrics_sha256=sha256_file(metrics_path),
                   inputs_sha256=receipt["inputs_sha256"], truth_sha256=receipt["truth_sha256"])
    outdir.mkdir(parents=True, exist_ok=False)
    write_jsonl(outdir / "per_case.jsonl", per_case)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    score(args.bundle, args.metrics, args.predictions, args.outdir)


if __name__ == "__main__":
    main()
