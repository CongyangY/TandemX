"""Score frozen Col-CEN v3 injected edits with explicit abstention costs and N/A gates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.scripts.build_b1_structure_cases import sha256_file, write_jsonl
from benchmarks.scripts.score_b1_structure_cases import finite_number, keyed, read_jsonl


def edited_breakpoints(truth: dict) -> list[int]:
    points = set()
    for segment in truth["coordinate_chain"]:
        operation = segment["operation"]
        if operation == "delete":
            points.add(segment["edited_start0"])
        elif operation in {"duplicate", "reorder", "reverse_complement"}:
            points.update((segment["edited_start0"], segment["edited_end0"]))
    return sorted(points)


def score(bundle: Path, metrics_path: Path, predictions_path: Path, outdir: Path,
          resource_receipt_path: Path | None = None) -> dict:
    receipt = json.loads((bundle / "receipt.json").read_text())
    metrics = json.loads(metrics_path.read_text())
    if metrics["case_count"] != receipt["case_count"] or metrics["schema_version"] != 2:
        raise ValueError("v3 frozen metrics/bundle mismatch")
    for field, filename in (("inputs_sha256", "inputs.jsonl"), ("truth_sha256", "truth.jsonl")):
        if receipt[field] != sha256_file(bundle / filename):
            raise ValueError(f"bundle {filename} hash mismatch")
    inputs = keyed(read_jsonl(bundle / "inputs.jsonl"), "inputs")
    truths = keyed(read_jsonl(bundle / "truth.jsonl"), "truth")
    predictions = keyed(read_jsonl(predictions_path), "predictions")
    ids = receipt["case_ids"]
    if len(ids) != metrics["case_count"] or set(ids) != set(inputs) or set(ids) != set(truths) or set(predictions) - set(ids):
        raise ValueError("case ID mismatch")
    equivalence = {}
    for case_id in ids:
        if truths[case_id]["edited_sequence_sha256"] != inputs[case_id]["assembly_sequence_sha256"]:
            raise ValueError("truth/input edited sequence hash mismatch")
        equivalence.setdefault(inputs[case_id]["assembly_sequence_sha256"], []).append(case_id)
    rows = []
    for case_id in ids:
        truth = truths[case_id]
        pred = predictions.get(case_id)
        status = "not_reported" if pred is None else pred.get("status")
        if status not in {"ok", "abstain", "unmapped", "failed", "not_run", "not_reported"}:
            raise ValueError(f"invalid prediction status: {case_id}")
        event_score = None if pred is None else pred.get("event_score")
        if status == "ok":
            if not finite_number(event_score) or not 0 <= event_score <= 1:
                raise ValueError(f"ok requires finite event_score in [0,1]: {case_id}")
        elif event_score is not None:
            raise ValueError(f"non-ok score must be null: {case_id}")
        positive = truth["event_status"] == "positive"
        similar = equivalence[inputs[case_id]["assembly_sequence_sha256"]]
        identifiable = len({truths[other]["event_type"] for other in similar}) == 1
        cell = ("fn" if positive else "unresolved_negative") if status != "ok" else (
            ("tp" if positive else "fp") if event_score >= metrics["decision_threshold"]
            else ("fn" if positive else "tn"))
        row = dict(case_id=case_id, array_id=truth["array_id"], truth_event_type=truth["event_type"],
                   truth_event_status=truth["event_status"], status=status,
                   event_score=event_score, intent_to_diagnose_cell=cell,
                   observable_equivalence_class_sha256=inputs[case_id]["assembly_sequence_sha256"],
                   observable_equivalence_class_size=len(similar),
                   event_type_identifiable_within_bundle=identifiable,
                   event_type_match=None, signed_bp_delta_abs_error=None,
                   truth_breakpoints_edited_bp=edited_breakpoints(truth),
                   breakpoint_localization_error_bp=None, breakpoint_status="not_reported")
        if status == "ok":
            if positive and identifiable and pred.get("event_type") is not None:
                row["event_type_match"] = pred["event_type"] == truth["event_type"]
            delta = pred.get("predicted_signed_bp_delta")
            if finite_number(delta):
                row["signed_bp_delta_abs_error"] = abs(delta - truth["signed_bp_delta"])
            points = pred.get("predicted_breakpoints_bp")
            expected = row["truth_breakpoints_edited_bp"]
            if not positive:
                row["breakpoint_status"] = "not_applicable_intact"
            elif not identifiable:
                row["breakpoint_status"] = "not_identifiable_exact_sequence_equivalence"
            elif points is not None:
                if (not isinstance(points, list) or len(points) != len(expected)
                        or any(not isinstance(x, int) or isinstance(x, bool) or x < 0 for x in points)):
                    row["breakpoint_status"] = "invalid_or_count_mismatch"
                elif len(expected) == 0:
                    row["breakpoint_status"] = "no_defined_junction"
                else:
                    row["breakpoint_localization_error_bp"] = sum(abs(a - b) for a, b in zip(sorted(points), expected)) / len(expected)
                    row["breakpoint_status"] = "scored"
        rows.append(row)
    counts = {cell: sum(row["intent_to_diagnose_cell"] == cell for row in rows)
              for cell in ("tp", "fn", "fp", "tn", "unresolved_negative")}
    positives = metrics["positive_count"]
    negatives = metrics["intact_negative_count"]
    if counts["tp"] + counts["fn"] != positives or counts["fp"] + counts["tn"] + counts["unresolved_negative"] != negatives:
        raise ValueError("intent-to-diagnose denominator conservation failed")
    ok = [row for row in rows if row["status"] == "ok"]
    ok_pos = sum(row["truth_event_status"] == "positive" for row in ok)
    ok_neg = len(ok) - ok_pos
    eligible_counts = {cell: sum(row["intent_to_diagnose_cell"] == cell for row in ok)
                       for cell in ("tp", "fn", "fp", "tn")}
    positive_type = [row for row in rows if row["truth_event_status"] == "positive" and row["event_type_identifiable_within_bundle"]]
    bp_rows = [row for row in rows if row["truth_event_status"] == "positive" and row["event_type_identifiable_within_bundle"]
               and row["truth_breakpoints_edited_bp"]]
    def complete_mean(selected: list[dict], field: str) -> float | None:
        values = [row[field] for row in selected]
        return sum(values) / len(values) if values and all(value is not None for value in values) else None
    resource = json.loads(resource_receipt_path.read_text()) if resource_receipt_path else {}
    runtime = resource.get("wall_elapsed_seconds")
    rss = resource.get("peak_rss_bytes")
    if runtime is not None and (not finite_number(runtime) or runtime < 0):
        raise ValueError("invalid measured wall elapsed seconds")
    if rss is not None and (not isinstance(rss, int) or rss < 0):
        raise ValueError("invalid measured peak RSS bytes")
    summary = dict(schema_version=2, benchmark="B1_structure_development_v3", split="development",
                   denominator=len(rows), source_array_count=receipt["source_array_count"],
                   independent_donor_count=receipt["independent_donor_count"],
                   positive_count=positives, intact_negative_count=negatives,
                   observable_equivalence_class_count=len(equivalence),
                   status_counts={status: sum(row["status"] == status for row in rows)
                                  for status in ("ok", "abstain", "unmapped", "failed", "not_run", "not_reported")},
                   primary_intent_to_diagnose=dict(**counts,
                       sensitivity=counts["tp"] / positives,
                       fpr=counts["fp"] / negatives,
                       negative_failure_rate=(counts["fp"] + counts["unresolved_negative"]) / negatives),
                   eligible_only=dict(denominator=len(ok), positive_count=ok_pos, intact_negative_count=ok_neg,
                                      coverage=len(ok) / len(rows), **eligible_counts,
                                      sensitivity=eligible_counts["tp"] / ok_pos if ok_pos else None,
                                      fpr=eligible_counts["fp"] / ok_neg if ok_neg else None),
                   ambiguity_rejected_fraction=(len(rows) - len(ok)) / len(rows),
                   event_type_accuracy=complete_mean(positive_type, "event_type_match"),
                   event_type_denominator=len(positive_type),
                   signed_bp_delta_mae=complete_mean(ok, "signed_bp_delta_abs_error") if len(ok) == len(rows) else None,
                   signed_bp_delta_status="complete" if len(ok) == len(rows) and all(row["signed_bp_delta_abs_error"] is not None for row in ok) else "NA_missing_or_non_ok_predictions",
                   breakpoint_localization_error_bp=complete_mean(bp_rows, "breakpoint_localization_error_bp"),
                   breakpoint_denominator=len(bp_rows),
                   breakpoint_status="complete" if bp_rows and all(row["breakpoint_localization_error_bp"] is not None for row in bp_rows) else "NA_missing_or_nonidentifiable_breakpoint_predictions",
                   hor_period_error_bp=None, hor_period_status="NA_unverified_native_HOR_period",
                   repeat_order_accuracy=None, repeat_order_status="NA_native_tile_identity_and_boundary_not_independently_validated",
                   runtime_seconds=runtime, runtime_status="measured" if runtime is not None else "NA_no_route_wall_time_receipt",
                   peak_rss_bytes=rss, peak_rss_status="measured" if rss is not None else "NA_no_route_peak_RSS_receipt",
                   auprc_status=metrics["auprc_status"],
                   confidence_interval_status=metrics["confidence_interval_status"],
                   biological_accuracy_status=metrics["physical_or_read_assembly_accuracy_status"],
                   prediction_sha256=sha256_file(predictions_path), metrics_sha256=sha256_file(metrics_path),
                   inputs_sha256=receipt["inputs_sha256"], truth_sha256=receipt["truth_sha256"])
    outdir.mkdir(parents=True, exist_ok=False)
    write_jsonl(outdir / "per_case.jsonl", rows)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--resource-receipt", type=Path)
    args = parser.parse_args()
    score(args.bundle, args.metrics, args.predictions, args.outdir, args.resource_receipt)


if __name__ == "__main__":
    main()
