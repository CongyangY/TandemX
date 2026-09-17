"""Score frozen synthetic B2 structure hold-out cases without biological claims."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from benchmarks.scripts.build_b1_structure_cases import sha256_file, write_jsonl
from benchmarks.scripts.score_b1_structure_cases import edit_distance, finite_number, keyed, read_jsonl
from benchmarks.scripts.score_b1_colcen_cases import edited_breakpoints


STATUSES = ("ok", "abstain", "unmapped", "failed", "not_run", "not_reported")
CELLS = ("tp", "fn", "fp", "tn", "unresolved_negative")


def counts(rows: list[dict]) -> dict[str, int]:
    return {cell: sum(row["intent_to_diagnose_cell"] == cell for row in rows) for cell in CELLS}


def complete_mean(rows: list[dict], field: str) -> float | None:
    values = [row[field] for row in rows]
    return sum(values) / len(values) if values and all(value is not None for value in values) else None


def score(bundle: Path, metrics_path: Path, predictions_path: Path, outdir: Path,
          resource_receipt_path: Path | None = None) -> dict:
    receipt = json.loads((bundle / "receipt.json").read_text())
    metrics = json.loads(metrics_path.read_text())
    if (receipt["split"] != metrics["split"] or metrics["split"] != "synthetic_held_out"
            or receipt["case_count"] != metrics["case_count"] or metrics["schema_version"] != 2):
        raise ValueError("B2 frozen metrics/bundle mismatch")
    if sha256_file(metrics_path.parent / "protocol.json") != receipt["protocol_sha256"]:
        raise ValueError("B2 protocol hash mismatch")
    for field, filename in (("inputs_sha256", "inputs.jsonl"), ("truth_sha256", "truth.jsonl")):
        if sha256_file(bundle / filename) != receipt[field]:
            raise ValueError(f"B2 {filename} hash mismatch")
    inputs = keyed(read_jsonl(bundle / "inputs.jsonl"), "inputs")
    truths = keyed(read_jsonl(bundle / "truth.jsonl"), "truth")
    predictions = keyed(read_jsonl(predictions_path), "predictions")
    ids = receipt["case_ids"]
    if len(ids) != metrics["case_count"] or set(ids) != set(inputs) or set(ids) != set(truths) or set(predictions) - set(ids):
        raise ValueError("B2 case IDs mismatch")
    equivalents: dict[str, list[str]] = {}
    for case_id in ids:
        if truths[case_id]["edited_sequence_sha256"] != inputs[case_id]["assembly_sequence_sha256"]:
            raise ValueError("B2 truth/input edited sequence mismatch")
        equivalents.setdefault(inputs[case_id]["assembly_sequence_sha256"], []).append(case_id)
    rows = []
    for case_id in ids:
        truth, public = truths[case_id], inputs[case_id]
        pred = predictions.get(case_id)
        status = "not_reported" if pred is None else pred.get("status")
        if status not in STATUSES:
            raise ValueError(f"invalid B2 status: {case_id}")
        event_score = None if pred is None else pred.get("event_score")
        if status == "ok":
            if not finite_number(event_score) or not 0 <= event_score <= 1:
                raise ValueError(f"ok requires finite event_score in [0,1]: {case_id}")
        elif event_score is not None:
            raise ValueError(f"non-ok score must be null: {case_id}")
        positive = truth["event_status"] == "positive"
        equivalent = equivalents[public["assembly_sequence_sha256"]]
        identifiable = len({truths[other]["event_type"] for other in equivalent}) == 1
        cell = ("fn" if positive else "unresolved_negative") if status != "ok" else (
            ("tp" if positive else "fp") if event_score >= metrics["decision_threshold"]
            else ("fn" if positive else "tn"))
        row = dict(case_id=case_id, material_id=public["material_id"],
                   monomer_length_bp=public["monomer_length_bp"],
                   synthetic_read_support=public["synthetic_read_support"],
                   error_mode="noisy" if public["error_profile"]["substitution_rate"] > 0 else "exact",
                   status=status, event_score=event_score, truth_event_status=truth["event_status"],
                   truth_event_type=truth["event_type"], intent_to_diagnose_cell=cell,
                   observable_equivalence_class_size=len(equivalent),
                   event_type_identifiable_within_bundle=identifiable,
                   event_type_match=None, signed_bp_delta_abs_error=None,
                   label_path_edit_distance=None, repeat_order_accuracy=None,
                   HOR_period_abs_error_monomers=None,
                   truth_breakpoints_edited_bp=edited_breakpoints(truth),
                   breakpoint_localization_error_bp=None,
                   breakpoint_status="not_reported")
        if status == "ok":
            if positive and identifiable and pred.get("event_type") is not None:
                row["event_type_match"] = pred["event_type"] == truth["event_type"]
            delta = pred.get("predicted_signed_bp_delta")
            if finite_number(delta):
                row["signed_bp_delta_abs_error"] = abs(delta - truth["signed_bp_delta"])
            labels = pred.get("predicted_edited_label_path")
            if isinstance(labels, list) and all(isinstance(x, str) for x in labels):
                target = truth["edited_label_path"]
                row["label_path_edit_distance"] = edit_distance(labels, target)
                if len(labels) == len(target):
                    row["repeat_order_accuracy"] = (sum(a == b for a, b in zip(labels, target)) / len(target)
                                                    if target else 1.0)
            period = pred.get("predicted_HOR_period_monomers")
            if isinstance(period, int) and not isinstance(period, bool) and period > 0:
                row["HOR_period_abs_error_monomers"] = abs(period - truth["source_HOR_period_monomers"])
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
                elif not expected:
                    row["breakpoint_status"] = "no_defined_junction"
                else:
                    row["breakpoint_localization_error_bp"] = sum(abs(a - b) for a, b in zip(sorted(points), expected)) / len(expected)
                    row["breakpoint_status"] = "scored"
        rows.append(row)
    overall = counts(rows)
    positives, negatives = metrics["positive_count"], metrics["intact_negative_count"]
    if (overall["tp"] + overall["fn"] != positives or
            overall["fp"] + overall["tn"] + overall["unresolved_negative"] != negatives):
        raise ValueError("B2 primary denominator conservation failed")
    ok = [row for row in rows if row["status"] == "ok"]
    ok_counts = counts(ok)
    ok_pos = sum(row["truth_event_status"] == "positive" for row in ok)
    ok_neg = len(ok) - ok_pos
    positive_identifiable = [row for row in rows if row["truth_event_status"] == "positive"
                             and row["event_type_identifiable_within_bundle"]]
    breakpoint_applicable = [row for row in positive_identifiable if row["truth_breakpoints_edited_bp"]]
    resource = json.loads(resource_receipt_path.read_text()) if resource_receipt_path else {}
    runtime = resource.get("wall_elapsed_seconds")
    rss = resource.get("peak_rss_bytes")
    if runtime is not None and (not finite_number(runtime) or runtime < 0):
        raise ValueError("invalid B2 runtime receipt")
    if rss is not None and (not isinstance(rss, int) or isinstance(rss, bool) or rss < 0):
        raise ValueError("invalid B2 RSS receipt")
    strata = {}
    for field in ("monomer_length_bp", "synthetic_read_support", "error_mode"):
        strata[field] = {}
        for key in sorted({row[field] for row in rows}):
            selected = [row for row in rows if row[field] == key]
            strata[field][str(key)] = dict(denominator=len(selected), **counts(selected),
                                           non_ok=sum(row["status"] != "ok" for row in selected))
    summary = dict(schema_version=2, benchmark="B2_synthetic_structural_held_out_v1",
                   split="synthetic_held_out", denominator=len(rows), positive_count=positives,
                   intact_negative_count=negatives,
                   synthetic_founder_count=receipt["synthetic_founder_count"], biological_donor_count=0,
                   observable_equivalence_class_count=len(equivalents),
                   status_counts={status: sum(row["status"] == status for row in rows) for status in STATUSES},
                   primary_intent_to_diagnose=dict(**overall,
                       sensitivity=overall["tp"] / positives, fpr=overall["fp"] / negatives,
                       negative_failure_rate=(overall["fp"] + overall["unresolved_negative"]) / negatives),
                   eligible_only=dict(denominator=len(ok), coverage=len(ok) / len(rows),
                                      positive_count=ok_pos, intact_negative_count=ok_neg,
                                      **{cell: ok_counts[cell] for cell in ("tp", "fn", "fp", "tn")},
                                      sensitivity=ok_counts["tp"] / ok_pos if ok_pos else None,
                                      fpr=ok_counts["fp"] / ok_neg if ok_neg else None),
                   ambiguity_rejected_fraction=(len(rows) - len(ok)) / len(rows),
                   strata=strata,
                   event_type_accuracy=complete_mean(positive_identifiable, "event_type_match"),
                   event_type_denominator=len(positive_identifiable),
                   signed_bp_delta_mae=complete_mean(rows, "signed_bp_delta_abs_error"),
                   signed_bp_delta_status="complete" if all(row["signed_bp_delta_abs_error"] is not None for row in rows) else "NA_missing_or_non_ok_predictions",
                   breakpoint_localization_error_bp=complete_mean(breakpoint_applicable, "breakpoint_localization_error_bp"),
                   breakpoint_denominator=len(breakpoint_applicable),
                   breakpoint_status="complete" if breakpoint_applicable and all(row["breakpoint_localization_error_bp"] is not None for row in breakpoint_applicable) else "NA_missing_or_nonidentifiable_breakpoint_predictions",
                   HOR_period_error_monomers=complete_mean(rows, "HOR_period_abs_error_monomers"),
                   HOR_period_status="complete" if all(row["HOR_period_abs_error_monomers"] is not None for row in rows) else "NA_missing_or_non_ok_HOR_predictions",
                   repeat_order_accuracy=complete_mean(rows, "repeat_order_accuracy"),
                   repeat_order_status="complete" if all(row["repeat_order_accuracy"] is not None for row in rows) else "NA_missing_or_length_mismatched_label_paths",
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
