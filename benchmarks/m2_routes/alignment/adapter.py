"""Adapter from frozen B1 development inputs to research route A predictions."""

from __future__ import annotations

import argparse
import json
import platform
import resource
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

from .prototype import AuditResult, audit


def _between_unique_flanks(sequence: str, left: str, right: str) -> tuple[str, tuple[int, int]]:
    if not left or not right or sequence.count(left) != 1 or sequence.count(right) != 1:
        raise ValueError("flank_absent_or_nonunique")
    start = sequence.index(left) + len(left)
    end = sequence.index(right)
    if end < start:
        raise ValueError("flank_order_or_overlap_invalid")
    return sequence[start:end], (start, end)


def predict_case(row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one B1 prediction without reading or requiring truth."""
    case_id = str(row["case_id"])
    prediction: dict[str, Any] = {
        "case_id": case_id,
        "status": "abstain",
        "event_score": None,
        "event_type": None,
        "predicted_edited_label_path": None,
        "predicted_edited_orientation_path": None,
        "predicted_copy_spans_bp": None,
        "predicted_edited_interval_bp": None,
        "predicted_signed_bp_delta": None,
        "reason": "",
        "audit_state": "AMBIGUOUS",
        "route": "local_edlib_monomer_path_dp_v1",
        "score_interpretation": "uncalibrated_binary_decision",
    }
    if row.get("input_status") != "ok":
        prediction["reason"] = "input_status_not_ok"
        return prediction
    if row.get("array_window_policy") != "between_exact_unique_synthetic_flanks":
        prediction["reason"] = "unsupported_array_window_policy"
        return prediction
    try:
        left = str(row["left_flank_sequence"])
        right = str(row["right_flank_sequence"])
        assembly, interval = _between_unique_flanks(str(row["assembly_sequence"]), left, right)
        raw_reads = row["raw_read_sequences"]
        if not isinstance(raw_reads, dict):
            raise ValueError("raw_reads_not_mapping")
        reads = {str(key): _between_unique_flanks(str(seq), left, right)[0]
                 for key, seq in raw_reads.items()}
        monomers = row["candidate_monomers"]
        if not isinstance(monomers, dict):
            raise ValueError("monomers_not_mapping")
        result: AuditResult = audit(assembly, reads, monomers,
                                    pairing_status=str(row["read_pairing_status"]))
    except (KeyError, TypeError, ValueError) as exc:
        prediction["reason"] = "invalid_or_untrimmed_input:" + str(exc)
        return prediction
    prediction["audit_state"] = result.state
    prediction["reason"] = result.reason
    prediction["support_count"] = result.support_count
    prediction["discordant_count"] = result.discordant_count
    prediction["ambiguous_count"] = result.ambiguous_count
    prediction["label_edit_distance"] = result.label_edit_distance
    prediction["breakpoint_label_interval"] = result.breakpoint_label_interval
    if result.assembly.state == "RESOLVED":
        prediction["predicted_edited_label_path"] = [c.label for c in result.assembly.copies]
        prediction["predicted_edited_orientation_path"] = [c.orientation for c in result.assembly.copies]
        prediction["predicted_copy_spans_bp"] = [
            [interval[0] + c.start, interval[0] + c.end] for c in result.assembly.copies
        ]
        prediction["predicted_edited_interval_bp"] = list(interval)
    if result.state == "DISCORDANT":
        prediction["status"] = "ok"
        prediction["event_score"] = 1.0
        prediction["event_type"] = result.candidate_event
    elif result.state == "SUPPORTED":
        prediction["status"] = "ok"
        prediction["event_score"] = 0.0
        prediction["event_type"] = "intact"
    elif result.state == "INSUFFICIENT_READ_SUPPORT":
        prediction["reason"] = "insufficient_read_support:" + result.reason
    return prediction


def run_cases(input_path: Path, output_path: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    with input_path.open() as source:
        for line in source:
            if not line.strip():
                continue
            row = json.loads(line)
            start = perf_counter()
            prediction = predict_case(row)
            prediction["elapsed_seconds"] = perf_counter() - start
            rows.append(prediction)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as target:
        for row in rows:
            target.write(json.dumps(row, sort_keys=True) + "\n")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {"cases": len(rows), "ok": sum(row["status"] == "ok" for row in rows),
            "abstain": sum(row["status"] == "abstain" for row in rows),
            "seconds": sum(row["elapsed_seconds"] for row in rows),
            "peak_process_rss_bytes": peak if platform.system() == "Darwin" else peak * 1024}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_cases(args.inputs, args.predictions), sort_keys=True))


if __name__ == "__main__":
    main()
