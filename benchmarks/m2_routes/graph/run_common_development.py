"""Run graph route on frozen B1 public sequence inputs without loading truth."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from time import perf_counter
import tracemalloc

from . import prototype
from .prototype import ReadPath, audit_transition_graph, to_common_prediction


_RC = str.maketrans("ACGT", "TGCA")


class ComponentUnresolved(ValueError):
    pass


def _extract_array(sequence: str, left: str, right: str) -> str:
    if (not left or not right or sequence.count(left) != 1 or sequence.count(right) != 1
            or not sequence.startswith(left) or not sequence.endswith(right)
            or len(sequence) < len(left) + len(right)):
        raise ComponentUnresolved("engineered_unique_flank_check_failed")
    return sequence[len(left):len(sequence) - len(right)]


def _exact_tile_labels(array: str, monomers: dict[str, str]) -> tuple[str, ...]:
    lengths = {len(sequence) for sequence in monomers.values()}
    if len(lengths) != 1:
        raise ComponentUnresolved("unequal_candidate_tile_lengths")
    tile_length = lengths.pop()
    if not tile_length or len(array) % tile_length:
        raise ComponentUnresolved("partial_or_nontiled_array_interval")
    index: dict[str, list[str]] = {}
    for label, sequence in monomers.items():
        index.setdefault(sequence, []).append(label)
    labels = []
    for start in range(0, len(array), tile_length):
        matches = index.get(array[start:start + tile_length], ())
        if len(matches) != 1:
            raise ComponentUnresolved("tile_not_uniquely_exact_in_catalogue")
        labels.append(matches[0])
    return tuple(labels)


def _abstain(case_id: str, reason: str) -> dict[str, object]:
    return {"case_id": case_id, "status": "abstain", "event_score": None,
            "event_type": None, "predicted_edited_label_path": None,
            "predicted_edited_interval_bp": None, "predicted_signed_bp_delta": None,
            "predicted_edited_orientation_path": None,
            "predicted_copy_spans_bp": None, "reason": reason,
            "audit_state": "COMPONENT_UNRESOLVED"}


def predict_case(row: dict[str, object]) -> tuple[dict[str, object], dict[str, object]]:
    case_id = row["case_id"]
    if (row.get("schema_version") != 2 or row.get("split") != "development"
            or row.get("input_status") != "ok"
            or row.get("read_pairing_status") != "synthetic_simulated"):
        return _abstain(case_id, "unsupported_or_unverified_input_contract"), {}
    error = row["error_profile"]
    if (error.get("name") != "exact_synthetic" or error.get("indel_rate") != 0
            or error.get("substitution_rate") != 0):
        return _abstain(case_id, "nonexact_read_error_profile_not_supported"), {}
    monomers: dict[str, str] = {}
    for label, sequence in row["candidate_monomers"].items():
        monomers[label + "+"] = sequence
        monomers[label + "-"] = sequence.translate(_RC)[::-1]
    try:
        assembly_array = _extract_array(row["assembly_sequence"],
                                        row["left_flank_sequence"], row["right_flank_sequence"])
        assembly = _exact_tile_labels(assembly_array, monomers)
        reads = []
        for molecule_id, sequence in sorted(row["raw_read_sequences"].items()):
            array = _extract_array(sequence, row["left_flank_sequence"],
                                   row["right_flank_sequence"])
            labels = _exact_tile_labels(array, monomers)
            reads.append(ReadPath(molecule_id, labels, (1.0,) * len(labels),
                                  "exact_synthetic", both_flanks_verified=True,
                                  haplotype=row.get("haplotype_id")))
    except ComponentUnresolved as exc:
        return _abstain(case_id, str(exc)), {}
    audit = audit_transition_graph(assembly, reads, monomer_sequences=monomers,
                                   min_error_profiles=1,
                                   synthetic_error_free_mode=True)
    return to_common_prediction(case_id, audit), {
        "assembly_copy_count": len(assembly),
        "read_copy_counts": [len(read.labels) for read in reads],
        "graph_audit": asdict(audit),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()
    input_bytes = args.inputs.read_bytes()
    rows = [json.loads(line) for line in input_bytes.splitlines() if line.strip()]
    if len({row["case_id"] for row in rows}) != len(rows):
        raise ValueError("Duplicate opaque case IDs in input")
    predictions = []
    audits = []
    for row in rows:
        tracemalloc.start()
        start = perf_counter()
        prediction, detail = predict_case(row)
        elapsed = perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        predictions.append(prediction)
        audits.append({"case_id": row["case_id"], "prediction_status": prediction["status"],
                       "audit_state": prediction["audit_state"],
                       "elapsed_seconds": elapsed,
                       "tracemalloc_peak_bytes": peak, **detail})
    args.predictions.parent.mkdir(parents=True, exist_ok=True)
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    prediction_text = "".join(json.dumps(row, sort_keys=True) + "\n"
                                   for row in predictions)
    args.predictions.write_text(prediction_text)
    receipt = {
        "status": "development_input_only_exact_synthetic_component",
        "inputs_sha256": hashlib.sha256(input_bytes).hexdigest(),
        "predictions_sha256": hashlib.sha256(prediction_text.encode()).hexdigest(),
        "prototype_sha256": hashlib.sha256(Path(prototype.__file__).read_bytes()).hexdigest(),
        "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "case_count": len(rows),
        "prediction_status_counts": dict(Counter(row["status"] for row in predictions)),
        "audit_state_counts": dict(Counter(row["audit_state"] for row in predictions)),
        "truth_file_opened": False,
        "limitations": [
            "24-bp equal-length exact tile decomposition is possible only by construction",
            "one named error profile is allowed only because all input reads are exact synthetic",
            "synthetic flanks do not prove biological same-locus anchoring",
            "no noisy reads, held-out donor, or real molecule validation",
            "event scores are binary decisions, not calibrated probabilities",
        ],
        "cases": audits,
    }
    args.audit.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt["audit_state_counts"], sort_keys=True))


if __name__ == "__main__":
    main()
