"""Independent scoring against planted genome copy counts and interval unions."""
from __future__ import annotations

import math
from pathlib import Path

from benchmarks.challenge.schema import read_table


def truth_by_family(truth: list[dict]) -> dict[str, dict]:
    grouped: dict[str, dict] = {}
    for row in truth:
        family = row["family_id"]
        period = int(row["period"])
        start, end = int(row["start"]), int(row["end"])
        if period <= 0 or start > end or int(row["repeat_bp"]) != end-start:
            raise ValueError("Invalid truth interval")
        record = grouped.setdefault(family, dict(
            family_id=family, period=period, copies=0, repeat_bp=0, intervals=[]
        ))
        if record["period"] != period:
            raise ValueError("One family cannot have multiple truth periods")
        record["copies"] += int(row["copies"])
        record["repeat_bp"] += int(row["repeat_bp"])
        if end > start:
            record["intervals"].append((start, end))
    if not grouped:
        raise ValueError("Truth must contain at least one family")
    for record in grouped.values():
        record["intervals"] = union(record["intervals"])
        if sum(end-start for start, end in record["intervals"]) != record["repeat_bp"]:
            raise ValueError("Overlapping truth intervals within a family")
    return grouped


def finite_nonnegative(value: str, field: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"Invalid nonnegative {field}: {value}")
    return result


def score_copy_number(path: Path, truth: list[dict], sampling: dict) -> list[dict]:
    observed = read_table(path, {"family_id", "estimated_copy_number", "copy_number_interval_low", "copy_number_interval_high"})
    by_family = {r["family_id"]: r for r in observed}
    expected = truth_by_family(truth)
    if len(by_family) != len(observed) or set(by_family) != set(expected):
        raise ValueError("Copy-number families differ from supplied known catalogue")
    scores = []
    for family_id, row in expected.items():
        value = by_family[family_id]
        est, lo, hi = [finite_nonnegative(value[f], f) for f in
                       ("estimated_copy_number", "copy_number_interval_low", "copy_number_interval_high")]
        if lo > hi:
            raise ValueError("Reversed copy-number interval")
        copies = int(row["copies"])
        oracle = sampling["sampled_repeat_bp"][family_id] / int(row["period"]) / sampling["actual_base_coverage"]
        scores.append(dict(family_id=family_id, truth_copies=copies, estimate=est,
                           signed_relative_error=(est-copies)/copies, absolute_relative_error=abs(est-copies)/copies,
                           interval_low=lo, interval_high=hi, interval_contains_truth=lo<=copies<=hi,
                           interval_relative_width=(hi-lo)/copies, sampling_oracle_copy_estimate=oracle,
                           oracle_relative_error=(oracle-copies)/copies,
                           estimator_minus_sampling_oracle=(est-oracle)/copies,
                           interval_kind="diagnostic_kmer_10_90_spread_not_sampling_CI"))
    return scores


def union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def score_localization(path: Path, truth: list[dict], genome_bp: int) -> list[dict]:
    expected = truth_by_family(truth)
    predictions = {family_id: [] for family_id in expected}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 4 or fields[0] != "chr_sim" or fields[3] not in predictions:
            raise ValueError("Unknown or malformed assembly prediction")
        start, end = map(int, fields[1:3])
        if not 0 <= start < end <= genome_bp:
            raise ValueError("Invalid assembly prediction coordinates")
        predictions[fields[3]].append((start, end))
    rows = []
    for family_id, record in expected.items():
        spans = union(predictions[family_id])
        bp = sum(b-a for a,b in spans)
        overlap = sum(max(0, min(pred_end, truth_end)-max(pred_start, truth_start))
                      for pred_start, pred_end in spans
                      for truth_start, truth_end in record["intervals"])
        true_bp = int(record["repeat_bp"])
        rows.append(dict(family_id=family_id, true_assembly_bp=true_bp, predicted_assembly_bp=bp,
                         overlap_bp=overlap, base_recall=overlap/true_bp if true_bp else None,
                         base_precision=overlap/bp if bp else None, fragments=len(spans),
                         truth_fragments=len(record["intervals"])))
    return rows


def score_comparison(path: Path, truth: list[dict], assembly_truth: list[dict], threshold: float=.6) -> list[dict]:
    observed = read_table(path, {"family_id", "assembly_read_ratio", "status"})
    by_family = {r["family_id"]: r for r in observed}
    full = truth_by_family(truth)
    assembly = truth_by_family(assembly_truth)
    if len(by_family) != len(observed) or set(by_family) != set(full):
        raise ValueError("Comparison families differ from known catalogue")
    if set(assembly) != set(full):
        raise ValueError("Assembly truth families differ from full truth")
    result = []
    for family, row in assembly.items():
        ratio = int(row["repeat_bp"]) / int(full[family]["repeat_bp"])
        status = by_family[family]["status"]
        # 'reads_only' is a separate observed absence category. Report an explicit
        # under-representation endpoint and retain native status for all rows.
        call = status in {"possible_collapse", "reads_only"}
        positive = ratio < threshold
        result.append(dict(family_id=family, truth_assembly_read_ratio=ratio, truth_underrepresented=positive,
                           predicted_assembly_read_ratio=finite_nonnegative(by_family[family]["assembly_read_ratio"], "ratio"),
                           native_status=status, predicted_underrepresented=call,
                           outcome="TP" if positive and call else "FN" if positive else "FP" if call else "TN"))
    return result
