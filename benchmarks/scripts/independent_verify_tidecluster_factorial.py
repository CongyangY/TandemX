#!/usr/bin/env python3
"""Independently recompute frozen TideCluster factorial accuracy endpoints."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


CHECKED_METRICS = (
    "truth_array_count",
    "predicted_array_count",
    "matched_array_count",
    "array_recall",
    "array_precision",
    "base_union_recall",
    "base_union_precision",
    "matched_boundary_mae_bp",
    "matched_period_mae_bp",
    "truth_family_count",
    "recovered_family_count",
    "cyclic_monomer_recall",
    "distinct_consensus_count",
    "homologous_consensus_fraction",
    "operational_family_count",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence).upper()
                name = value[1:].split()[0]
                if not name or name in records:
                    raise ValueError(f"invalid FASTA header: {path}:{line_number}")
                sequence = []
            elif name is None:
                raise ValueError(f"sequence before FASTA header: {path}:{line_number}")
            else:
                sequence.append(value)
    if name is not None:
        records[name] = "".join(sequence).upper()
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def maximum_matching(edges: list[list[int]]) -> dict[int, int]:
    owners: dict[int, int] = {}

    def augment(left: int, visited: set[int]) -> bool:
        for right in edges[left]:
            if right in visited:
                continue
            visited.add(right)
            if right not in owners or augment(owners[right], visited):
                owners[right] = left
                return True
        return False

    for left in range(len(edges)):
        augment(left, set())
    return {left: right for right, left in owners.items()}


def interval_iou(left: dict[str, Any], right: dict[str, Any]) -> float:
    if left["chrom"] != right["chrom"]:
        return 0.0
    intersection = max(0, min(left["end"], right["end"]) - max(left["start"], right["start"]))
    union = left["end"] - left["start"] + right["end"] - right["start"] - intersection
    return intersection / union


def merged_intervals(rows: list[dict[str, Any]]) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in sorted(rows, key=lambda item: (item["chrom"], item["start"], item["end"])):
        intervals = result[row["chrom"]]
        if intervals and row["start"] <= intervals[-1][1]:
            intervals[-1] = (intervals[-1][0], max(intervals[-1][1], row["end"]))
        else:
            intervals.append((row["start"], row["end"]))
    return result


def base_union_metrics(
    predicted: list[dict[str, Any]], truth: list[dict[str, Any]]
) -> tuple[float, float]:
    p_merged, t_merged = merged_intervals(predicted), merged_intervals(truth)
    p_total = sum(end - start for values in p_merged.values() for start, end in values)
    t_total = sum(end - start for values in t_merged.values() for start, end in values)
    overlap = 0
    for chrom in p_merged.keys() & t_merged.keys():
        p_values, t_values = p_merged[chrom], t_merged[chrom]
        i = j = 0
        while i < len(p_values) and j < len(t_values):
            overlap += max(
                0,
                min(p_values[i][1], t_values[j][1])
                - max(p_values[i][0], t_values[j][0]),
            )
            if p_values[i][1] < t_values[j][1]:
                i += 1
            else:
                j += 1
    return (
        overlap / t_total if t_total else math.nan,
        overlap / p_total if p_total else math.nan,
    )


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_monomer(sequence: str) -> str:
    sequence = sequence.upper()
    reverse = reverse_complement(sequence)
    return min(
        strand[offset:] + strand[:offset]
        for strand in (sequence, reverse)
        for offset in range(len(sequence))
    )


def reaches_cyclic_threshold(truth: str, prediction: str, threshold: float = 0.9) -> bool:
    import edlib

    truth, prediction = truth.upper(), prediction.upper()
    maximum = max(len(truth), len(prediction))
    limit = math.floor((1 - threshold + 1e-12) * maximum)
    if abs(len(truth) - len(prediction)) > limit:
        return False
    reference = truth.replace("N", "X")
    for strand in (prediction, reverse_complement(prediction)):
        strand = strand.replace("N", "Y")
        for offset in range(len(strand)):
            rotated = strand[offset:] + strand[:offset]
            result = edlib.align(reference, rotated, mode="NW", task="distance", k=limit)
            if result["editDistance"] >= 0:
                return True
    return False


def recompute_run(genome_dir: Path, evaluation_dir: Path) -> dict[str, float | int]:
    truth_rows = read_tsv(genome_dir / "truth_copy_number.tsv")
    truth = [
        {
            "chrom": row["chrom"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "period": int(row["period"]),
        }
        for row in truth_rows
    ]
    normalized_rows = read_tsv(evaluation_dir / "normalized_arrays.tsv")
    predicted = [
        {
            "chrom": row["sequence_id"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "period": int(row["period"]),
            "family_id": row["family_id"],
            "consensus": row["consensus_sequence"],
        }
        for row in normalized_rows
    ]
    truth_by_chrom: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(truth):
        truth_by_chrom[row["chrom"]].append(index)
    edges = [
        [
            index
            for index in truth_by_chrom[row["chrom"]]
            if interval_iou(row, truth[index]) >= 0.5
            and abs(row["period"] - truth[index]["period"])
            <= max(2, round(0.02 * truth[index]["period"]))
        ]
        for row in predicted
    ]
    matching = maximum_matching(edges)
    base_recall, base_precision = base_union_metrics(predicted, truth)
    if matching:
        boundary_mae = statistics.fmean(
            (
                abs(predicted[left]["start"] - truth[right]["start"])
                + abs(predicted[left]["end"] - truth[right]["end"])
            )
            / 2
            for left, right in matching.items()
        )
        period_mae = statistics.fmean(
            abs(predicted[left]["period"] - truth[right]["period"])
            for left, right in matching.items()
        )
    else:
        boundary_mae = period_mae = math.nan

    founders = read_fasta(genome_dir / "catalogue.fa")
    sequences = sorted(
        {canonical_monomer(row["consensus"]) for row in predicted if row["consensus"]}
    )
    identifiers = sorted(founders)
    family_edges = [
        [
            index
            for index, sequence in enumerate(sequences)
            if reaches_cyclic_threshold(founders[name], sequence)
        ]
        for name in identifiers
    ]
    family_matching = maximum_matching(family_edges)
    supported = {index for values in family_edges for index in values}
    return {
        "truth_array_count": len(truth),
        "predicted_array_count": len(predicted),
        "matched_array_count": len(matching),
        "array_recall": len(matching) / len(truth) if truth else math.nan,
        "array_precision": len(matching) / len(predicted) if predicted else math.nan,
        "base_union_recall": base_recall,
        "base_union_precision": base_precision,
        "matched_boundary_mae_bp": boundary_mae,
        "matched_period_mae_bp": period_mae,
        "truth_family_count": len(identifiers),
        "recovered_family_count": len(family_matching),
        "cyclic_monomer_recall": (
            len(family_matching) / len(identifiers) if identifiers else math.nan
        ),
        "distinct_consensus_count": len(sequences),
        "homologous_consensus_fraction": (
            len(supported) / len(sequences) if sequences else math.nan
        ),
        "operational_family_count": len({row["family_id"] for row in predicted}),
    }


def same_value(observed: Any, expected: Any) -> bool:
    if observed is None:
        return isinstance(expected, float) and math.isnan(expected)
    if isinstance(expected, int):
        return int(observed) == expected
    return math.isclose(float(observed), float(expected), rel_tol=1e-12, abs_tol=1e-12)


def verify(config_path: Path, result_dir: Path, output: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    receipt = json.loads((result_dir / "run_receipt.json").read_text(encoding="utf-8"))
    summary_rows = read_tsv(result_dir / "summary.tsv")
    expected_keys = {
        (int(dataset["seed"]), setting)
        for dataset in config["datasets"]
        for setting in config["settings"]
    }
    observed_keys = {(int(row["seed"]), row["setting"]) for row in summary_rows}
    failures: list[str] = []
    if receipt.get("complete") is not True:
        failures.append("run_receipt_not_complete")
    actual_summary_sha256 = sha256(result_dir / "summary.tsv")
    if receipt.get("summary_sha256") != actual_summary_sha256:
        failures.append("run_receipt_summary_hash_changed")
    if receipt.get("run_count") != len(expected_keys):
        failures.append("run_count_changed")
    if observed_keys != expected_keys:
        failures.append("summary_run_keys_changed")
    cell_fates_path = result_dir / "cell_fates.tsv"
    stage_fates_path = result_dir / "stage_fates.tsv"
    cell_fates = read_tsv(cell_fates_path) if cell_fates_path.is_file() else []
    stage_fates = read_tsv(stage_fates_path) if stage_fates_path.is_file() else []
    has_unavailable = any(row.get("status") != "ok" for row in summary_rows)
    if has_unavailable and not cell_fates:
        failures.append("missing_cell_fates_for_unavailable_runs")
    if has_unavailable and not stage_fates:
        failures.append("missing_stage_fates_for_unavailable_runs")
    if cell_fates:
        fate_keys = {(int(row["seed"]), row["setting"]) for row in cell_fates}
        if fate_keys != expected_keys or len(cell_fates) != len(expected_keys):
            failures.append("cell_fate_keys_changed")
        fate_by_key = {
            (int(row["seed"]), row["setting"]): row for row in cell_fates
        }
        for row in summary_rows:
            key = (int(row["seed"]), row["setting"])
            if fate_by_key.get(key, {}).get("status") != row.get("status"):
                failures.append(f"seed{key[0]}/{key[1]}:cell_fate_status_mismatch")
        if receipt.get("cell_fates_sha256") != sha256(cell_fates_path):
            failures.append("run_receipt_cell_fates_hash_changed")
    if stage_fates:
        expected_stage_keys = {
            (int(dataset["seed"]), setting, component)
            for dataset in config["datasets"]
            for setting in config["settings"]
            for component in ("tidehunter", "clustering")
        }
        observed_stage_keys = {
            (int(row["seed"]), row["setting"], row["component"])
            for row in stage_fates
        }
        if observed_stage_keys != expected_stage_keys or len(stage_fates) != len(
            expected_stage_keys
        ):
            failures.append("stage_fate_keys_changed")
        if receipt.get("stage_fates_sha256") != sha256(stage_fates_path):
            failures.append("run_receipt_stage_fates_hash_changed")
    environment_path = result_dir / "environment.json"
    if environment_path.is_file():
        environment = json.loads(environment_path.read_text(encoding="utf-8"))
        for name, expected_hash in environment.get("parent_failure", {}).get(
            "artifacts", {}
        ).items():
            parent_copy = result_dir / "parent_failure_snapshot" / name
            if not parent_copy.is_file() or sha256(parent_copy) != expected_hash:
                failures.append(f"parent_failure_snapshot_changed:{name}")
    dataset_by_seed = {int(row["seed"]): row for row in config["datasets"]}
    for seed, dataset in dataset_by_seed.items():
        genome_dir = Path(dataset["genome_dir"])
        for filename, key in (
            ("genome.fa", "genome_sha256"),
            ("catalogue.fa", "catalogue_sha256"),
            ("truth_copy_number.tsv", "truth_sha256"),
        ):
            if sha256(genome_dir / filename) != dataset[key]:
                failures.append(f"seed{seed}:frozen_input_hash_changed:{filename}")
    checks: list[dict[str, Any]] = []
    for summary in summary_rows:
        seed, setting = int(summary["seed"]), summary["setting"]
        if summary.get("status") != "ok":
            nonmissing = [
                metric for metric in CHECKED_METRICS if summary.get(metric, "") != ""
            ]
            mismatches = []
            if nonmissing:
                mismatches.append("unavailable_run_has_accuracy_values:" + ",".join(nonmissing))
            if "unavailable" not in summary.get("warning", ""):
                mismatches.append("unavailable_run_missing_explicit_warning")
            if mismatches:
                failures.append(f"seed{seed}/{setting}:" + ",".join(mismatches))
            checks.append(
                {
                    "seed": seed,
                    "setting": setting,
                    "status": summary.get("status"),
                    "verification_passed": not mismatches,
                    "mismatches": mismatches,
                    "recomputed": None,
                }
            )
            continue
        evaluation_dir = result_dir / f"seed{seed}" / setting / "evaluation"
        expected = recompute_run(Path(dataset_by_seed[seed]["genome_dir"]), evaluation_dir)
        stored = json.loads((evaluation_dir / "metrics.json").read_text(encoding="utf-8"))
        mismatches = []
        for metric in CHECKED_METRICS:
            if not same_value(stored.get(metric), expected[metric]):
                mismatches.append(f"metrics.json:{metric}")
            summary_value: Any = summary[metric]
            if summary_value == "":
                summary_value = None
            if not same_value(summary_value, expected[metric]):
                mismatches.append(f"summary.tsv:{metric}")
        evaluation_receipt = json.loads(
            (evaluation_dir / "evaluation_receipt.json").read_text(encoding="utf-8")
        )
        for name, digest in evaluation_receipt.get("outputs", {}).items():
            if sha256(evaluation_dir / name) != digest:
                mismatches.append(f"output_hash:{name}")
        if mismatches:
            failures.append(f"seed{seed}/{setting}:" + ",".join(mismatches))
        checks.append(
            {
                "seed": seed,
                "setting": setting,
                "status": "ok",
                "verification_passed": not mismatches,
                "mismatches": mismatches,
                "recomputed": expected,
            }
        )
    successful_run_count = sum(row.get("status") == "ok" for row in summary_rows)
    if receipt.get("successful_run_count") != successful_run_count:
        failures.append("successful_run_count_changed")
    if "accuracy_complete" in receipt and receipt.get("accuracy_complete") != (
        successful_run_count == len(expected_keys)
    ):
        failures.append("accuracy_complete_flag_changed")
    payload = {
        "schema_version": 1,
        "verification_passed": not failures,
        "config_sha256": sha256(config_path),
        "summary_sha256": actual_summary_sha256,
        "checked_run_count": len(checks),
        "successful_accuracy_run_count": successful_run_count,
        "unavailable_accuracy_run_count": len(checks) - successful_run_count,
        "checked_metrics": list(CHECKED_METRICS),
        "failures": failures,
        "runs": checks,
        "implementation_boundary": (
            "does_not_import_TideCluster_normalization_or_primary_evaluation_modules;"
            "uses_independent_TSV_FASTA_parsers_interval_matching_union_and_cyclic_edlib"
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = verify(args.config, args.result_dir, args.output)
    return 0 if payload["verification_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
