"""Combine frozen per-k and direction-only orthogonal abundance results."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file


RESIDUAL_K = "orthogonal_supports_residual_collapse_at_this_k"
BIAS_K = "orthogonal_supports_quantification_bias_at_this_k"
ONT_RESIDUAL = "ONT_direction_supports_residual_collapse"


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader)


def _candidate_rows(path: Path, expected_species: str) -> dict[str, dict[str, str]]:
    rows = read_table(path)
    if not rows:
        raise ValueError(f"Empty input: {path}")
    species = {row["species_key"] for row in rows}
    if species != {expected_species}:
        raise ValueError(f"Unexpected species in {path}: {sorted(species)}")
    candidates = {
        row["family_id"]: row
        for row in rows
        if row["preorthogonal_candidate"] == "true"
        and row["eligibility"] == "eligible"
    }
    if not candidates:
        raise ValueError(f"No eligible frozen candidates in {path}")
    return candidates


def finalize(
    species_key: str,
    k21_metrics: Path,
    k31_metrics: Path,
    output_tsv: Path,
    summary_json: Path,
    *,
    ont_direction: Path | None = None,
) -> dict[str, object]:
    """Apply the frozen cross-k/cross-platform interpretation without retuning."""
    if species_key not in {"ey15", "macadamia"}:
        raise ValueError("species_key must be ey15 or macadamia")
    if species_key == "macadamia" and ont_direction is None:
        raise ValueError("Macadamia finalization requires direction-only ONT results")
    if species_key == "ey15" and ont_direction is not None:
        raise ValueError("Ey15 finalization does not accept an ONT result")
    for path in (output_tsv, summary_json):
        if path.exists():
            raise ValueError(f"Refusing to overwrite existing output: {path}")

    k21 = _candidate_rows(k21_metrics, species_key)
    k31 = _candidate_rows(k31_metrics, species_key)
    if set(k21) != set(k31):
        raise ValueError("Frozen candidate sets differ between k=21 and k=31")
    ont: dict[str, dict[str, str]] = {}
    if ont_direction is not None:
        ont = {
            row["family_id"]: row
            for row in read_table(ont_direction)
            if row["preorthogonal_candidate"] == "true"
        }
        if set(ont) != set(k21):
            raise ValueError("Frozen candidate sets differ between Illumina and ONT")

    rows: list[dict[str, object]] = []
    for family_id in sorted(k21):
        row21, row31 = k21[family_id], k31[family_id]
        states = (row21["k_result"], row31["k_result"])
        illumina = (
            "stable_illumina_supports_residual_collapse"
            if states == (RESIDUAL_K, RESIDUAL_K)
            else "stable_illumina_supports_quantification_bias"
            if states == (BIAS_K, BIAS_K)
            else "illumina_unresolved_across_k"
        )
        ont_state = ont.get(family_id, {}).get("ONT_direction_result", "not_available")
        if illumina == "stable_illumina_supports_quantification_bias":
            final = "orthogonal_supports_quantification_bias"
        elif illumina == "stable_illumina_supports_residual_collapse":
            if species_key == "ey15" or ont_state == ONT_RESIDUAL:
                final = "orthogonal_supports_residual_collapse"
            else:
                final = "unresolved"
        else:
            final = "unresolved"

        estimates = [float(row21["orthogonal_estimated_bp"]), float(row31["orthogonal_estimated_bp"])]
        deficits = [
            float(row21["estimated_under_representation_bp"]),
            float(row31["estimated_under_representation_bp"]),
        ]
        assembly_ratios = [
            float(row21["newer_assembly_orthogonal_ratio"]),
            float(row31["newer_assembly_orthogonal_ratio"]),
        ]
        hifi_ratios = [
            float(row21["hifi_orthogonal_ratio"]),
            float(row31["hifi_orthogonal_ratio"]),
        ]
        rows.append(
            {
                "species_key": species_key,
                "family_id": family_id,
                "newer_assembly_bp": row21["newer_assembly_bp"],
                "frozen_hifi_estimated_bp": row21["frozen_hifi_estimated_bp"],
                "k21_orthogonal_estimated_bp": row21["orthogonal_estimated_bp"],
                "k31_orthogonal_estimated_bp": row31["orthogonal_estimated_bp"],
                "orthogonal_estimated_bp_min": min(estimates),
                "orthogonal_estimated_bp_max": max(estimates),
                "estimated_under_representation_bp_min": min(deficits),
                "estimated_under_representation_bp_max": max(deficits),
                "newer_assembly_orthogonal_ratio_min": min(assembly_ratios),
                "newer_assembly_orthogonal_ratio_max": max(assembly_ratios),
                "hifi_orthogonal_ratio_min": min(hifi_ratios),
                "hifi_orthogonal_ratio_max": max(hifi_ratios),
                "k21_result": row21["k_result"],
                "k31_result": row31["k_result"],
                "illumina_cross_k_result": illumina,
                "ont_direction_result": ont_state,
                "final_interpretation": final,
                "magnitude_boundary": "read_assembly_deficit_not_physical_missing_bp_truth",
            }
        )

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    states = sorted({str(row["final_interpretation"]) for row in rows})
    summary: dict[str, object] = {
        "schema_version": 1,
        "status": "complete_frozen_cross_platform_interpretation",
        "species_key": species_key,
        "candidate_families": len(rows),
        "final_interpretation_counts": {
            state: sum(row["final_interpretation"] == state for row in rows)
            for state in states
        },
        "classification_boundary": "binary_interpretation_is_separate_from_continuous_magnitude",
        "magnitude_boundary": "ranges_are_read_assembly_deficits_not_missing_bp_truth_or_prediction_error",
        "inputs": {
            "k21_metrics_sha256": digest_file(k21_metrics),
            "k31_metrics_sha256": digest_file(k31_metrics),
            "ont_direction_sha256": digest_file(ont_direction) if ont_direction else None,
        },
        "output_tsv_sha256": digest_file(output_tsv),
    }
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-key", required=True, choices=("ey15", "macadamia"))
    parser.add_argument("--k21-metrics", required=True, type=Path)
    parser.add_argument("--k31-metrics", required=True, type=Path)
    parser.add_argument("--ont-direction", type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    finalize(
        args.species_key,
        args.k21_metrics,
        args.k31_metrics,
        args.output_tsv,
        args.summary,
        ont_direction=args.ont_direction,
    )


if __name__ == "__main__":
    main()
