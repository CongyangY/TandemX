"""Prepare and evaluate frozen Illumina k-mer abundance validation assets."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median

from benchmarks.challenge.schema import digest_file
from tandemx.quantify.mvp import (
    empirical_quantile,
    family_kmer_membership,
    monomer_kmer_counts,
    read_monomer_fasta,
    read_single_copy_kmers,
)
from tandemx.utils.kmers import is_low_complexity_kmer


def diagnostic_kmers(catalogue: Path, k: int) -> tuple[list, dict[str, dict[str, int]]]:
    monomers = list(read_monomer_fasta(catalogue))
    membership = family_kmer_membership(monomers, k)
    diagnostic = {
        monomer.family_id: {
            word: multiplicity
            for word, multiplicity in monomer_kmer_counts(monomer.sequence, k).items()
            if len(membership[word]) == 1 and not is_low_complexity_kmer(word)
        }
        for monomer in monomers
    }
    return monomers, diagnostic


def prepare_targets(
    catalogue: Path,
    controls_tsv: Path,
    target_fasta: Path,
    target_map_tsv: Path,
    receipt_json: Path,
    *,
    k: int,
) -> dict[str, object]:
    """Write one unique FASTA record per diagnostic or control target."""
    for path in (target_fasta, target_map_tsv, receipt_json):
        if path.exists():
            raise ValueError(f"Refusing to overwrite existing output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    monomers, diagnostic = diagnostic_kmers(catalogue, k)
    controls = read_single_copy_kmers(controls_tsv, k)
    repeat_words = {word for family in diagnostic.values() for word in family}
    overlap = repeat_words.intersection(controls)
    if overlap:
        raise ValueError(f"Control/diagnostic overlap: {sorted(overlap)[:3]}")
    targets = sorted(repeat_words | set(controls))

    with target_fasta.open("w", encoding="utf-8") as fasta:
        for index, word in enumerate(targets, start=1):
            fasta.write(f">target_{index:07d}\n{word}\n")
    with target_map_tsv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["kmer", "target_type", "family_id", "monomer_multiplicity", "expected_copy_number"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for family_id in sorted(diagnostic):
            for word, multiplicity in sorted(diagnostic[family_id].items()):
                writer.writerow(
                    {
                        "kmer": word,
                        "target_type": "family_diagnostic",
                        "family_id": family_id,
                        "monomer_multiplicity": multiplicity,
                        "expected_copy_number": "NA",
                    }
                )
        for word, expected in sorted(controls.items()):
            writer.writerow(
                {
                    "kmer": word,
                    "target_type": "single_copy_control",
                    "family_id": "NA",
                    "monomer_multiplicity": "NA",
                    "expected_copy_number": f"{expected:.12g}",
                }
            )
    receipt: dict[str, object] = {
        "schema_version": 1,
        "status": "complete",
        "k": k,
        "family_count": len(monomers),
        "families_with_diagnostic_kmers": sum(bool(words) for words in diagnostic.values()),
        "diagnostic_kmers": len(repeat_words),
        "single_copy_controls": len(controls),
        "total_unique_targets": len(targets),
        "catalogue_sha256": digest_file(catalogue),
        "controls_sha256": digest_file(controls_tsv),
        "target_fasta_sha256": digest_file(target_fasta),
        "target_map_tsv_sha256": digest_file(target_map_tsv),
        "boundary": "Target preparation fixes the existing TandemX family-exclusive k-mer definition; it does not inspect orthogonal counts.",
    }
    receipt_json.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def read_kmc_dump(path: Path, k: int) -> dict[str, int]:
    counts: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.split()
            if len(parts) != 2 or len(parts[0]) != k or set(parts[0]) - set("ACGT"):
                raise ValueError(f"Malformed KMC dump at line {line_number}")
            word = parts[0]
            value = int(parts[1])
            if value < 2:
                raise ValueError("Frozen KMC extraction expects input database minimum count 2")
            if word in counts:
                raise ValueError(f"Duplicate KMC k-mer: {word}")
            counts[word] = value
    return counts


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader)


def median_absolute_deviation(values: list[float], center: float) -> float:
    return float(median(abs(value - center) for value in values)) if values else 0.0


def evaluate(
    species_key: str,
    config_json: Path,
    catalogue: Path,
    controls_tsv: Path,
    kmc_dump: Path,
    prior_family_metrics: Path,
    output_tsv: Path,
    summary_json: Path,
    *,
    k: int,
) -> dict[str, object]:
    """Recompute frozen abundance from independent KMC target counts."""
    for path in (output_tsv, summary_json):
        if path.exists():
            raise ValueError(f"Refusing to overwrite existing output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    config = json.loads(config_json.read_text(encoding="utf-8"))
    candidates = set(config["preorthogonal_hifi_newer_assembly_deficit_candidates"][species_key])
    threshold = float(config["binary_threshold"])
    bias_fold = float(config["quantification_bias_fold_threshold"])
    monomers, diagnostic = diagnostic_kmers(catalogue, k)
    controls = read_single_copy_kmers(controls_tsv, k)
    counts = read_kmc_dump(kmc_dump, k)
    control_depths = [counts.get(word, 0) / expected for word, expected in controls.items()]
    if not control_depths:
        raise ValueError("No single-copy controls")
    control_mean = math.fsum(control_depths) / len(control_depths)
    control_median = float(median(control_depths))
    control_mad = median_absolute_deviation(control_depths, control_median)
    control_zero_fraction = sum(value == 0 for value in control_depths) / len(control_depths)
    if control_mean <= 0:
        raise ValueError("Single-copy control mean is zero")

    prior_rows = {row["family_id"]: row for row in read_table(prior_family_metrics)}
    if set(prior_rows) != {monomer.family_id for monomer in monomers}:
        raise ValueError("Prior family metrics and catalogue family sets differ")
    rows: list[dict[str, object]] = []
    for monomer in monomers:
        words = diagnostic[monomer.family_id]
        depths = [counts.get(word, 0) / multiplicity for word, multiplicity in words.items()]
        family_median = float(median(depths)) if depths else 0.0
        family_mad = median_absolute_deviation(depths, family_median)
        copy_number = family_median / control_mean
        orthogonal_bp = copy_number * len(monomer.sequence)
        prior = prior_rows[monomer.family_id]
        hifi_bp = float(prior["read_estimated_bp"])
        new_bp = float(prior["new_assembly_bp"])
        new_orth_ratio = new_bp / orthogonal_bp if orthogonal_bp > 0 else math.inf
        hifi_orth_ratio = hifi_bp / orthogonal_bp if orthogonal_bp > 0 else math.inf
        is_candidate = monomer.family_id in candidates
        if is_candidate and orthogonal_bp <= 0:
            k_result = "technical_failure_no_orthogonal_estimate"
        elif is_candidate and new_orth_ratio < threshold:
            k_result = "orthogonal_supports_residual_collapse_at_this_k"
        elif is_candidate and hifi_orth_ratio >= bias_fold:
            k_result = "orthogonal_supports_quantification_bias_at_this_k"
        elif is_candidate:
            k_result = "unresolved_at_this_k"
        elif new_orth_ratio >= threshold:
            k_result = "negative_context_no_residual_deficit"
        else:
            k_result = "negative_context_orthogonal_deficit"
        rows.append(
            {
                "species_key": species_key,
                "k": k,
                "family_id": monomer.family_id,
                "monomer_length": len(monomer.sequence),
                "diagnostic_kmer_count": len(words),
                "diagnostic_kmers_observed_at_least_twice": sum(counts.get(word, 0) >= 2 for word in words),
                "median_diagnostic_kmer_depth": family_median,
                "diagnostic_depth_mad": family_mad,
                "empirical_haploid_depth_mean": control_mean,
                "orthogonal_estimated_copy_number": copy_number,
                "orthogonal_estimated_bp": orthogonal_bp,
                "frozen_hifi_estimated_bp": hifi_bp,
                "newer_assembly_bp": new_bp,
                "newer_assembly_orthogonal_ratio": new_orth_ratio,
                "hifi_orthogonal_ratio": hifi_orth_ratio,
                "estimated_under_representation_bp": max(orthogonal_bp - new_bp, 0.0),
                "preorthogonal_candidate": str(is_candidate).lower(),
                "eligibility": prior["eligibility"],
                "k_result": k_result,
                "warning": "KMC_ci2_treats_unreported_counts_as_zero_or_one;assembly_derived_controls_do_not_make_assembly_copy_truth",
            }
        )
    fieldnames = list(rows[0])
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    eligible_rows = [row for row in rows if row["eligibility"] == "eligible"]
    candidate_rows = [row for row in eligible_rows if row["preorthogonal_candidate"] == "true"]
    summary: dict[str, object] = {
        "schema_version": 1,
        "status": "complete_single_k_result_not_final_cross_k_interpretation",
        "species_key": species_key,
        "k": k,
        "catalogue_families": len(rows),
        "eligible_families": len(eligible_rows),
        "preorthogonal_candidates": len(candidate_rows),
        "control_count": len(controls),
        "control_mean_depth": control_mean,
        "control_median_depth": control_median,
        "control_depth_mad": control_mad,
        "control_zero_fraction": control_zero_fraction,
        "kmc_reported_target_kmers": len(counts),
        "candidate_k_result_counts": {
            state: sum(row["k_result"] == state for row in candidate_rows)
            for state in sorted({str(row["k_result"]) for row in candidate_rows})
        },
        "eligible_negative_context_orthogonal_deficits": sum(
            row["k_result"] == "negative_context_orthogonal_deficit" for row in eligible_rows
        ),
        "magnitude_boundary": "estimated_under_representation_is_a_read_assembly_deficit_not_missing_bp_truth_or_prediction_error",
        "inputs": {
            "config_sha256": digest_file(config_json),
            "catalogue_sha256": digest_file(catalogue),
            "controls_sha256": digest_file(controls_tsv),
            "kmc_dump_sha256": digest_file(kmc_dump),
            "prior_family_metrics_sha256": digest_file(prior_family_metrics),
        },
        "output_tsv_sha256": digest_file(output_tsv),
    }
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--catalogue", required=True, type=Path)
    prepare.add_argument("--controls", required=True, type=Path)
    prepare.add_argument("--target-fasta", required=True, type=Path)
    prepare.add_argument("--target-map", required=True, type=Path)
    prepare.add_argument("--receipt", required=True, type=Path)
    prepare.add_argument("--k", required=True, type=int)

    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument("--species-key", required=True, choices=("ey15", "macadamia"))
    evaluate_parser.add_argument("--config", required=True, type=Path)
    evaluate_parser.add_argument("--catalogue", required=True, type=Path)
    evaluate_parser.add_argument("--controls", required=True, type=Path)
    evaluate_parser.add_argument("--kmc-dump", required=True, type=Path)
    evaluate_parser.add_argument("--prior-family-metrics", required=True, type=Path)
    evaluate_parser.add_argument("--output-tsv", required=True, type=Path)
    evaluate_parser.add_argument("--summary", required=True, type=Path)
    evaluate_parser.add_argument("--k", required=True, type=int)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare_targets(
            args.catalogue,
            args.controls,
            args.target_fasta,
            args.target_map,
            args.receipt,
            k=args.k,
        )
    else:
        evaluate(
            args.species_key,
            args.config,
            args.catalogue,
            args.controls,
            args.kmc_dump,
            args.prior_family_metrics,
            args.output_tsv,
            args.summary,
            k=args.k,
        )


if __name__ == "__main__":
    main()
