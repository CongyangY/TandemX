"""Prepare and evaluate the frozen direction-only Macadamia ONT occupancy audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median

from benchmarks.challenge.schema import digest_file
from tandemx.quantify.mvp import read_monomer_fasta


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader)


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Refusing to write an empty table")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ValueError(f"Refusing to overwrite existing output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def prepare_templates(
    catalogue: Path,
    prior_family_metrics: Path,
    output_fasta: Path,
    receipt_json: Path,
    *,
    template_length: int,
) -> dict[str, object]:
    """Tandemize all previously eligible representatives without result selection."""
    if template_length < 1:
        raise ValueError("template_length must be positive")
    if output_fasta.exists() or receipt_json.exists():
        raise ValueError("Refusing to overwrite existing template outputs")
    eligible = {
        row["family_id"]
        for row in read_table(prior_family_metrics)
        if row["eligibility"] == "eligible"
    }
    monomers = {record.family_id: record.sequence for record in read_monomer_fasta(catalogue)}
    if not eligible <= set(monomers):
        raise ValueError("Eligible family is missing from catalogue")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as handle:
        for family_id in sorted(eligible):
            sequence = monomers[family_id]
            repeats = (template_length + len(sequence) - 1) // len(sequence)
            tandem = (sequence * repeats)[:template_length]
            handle.write(f">{family_id}\n{tandem}\n")
    receipt: dict[str, object] = {
        "schema_version": 1,
        "status": "complete",
        "eligible_families": len(eligible),
        "template_length_bp_per_family": template_length,
        "total_template_bp": template_length * len(eligible),
        "catalogue_sha256": digest_file(catalogue),
        "prior_family_metrics_sha256": digest_file(prior_family_metrics),
        "template_fasta_sha256": digest_file(output_fasta),
        "boundary": "Templates include all frozen eligible families and do not use ONT outcomes.",
    }
    receipt_json.parent.mkdir(parents=True, exist_ok=True)
    receipt_json.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start > end:
            total += end - start
            start, end = next_start, next_end
        else:
            end = max(end, next_end)
    return total + end - start


def summarize_paf(
    paf: Path,
    template_fasta: Path,
    output_tsv: Path,
    summary_json: Path,
    *,
    total_library_bases: int,
    genome_size_bp: int,
    minimum_block_bp: int,
    minimum_identity: float,
) -> dict[str, object]:
    """Sum accepted primary query intervals, excluding multi-family assignments."""
    if total_library_bases < 1 or genome_size_bp < 1 or minimum_block_bp < 1:
        raise ValueError("Library, genome, and block sizes must be positive")
    if not 0 < minimum_identity <= 1:
        raise ValueError("minimum_identity must be in (0,1]")
    families = [record.family_id for record in read_monomer_fasta(template_fasta)]
    family_set = set(families)
    accepted_bp = {family: 0 for family in families}
    accepted_reads = {family: 0 for family in families}
    ambiguous_reads = accepted_paf_rows = rejected_paf_rows = 0
    current_query: str | None = None
    current: dict[str, list[tuple[int, int]]] = {}
    completed_queries: set[str] = set()

    def flush() -> None:
        nonlocal ambiguous_reads
        if not current:
            return
        if len(current) != 1:
            ambiguous_reads += 1
            return
        family, intervals = next(iter(current.items()))
        accepted_bp[family] += union_length(intervals)
        accepted_reads[family] += 1

    with paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number}")
            query = fields[0]
            if current_query is None:
                current_query = query
            elif query != current_query:
                flush()
                completed_queries.add(current_query)
                if query in completed_queries:
                    raise ValueError("PAF is not grouped by query name")
                current_query = query
                current = {}
            target = fields[5]
            if target not in family_set:
                raise ValueError(f"Unexpected target family at line {line_number}: {target}")
            qstart, qend = int(fields[2]), int(fields[3])
            nmatch, block = int(fields[9]), int(fields[10])
            tags = fields[12:]
            is_primary = "tp:A:P" in tags
            if (
                not is_primary
                or block < minimum_block_bp
                or block <= 0
                or nmatch / block < minimum_identity
            ):
                rejected_paf_rows += 1
                continue
            accepted_paf_rows += 1
            current.setdefault(target, []).append((qstart, qend))
    flush()
    rows = [
        {
            "family_id": family,
            "accepted_query_bp": accepted_bp[family],
            "accepted_read_count": accepted_reads[family],
            "library_fraction": accepted_bp[family] / total_library_bases,
            "raw_mapping_estimated_bp": accepted_bp[family]
            / total_library_bases
            * genome_size_bp,
        }
        for family in sorted(families)
    ]
    write_table(output_tsv, rows)
    summary: dict[str, object] = {
        "schema_version": 1,
        "status": "complete",
        "families": len(families),
        "total_library_bases": total_library_bases,
        "genome_size_bp": genome_size_bp,
        "minimum_alignment_block_bp": minimum_block_bp,
        "minimum_identity": minimum_identity,
        "accepted_primary_paf_rows": accepted_paf_rows,
        "rejected_paf_rows": rejected_paf_rows,
        "ambiguous_multi_family_reads_excluded": ambiguous_reads,
        "accepted_query_bp": sum(accepted_bp.values()),
        "paf_sha256": digest_file(paf),
        "template_fasta_sha256": digest_file(template_fasta),
        "output_tsv_sha256": digest_file(output_tsv),
        "boundary": "Mapping occupancy is direction-only and not physical copy-number truth.",
    }
    if summary_json.exists():
        raise ValueError(f"Refusing to overwrite existing output: {summary_json}")
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def calibrate_hifi_mapping(
    hifi_occupancy_tsv: Path,
    prior_family_metrics: Path,
    output_json: Path,
) -> dict[str, object]:
    """Freeze one global mapping-efficiency factor before ONT inspection."""
    if output_json.exists():
        raise ValueError(f"Refusing to overwrite existing output: {output_json}")
    occupancy = {row["family_id"]: row for row in read_table(hifi_occupancy_tsv)}
    efficiencies: list[tuple[str, float]] = []
    for row in read_table(prior_family_metrics):
        if row["eligibility"] != "eligible":
            continue
        family = row["family_id"]
        mapping_bp = float(occupancy[family]["raw_mapping_estimated_bp"])
        kmer_bp = float(row["read_estimated_bp"])
        if mapping_bp > 0 and kmer_bp > 0:
            efficiencies.append((family, mapping_bp / kmer_bp))
    if not efficiencies:
        raise ValueError("No positive eligible families for HiFi mapping calibration")
    factor = float(median(value for _, value in efficiencies))
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError("Invalid mapping-efficiency factor")
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "frozen_before_ONT_result_inspection",
        "positive_eligible_families": len(efficiencies),
        "global_median_raw_HiFi_mapping_over_frozen_HiFi_kmer": factor,
        "family_efficiencies": [
            {"family_id": family, "efficiency": value}
            for family, value in efficiencies
        ],
        "hifi_occupancy_tsv_sha256": digest_file(hifi_occupancy_tsv),
        "prior_family_metrics_sha256": digest_file(prior_family_metrics),
        "boundary": "One global technical mapping correction; no family-specific or ONT-result fitting.",
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def evaluate_ont(
    ont_occupancy_tsv: Path,
    calibration_json: Path,
    prior_family_metrics: Path,
    config_json: Path,
    output_tsv: Path,
    summary_json: Path,
) -> dict[str, object]:
    """Apply frozen calibration and report direction across genome-size bounds."""
    config = json.loads(config_json.read_text(encoding="utf-8"))
    sizes = [int(value) for value in config["macadamia_genome_size_sensitivity_bp"]]
    threshold = float(config["binary_threshold"])
    candidates = set(config["preorthogonal_hifi_newer_assembly_deficit_candidates"]["macadamia"])
    calibration = json.loads(calibration_json.read_text(encoding="utf-8"))
    if calibration["status"] != "frozen_before_ONT_result_inspection":
        raise ValueError("HiFi mapping calibration was not frozen before ONT inspection")
    efficiency = float(calibration["global_median_raw_HiFi_mapping_over_frozen_HiFi_kmer"])
    occupancy = {row["family_id"]: row for row in read_table(ont_occupancy_tsv)}
    prior = {
        row["family_id"]: row
        for row in read_table(prior_family_metrics)
        if row["eligibility"] == "eligible"
    }
    rows: list[dict[str, object]] = []
    for family_id in sorted(prior):
        fraction = float(occupancy[family_id]["library_fraction"])
        new_bp = float(prior[family_id]["new_assembly_bp"])
        estimates = {size: fraction * size / efficiency for size in sizes}
        ratios = {
            size: new_bp / estimate if estimate > 0 else math.inf
            for size, estimate in estimates.items()
        }
        is_candidate = family_id in candidates
        if is_candidate and all(value < threshold for value in ratios.values()):
            result = "ONT_direction_supports_residual_collapse"
        elif is_candidate and all(value >= threshold for value in ratios.values()):
            result = "ONT_direction_does_not_support_residual_collapse"
        elif is_candidate:
            result = "ONT_direction_unresolved_across_genome_sizes"
        elif all(value >= threshold for value in ratios.values()):
            result = "negative_context_no_residual_deficit"
        else:
            result = "negative_context_ONT_deficit_or_unresolved"
        row: dict[str, object] = {
            "family_id": family_id,
            "preorthogonal_candidate": str(is_candidate).lower(),
            "accepted_query_bp": occupancy[family_id]["accepted_query_bp"],
            "accepted_read_count": occupancy[family_id]["accepted_read_count"],
            "library_fraction": fraction,
            "mapping_efficiency_correction": efficiency,
            "newer_assembly_bp": new_bp,
        }
        for size in sizes:
            row[f"calibrated_ONT_bp_at_{size}"] = estimates[size]
            row[f"newer_assembly_ONT_ratio_at_{size}"] = ratios[size]
        row["ONT_direction_result"] = result
        row["warning"] = "direction_only;global_HiFi_mapping_calibration;not_copy_number_truth"
        rows.append(row)
    write_table(output_tsv, rows)
    candidate_rows = [row for row in rows if row["preorthogonal_candidate"] == "true"]
    summary: dict[str, object] = {
        "schema_version": 1,
        "status": "complete_direction_only",
        "eligible_families": len(rows),
        "candidate_families": len(candidate_rows),
        "mapping_efficiency_correction": efficiency,
        "genome_size_sensitivity_bp": sizes,
        "candidate_result_counts": {
            state: sum(row["ONT_direction_result"] == state for row in candidate_rows)
            for state in sorted({str(row["ONT_direction_result"]) for row in candidate_rows})
        },
        "inputs": {
            "ont_occupancy_tsv_sha256": digest_file(ont_occupancy_tsv),
            "calibration_json_sha256": digest_file(calibration_json),
            "prior_family_metrics_sha256": digest_file(prior_family_metrics),
            "config_json_sha256": digest_file(config_json),
        },
        "output_tsv_sha256": digest_file(output_tsv),
        "boundary": "ONT supplies secondary direction only; magnitude is estimated under-representation, not missing-bp truth.",
    }
    if summary_json.exists():
        raise ValueError(f"Refusing to overwrite existing output: {summary_json}")
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--catalogue", required=True, type=Path)
    prepare.add_argument("--prior-family-metrics", required=True, type=Path)
    prepare.add_argument("--output-fasta", required=True, type=Path)
    prepare.add_argument("--receipt", required=True, type=Path)
    prepare.add_argument("--template-length", type=int, default=250_000)

    summarize = subparsers.add_parser("summarize")
    summarize.add_argument("--paf", required=True, type=Path)
    summarize.add_argument("--template-fasta", required=True, type=Path)
    summarize.add_argument("--output-tsv", required=True, type=Path)
    summarize.add_argument("--summary", required=True, type=Path)
    summarize.add_argument("--total-library-bases", required=True, type=int)
    summarize.add_argument("--genome-size", required=True, type=int)
    summarize.add_argument("--minimum-block", type=int, default=500)
    summarize.add_argument("--minimum-identity", type=float, default=0.75)

    calibrate = subparsers.add_parser("calibrate")
    calibrate.add_argument("--hifi-occupancy", required=True, type=Path)
    calibrate.add_argument("--prior-family-metrics", required=True, type=Path)
    calibrate.add_argument("--output", required=True, type=Path)

    ont = subparsers.add_parser("evaluate-ont")
    ont.add_argument("--ont-occupancy", required=True, type=Path)
    ont.add_argument("--calibration", required=True, type=Path)
    ont.add_argument("--prior-family-metrics", required=True, type=Path)
    ont.add_argument("--config", required=True, type=Path)
    ont.add_argument("--output-tsv", required=True, type=Path)
    ont.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    if args.command == "prepare":
        prepare_templates(
            args.catalogue,
            args.prior_family_metrics,
            args.output_fasta,
            args.receipt,
            template_length=args.template_length,
        )
    elif args.command == "summarize":
        summarize_paf(
            args.paf,
            args.template_fasta,
            args.output_tsv,
            args.summary,
            total_library_bases=args.total_library_bases,
            genome_size_bp=args.genome_size,
            minimum_block_bp=args.minimum_block,
            minimum_identity=args.minimum_identity,
        )
    elif args.command == "calibrate":
        calibrate_hifi_mapping(
            args.hifi_occupancy, args.prior_family_metrics, args.output
        )
    else:
        evaluate_ont(
            args.ont_occupancy,
            args.calibration,
            args.prior_family_metrics,
            args.config,
            args.output_tsv,
            args.summary,
        )


if __name__ == "__main__":
    main()
