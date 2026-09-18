"""Explain frozen Macadamia M2 abstentions without changing the caller."""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from math import floor
from pathlib import Path

import edlib

from benchmarks.m2_routes.native_read_pilot import fasta_records, fastq_records, reverse_complement
from benchmarks.scripts.build_native_read_collapse import sha256_file


ROOT = Path(__file__).resolve().parents[2]
EDIT = ROOT / "benchmarks/controlled_collapse/macadamia_bp_provisional_v1"
MAPPING = ROOT / "benchmarks/controlled_collapse/macadamia_native_mapping_v1_20260918"
PERIOD = 144
MAX_EDIT = floor(PERIOD * 0.15)
WIGGLE = floor(PERIOD * 0.08)


def primary_overlap_counts(paf: Path, start: int, end: int) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with gzip.open(paf, "rt") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if "tp:A:P" not in fields[12:]:
                continue
            counts["all_primary"] += 1
            left, right = int(fields[7]), int(fields[8])
            if left < end and right > start:
                counts["array_overlap"] += 1
                if left <= start and right >= end:
                    counts["whole_array_span"] += 1
                else:
                    counts["partial_array_overlap"] += 1
            else:
                counts["no_array_overlap"] += 1
    return dict(counts)


def longest_tiled_prefix(sequence: str, motif: str) -> int:
    """Reachability using the unchanged caller's width and edit limits."""
    if not sequence:
        return 0
    motifs = (motif, reverse_complement(motif))
    reachable = bytearray(len(sequence) + 1)
    reachable[0] = 1
    furthest = 0
    for start in range(len(sequence)):
        if not reachable[start]:
            continue
        for unit in motifs:
            for width in range(PERIOD - WIGGLE, PERIOD + WIGGLE + 1):
                end = start + width
                if end > len(sequence):
                    continue
                distance = edlib.align(unit, sequence[start:end], mode="NW",
                                       task="distance", k=MAX_EDIT)["editDistance"]
                if distance >= 0:
                    reachable[end] = 1
                    furthest = max(furthest, end)
    return furthest


def array_diagnostic(sequence: str, motif: str) -> dict:
    if not sequence:
        return {"length_bp": 0, "length_divided_by_144": 0.0,
                "max_tiled_prefix_bp": 0, "full_tiling_possible": True,
                "first_tile_min_edit": None, "first_144bp_best_rotation_edit": None,
                "best_local_template_edit": None,
                "no_full_tiling_interpretation": "empty_interval"}
    orientations = (motif, reverse_complement(motif))
    first = min(edlib.align(unit, sequence[:width], mode="NW", task="distance")["editDistance"]
                for unit in orientations for width in range(PERIOD - WIGGLE, PERIOD + WIGGLE + 1)
                if width <= len(sequence))
    rotated = min(edlib.align(unit[offset:] + unit[:offset], sequence[:PERIOD],
                              mode="NW", task="distance")["editDistance"]
                  for unit in orientations for offset in range(PERIOD)) if len(sequence) >= PERIOD else None
    local = min(edlib.align(unit, sequence, mode="HW", task="distance")["editDistance"]
                for unit in orientations)
    prefix = longest_tiled_prefix(sequence, motif)
    if prefix == len(sequence):
        cause = "full_tiling_possible_by_frozen_emissions"
    elif first > MAX_EDIT and rotated is not None and rotated <= MAX_EDIT:
        cause = "first_boundary_phase_incompatible_with_unrotated_template"
    elif first > MAX_EDIT:
        cause = "first_segment_exceeds_frozen_edit_limit_even_before_path_search"
    else:
        cause = "later_segment_or_terminal_boundary_exceeds_frozen_tiling_constraints"
    return {"length_bp": len(sequence), "length_divided_by_144": len(sequence) / PERIOD,
            "max_tiled_prefix_bp": prefix, "full_tiling_possible": prefix == len(sequence),
            "first_tile_min_edit": first, "first_144bp_best_rotation_edit": rotated,
            "best_local_template_edit": local, "frozen_per_tile_max_edit": MAX_EDIT,
            "frozen_indel_width_range_bp": [PERIOD - WIGGLE, PERIOD + WIGGLE],
            "no_full_tiling_interpretation": cause}


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError("Preserve prior abstention diagnosis")
    m2_protocol = EDIT / "m2_protocol.json"
    m2_summary = EDIT / "m2_score/summary.json"
    m2_rows = EDIT / "m2_score/per_case.jsonl"
    read_path = EDIT / "selected_original_reads.fastq.gz"
    trims_path = EDIT / "span_score/read_trims.json"
    paf = MAPPING / "score/context_alignments.paf.gz"
    template_path = EDIT / "monomer_template.fa"
    protocol = json.loads(m2_protocol.read_text())
    summary = json.loads(m2_summary.read_text())
    if (sha256_file(template_path) != protocol["single_family_monomer_fasta_sha256"] or
            sha256_file(read_path) != protocol["original_selected_fastq_sha256"] or
            sha256_file(paf) != json.loads((EDIT / "score_protocol.json").read_text())[
                "context_paf_gzip_sha256"] or
            sha256_file(m2_rows) != summary["per_case_sha256"]):
        raise ValueError("Frozen Macadamia M2 evidence changed")
    motif = next(iter(fasta_records(template_path).values()))
    if len(motif) != PERIOD:
        raise ValueError("Unexpected operational period")
    reads = fastq_records(read_path)
    trims = json.loads(trims_path.read_text())
    frozen_paths = json.loads((EDIT / "m2_score/read_paths.json").read_text())
    if set(reads) != set(frozen_paths) or len(trims) != 7:
        raise ValueError("Original record denominator changed")
    per_read = []
    for row in trims:
        read_id = row["read_id"]
        sequence = reads[read_id]
        if hashlib.sha256(sequence.encode()).hexdigest() != row["original_read_sequence_sha256"]:
            raise ValueError("Original read changed")
        oriented = sequence if row["strand"] == "+" else reverse_complement(sequence)
        left, right = row["trim_interval_0based"]
        diagnostic = array_diagnostic(oriented[left:right], motif)
        diagnostic.update(read_id=read_id, frozen_m2_state=frozen_paths[read_id]["state"],
                          frozen_m2_reason=frozen_paths[read_id]["reason"],
                          flank_trim_status=row["status"], original_read_length_bp=len(sequence))
        per_read.append(diagnostic)
    frozen_cases = [json.loads(line) for line in m2_rows.read_text().splitlines()]
    edit_receipt = json.loads((EDIT / "generated/receipt.json").read_text())
    if len(frozen_cases) != 9 or len(edit_receipt["cases"]) != 9:
        raise ValueError("Frozen case denominator changed")
    per_case = []
    for frozen, edit in zip(frozen_cases, edit_receipt["cases"]):
        if frozen["case_id"] != edit["case_id"]:
            raise ValueError("Case order changed")
        fasta = EDIT / "generated" / (edit["case_id"] + ".fa")
        if sha256_file(fasta) != edit["edited_fasta_sha256"]:
            raise ValueError("Edited case changed")
        seq = fasta_records(fasta)[edit["case_id"]]
        left, right = edit["edited_array_context_interval"]
        assembly = array_diagnostic(seq[left:right], motif)
        case = {"case_id": frozen["case_id"], "injected_deleted_bp": edit["injected_deleted_bp"],
                "qualified_original_spanning_records": len(per_read),
                "assembly": assembly,
                "assembly_supported_path": frozen["assembly_path"]["labels"] if
                    frozen["assembly_path"]["state"] == "RESOLVED" else None,
                "read_supported_path": None,
                "usable_original_read_decompositions": sum(
                    x["frozen_m2_state"] == "RESOLVED" for x in per_read),
                "read_orientation_consistency": "not_evaluable_no_resolved_read_paths",
                "HOR_candidate_count_under_supplied_one_template": 0,
                "HOR_candidate_scope": "no multi-label HOR search possible with one supplied family template",
                "path_ambiguity": "no_full_monomer_tiling",
                "posterior_probability": None,
                "edit_cost_is_calibrated_confidence": False,
                "frozen_abstention_state": frozen["m2_technical"]["state"],
                "frozen_abstention_rule": frozen["m2_technical"]["reason"],
                "cause_category": "C_monomer_decomposition_failure;D_multi_label_HOR_not_identifiable_from_one_template"}
        per_case.append(case)
    counts = primary_overlap_counts(paf, 3000, 6155)
    result = {"case_count": len(per_case), "qualified_spanning_record_count": len(per_read),
              "all_primary_context_overlap_counts": counts,
              "record_identity_scope": "distinct SRA records; original ZMW unavailable",
              "frozen_M2_prototype_sha256": protocol["frozen_m2_prototype_sha256"],
              "m2_protocol_sha256": sha256_file(m2_protocol),
              "m2_summary_sha256": sha256_file(m2_summary),
              "no_algorithm_or_threshold_change": True,
              "biological_accuracy": "not_evaluated"}
    output.mkdir(parents=True, exist_ok=False)
    (output / "per_read.json").write_text(json.dumps(per_read, indent=2, sort_keys=True) + "\n")
    with (output / "per_case.jsonl").open("w") as stream:
        for row in per_case:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    result["per_read_sha256"] = sha256_file(output / "per_read.json")
    result["per_case_sha256"] = sha256_file(output / "per_case.jsonl")
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
