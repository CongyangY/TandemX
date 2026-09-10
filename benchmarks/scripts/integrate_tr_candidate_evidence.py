#!/usr/bin/env python3
"""Integrate read, assembly, annotation, and known-library evidence for TR families.

This is a post hoc candidate-triage utility.  A shortlisted row means that it
passes the explicitly supplied follow-up rules; it is not proof that the repeat
is novel.  Coordinates emitted by this script are BED-style, zero-based and
half-open.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import zlib


OUTPUT_FIELDS = (
    "family_id",
    "monomer_length_bp",
    "support_read_count",
    "support_span_bp",
    "mean_identity",
    "confidence",
    "low_complexity_flag",
    "gc_fraction",
    "shannon_entropy_bits",
    "normalized_base_entropy",
    "max_homopolymer_bp",
    "best_short_period_bp",
    "best_short_period_match_fraction",
    "zlib_ratio",
    "hierarchy_degree",
    "redundancy_degree",
    "known_exclusion_state",
    "best_known_id",
    "best_known_identity",
    "array_count",
    "merged_array_count",
    "total_array_bp",
    "selected_seqid",
    "selected_start0",
    "selected_end0",
    "selected_span_bp",
    "selected_span_monomer_ratio",
    "selected_array_score",
    "selected_distance_to_nearest_end_bp",
    "selected_overlaps_primary_centromere",
    "selected_overlaps_quartet_centromere",
    "selected_overlaps_any_centromere_definition",
    "selected_gene_overlap_count",
    "nearest_upstream_gene_id",
    "nearest_upstream_gene_distance_bp",
    "nearest_downstream_gene_id",
    "nearest_downstream_gene_distance_bp",
    "passes_confidence_gate",
    "passes_complexity_gate",
    "passes_relation_gate",
    "passes_known_library_gate",
    "passes_read_span_gate",
    "passes_array_span_gate",
    "passes_centromere_exclusion_gate",
    "passes_terminal_exclusion_gate",
    "passes_gene_overlap_gate",
    "shortlist_pass",
    "candidate_state",
    "interpretation_boundary",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                header_token = line[1:].split()[0]
                structured = dict(
                    field.split("=", 1)
                    for field in header_token.split(";")
                    if "=" in field
                )
                current = structured.get("family_id", header_token)
                if not current or current in records:
                    raise ValueError(f"invalid or duplicate FASTA identifier: {current!r}")
                records[current] = []
            elif current is None:
                raise ValueError("FASTA sequence occurs before its header")
            else:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def shannon_entropy(sequence: str) -> float:
    counts = Counter(base for base in sequence.upper() if base in "ACGT")
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def max_homopolymer(sequence: str) -> int:
    best = current = 0
    previous = ""
    for base in sequence.upper():
        if base == previous and base in "ACGT":
            current += 1
        elif base in "ACGT":
            current = 1
        else:
            current = 0
        previous = base
        best = max(best, current)
    return best


def best_short_period(sequence: str, maximum: int = 50) -> tuple[int, float]:
    upper = sequence.upper()
    if len(upper) < 2:
        return 0, 0.0
    best_period, best_fraction = 1, -1.0
    for period in range(1, min(maximum, len(upper) - 1) + 1):
        comparable = len(upper) - period
        matches = sum(upper[index] == upper[index + period] for index in range(comparable))
        fraction = matches / comparable
        if fraction > best_fraction:
            best_period, best_fraction = period, fraction
    return best_period, best_fraction


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"expected true/false, observed {value!r}")
    return normalized == "true"


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def merge_intervals(rows: list[dict[str, object]], maximum_gap: int) -> list[dict[str, object]]:
    merged: list[dict[str, object]] = []
    ordered = sorted(rows, key=lambda row: (str(row["seqid"]), int(row["start0"]), int(row["end0"])))
    for row in ordered:
        if (
            merged
            and row["seqid"] == merged[-1]["seqid"]
            and int(row["start0"]) <= int(merged[-1]["end0"]) + maximum_gap
        ):
            merged[-1]["end0"] = max(int(merged[-1]["end0"]), int(row["end0"]))
            merged[-1]["score"] = max(float(merged[-1]["score"]), float(row["score"]))
            merged[-1]["part_count"] = int(merged[-1]["part_count"]) + 1
        else:
            merged.append(
                {
                    "seqid": row["seqid"],
                    "start0": int(row["start0"]),
                    "end0": int(row["end0"]),
                    "score": float(row["score"]),
                    "part_count": 1,
                }
            )
    return merged


def parse_gff_genes(path: Path, aliases: dict[str, str]) -> dict[str, list[dict[str, object]]]:
    genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            seqid = aliases.get(fields[0], fields[0])
            match = re.search(r"(?:^|;)ID=([^;]+)", fields[8])
            gene_id = match.group(1) if match else "--"
            genes[seqid].append(
                {"start0": int(fields[3]) - 1, "end0": int(fields[4]), "gene_id": gene_id}
            )
    for rows in genes.values():
        rows.sort(key=lambda row: (int(row["start0"]), int(row["end0"])))
    return genes


def gene_context(
    seqid: str, start0: int, end0: int, genes: dict[str, list[dict[str, object]]]
) -> tuple[int, str, int | None, str, int | None]:
    overlaps = 0
    upstream: tuple[int, str] | None = None
    downstream: tuple[int, str] | None = None
    for gene in genes.get(seqid, []):
        gene_start, gene_end = int(gene["start0"]), int(gene["end0"])
        gene_id = str(gene["gene_id"])
        if interval_overlap(start0, end0, gene_start, gene_end):
            overlaps += 1
        elif gene_end <= start0:
            distance = start0 - gene_end
            if upstream is None or distance < upstream[0]:
                upstream = (distance, gene_id)
        elif gene_start >= end0:
            distance = gene_start - end0
            if downstream is None or distance < downstream[0]:
                downstream = (distance, gene_id)
    return (
        overlaps,
        upstream[1] if upstream else "--",
        upstream[0] if upstream else None,
        downstream[1] if downstream else "--",
        downstream[0] if downstream else None,
    )


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def gate_columns(gates: dict[str, bool]) -> dict[str, str]:
    columns = {f"passes_{name}_gate": str(value).lower() for name, value in gates.items()}
    unexpected = sorted(set(columns).difference(OUTPUT_FIELDS))
    if unexpected:
        raise ValueError(f"gate output columns are not declared: {', '.join(unexpected)}")
    return columns


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", required=True, type=Path)
    parser.add_argument("--monomers", required=True, type=Path)
    parser.add_argument("--hierarchy", required=True, type=Path)
    parser.add_argument("--similarity", required=True, type=Path)
    parser.add_argument("--known-summary", required=True, type=Path)
    parser.add_argument("--arrays", required=True, type=Path)
    parser.add_argument("--sequence-lengths", required=True, type=Path)
    parser.add_argument("--sequence-aliases", required=True, type=Path)
    parser.add_argument("--centromeres", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--minimum-confidence", choices=("high",), default="high")
    parser.add_argument("--minimum-read-span", type=int, default=50_000)
    parser.add_argument("--minimum-array-monomer-ratio", type=float, default=10.0)
    parser.add_argument("--terminal-exclusion-bp", type=int, default=100_000)
    parser.add_argument("--maximum-array-merge-gap", type=int, default=100)
    args = parser.parse_args()

    inputs = [
        args.families,
        args.monomers,
        args.hierarchy,
        args.similarity,
        args.known_summary,
        args.arrays,
        args.sequence_lengths,
        args.sequence_aliases,
        args.centromeres,
        args.gff,
    ]
    for path in inputs:
        if not path.is_file():
            parser.error(f"input does not exist: {path}")

    families = read_tsv(args.families)
    monomers = read_fasta(args.monomers)
    known_rows = {row["family_id"]: row for row in read_tsv(args.known_summary)}
    lengths = {row["sequence_id"]: int(row["length_bp"]) for row in read_tsv(args.sequence_lengths)}
    aliases = {row["annotation_seqid"]: row["assembly_seqid"] for row in read_tsv(args.sequence_aliases)}
    genes = parse_gff_genes(args.gff, aliases)

    hierarchy_degree: Counter[str] = Counter()
    for row in read_tsv(args.hierarchy):
        hierarchy_degree[row["shorter_family_id"]] += 1
        hierarchy_degree[row["longer_family_id"]] += 1
    redundancy_degree: Counter[str] = Counter()
    for row in read_tsv(args.similarity):
        if parse_bool(row["redundant_candidate"]):
            redundancy_degree[row["family_a"]] += 1
            redundancy_degree[row["family_b"]] += 1

    arrays: dict[str, list[dict[str, object]]] = defaultdict(list)
    with args.arrays.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{args.arrays}:{line_number}: expected at least 5 BED fields")
            arrays[fields[3]].append(
                {"seqid": fields[0], "start0": int(fields[1]), "end0": int(fields[2]), "score": float(fields[4])}
            )

    centromeres: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for row in read_tsv(args.centromeres):
        centromeres[row["definition"]][row["assembly_seqid"]].append(
            (int(row["start0"]), int(row["end0"]))
        )

    output_rows: list[dict[str, object]] = []
    shortlisted_sequences: dict[str, str] = {}
    for family in families:
        family_id = family["family_id"]
        if family_id not in monomers or family_id not in known_rows:
            raise ValueError(f"missing monomer or known-library summary for {family_id}")
        sequence = monomers[family_id]
        monomer_length = int(family["monomer_length_bp"])
        if len(sequence) != monomer_length:
            raise ValueError(f"length mismatch for {family_id}: FASTA={len(sequence)} TSV={monomer_length}")
        raw_arrays = arrays.get(family_id, [])
        merged = merge_intervals(raw_arrays, args.maximum_array_merge_gap)
        selected = max(
            merged,
            key=lambda row: (int(row["end0"]) - int(row["start0"]), float(row["score"])),
            default=None,
        )
        known = known_rows[family_id]
        entropy = shannon_entropy(sequence)
        period, period_fraction = best_short_period(sequence)
        if selected:
            seqid, start0, end0 = str(selected["seqid"]), int(selected["start0"]), int(selected["end0"])
            span = end0 - start0
            if seqid not in lengths:
                raise ValueError(f"assembly sequence is missing from length table: {seqid}")
            distance_to_end = min(start0, lengths[seqid] - end0)
            primary_overlap = any(interval_overlap(start0, end0, left, right) for left, right in centromeres.get("primary", {}).get(seqid, []))
            quartet_overlap = any(interval_overlap(start0, end0, left, right) for left, right in centromeres.get("QuarTeT", {}).get(seqid, []))
            any_centromere_overlap = any(
                interval_overlap(start0, end0, left, right)
                for by_seqid in centromeres.values()
                for left, right in by_seqid.get(seqid, [])
            )
            overlap_count, upstream_id, upstream_distance, downstream_id, downstream_distance = gene_context(
                seqid, start0, end0, genes
            )
        else:
            seqid, start0, end0, span, distance_to_end = "--", None, None, 0, None
            primary_overlap = quartet_overlap = any_centromere_overlap = False
            overlap_count, upstream_id, upstream_distance, downstream_id, downstream_distance = 0, "--", None, "--", None

        gates = {
            "confidence": family["confidence"] == args.minimum_confidence,
            "complexity": not parse_bool(family["low_complexity_flag"]),
            "relation": hierarchy_degree[family_id] == 0 and redundancy_degree[family_id] == 0,
            "known_library": known["exclusion_state"] == "no_match_in_limited_library",
            "read_span": int(family["support_span_bp"]) >= args.minimum_read_span,
            "array_span": bool(selected) and span / monomer_length >= args.minimum_array_monomer_ratio,
            "centromere_exclusion": bool(selected) and not any_centromere_overlap,
            "terminal_exclusion": bool(selected) and distance_to_end is not None and distance_to_end >= args.terminal_exclusion_bp,
            "gene_overlap": bool(selected) and overlap_count == 0,
        }
        shortlist = all(gates.values())
        if shortlist:
            shortlisted_sequences[family_id] = sequence
        output_rows.append(
            {
                "family_id": family_id,
                "monomer_length_bp": monomer_length,
                "support_read_count": family["support_read_count"],
                "support_span_bp": family["support_span_bp"],
                "mean_identity": family["mean_identity"],
                "confidence": family["confidence"],
                "low_complexity_flag": family["low_complexity_flag"],
                "gc_fraction": family["gc_fraction"],
                "shannon_entropy_bits": f"{entropy:.6f}",
                "normalized_base_entropy": f"{entropy / 2:.6f}",
                "max_homopolymer_bp": max_homopolymer(sequence),
                "best_short_period_bp": period,
                "best_short_period_match_fraction": f"{period_fraction:.6f}",
                "zlib_ratio": f"{len(zlib.compress(sequence.encode('ascii'), 9)) / len(sequence):.6f}",
                "hierarchy_degree": hierarchy_degree[family_id],
                "redundancy_degree": redundancy_degree[family_id],
                "known_exclusion_state": known["exclusion_state"],
                "best_known_id": known["best_known_id"],
                "best_known_identity": known["glocal_edit_identity"],
                "array_count": len(raw_arrays),
                "merged_array_count": len(merged),
                "total_array_bp": sum(int(row["end0"]) - int(row["start0"]) for row in raw_arrays),
                "selected_seqid": seqid,
                "selected_start0": "--" if start0 is None else start0,
                "selected_end0": "--" if end0 is None else end0,
                "selected_span_bp": span,
                "selected_span_monomer_ratio": f"{span / monomer_length:.6f}",
                "selected_array_score": "--" if selected is None else f"{float(selected['score']):.6f}",
                "selected_distance_to_nearest_end_bp": "--" if distance_to_end is None else distance_to_end,
                "selected_overlaps_primary_centromere": str(primary_overlap).lower(),
                "selected_overlaps_quartet_centromere": str(quartet_overlap).lower(),
                "selected_overlaps_any_centromere_definition": str(any_centromere_overlap).lower(),
                "selected_gene_overlap_count": overlap_count,
                "nearest_upstream_gene_id": upstream_id,
                "nearest_upstream_gene_distance_bp": "--" if upstream_distance is None else upstream_distance,
                "nearest_downstream_gene_id": downstream_id,
                "nearest_downstream_gene_distance_bp": "--" if downstream_distance is None else downstream_distance,
                **gate_columns(gates),
                "shortlist_pass": str(shortlist).lower(),
                "candidate_state": "preliminary_previously_unreported_candidate" if shortlist else "not_shortlisted",
                "interpretation_boundary": "post hoc triage; limited known library; depth and orthogonal validation pending; not a novelty claim",
            }
        )

    output_rows.sort(key=lambda row: (row["shortlist_pass"] != "true", -int(row["selected_span_bp"]), str(row["family_id"])))
    args.outdir.mkdir(parents=True, exist_ok=True)
    evidence_path = args.outdir / "family_evidence.tsv"
    shortlist_path = args.outdir / "preliminary_shortlist.tsv"
    fasta_path = args.outdir / "preliminary_shortlist.fa"
    write_tsv(evidence_path, output_rows)
    write_tsv(shortlist_path, [row for row in output_rows if row["shortlist_pass"] == "true"])
    with fasta_path.open("w", encoding="utf-8") as handle:
        for family_id in sorted(shortlisted_sequences):
            sequence = shortlisted_sequences[family_id]
            handle.write(f">{family_id} preliminary_candidate_not_novel\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")

    receipt = {
        "schema_version": 1,
        "analysis_role": "post_hoc_candidate_triage",
        "family_count": len(output_rows),
        "shortlist_count": len(shortlisted_sequences),
        "shortlist_family_ids": sorted(shortlisted_sequences),
        "thresholds": {
            "minimum_confidence": args.minimum_confidence,
            "reject_tandemx_low_complexity": True,
            "maximum_hierarchy_degree": 0,
            "maximum_redundancy_degree": 0,
            "required_known_exclusion_state": "no_match_in_limited_library",
            "minimum_read_span_bp": args.minimum_read_span,
            "minimum_selected_array_monomer_ratio": args.minimum_array_monomer_ratio,
            "terminal_exclusion_bp": args.terminal_exclusion_bp,
            "maximum_array_merge_gap_bp": args.maximum_array_merge_gap,
            "reject_any_published_centromere_overlap": True,
            "maximum_gene_overlap_count": 0,
        },
        "input_sha256": {path.name: sha256(path) for path in inputs},
        "output_sha256": {
            evidence_path.name: sha256(evidence_path),
            shortlist_path.name: sha256(shortlist_path),
            fasta_path.name: sha256(fasta_path),
        },
        "interpretation_boundary": "shortlisted means follow-up candidate only; no-match is limited to the fixed known library; novelty is unresolved",
    }
    (args.outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
