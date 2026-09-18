"""Freeze a Macadamia assembly-only native interval before reading HiFi data.

The selection rule is the already committed 2026-09-17 Macadamia protocol.
FASTA is scanned once; only nominated intervals plus natural flanks are held.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def eligible_intervals(bed_path: Path, families_path: Path, protocol: dict) -> tuple[list[dict], int, list[dict]]:
    family_period = {}
    with families_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            family_period[row["family_id"]] = int(row["monomer_length_bp"])
    rule = protocol["interval_selection_before_read_mapping"]
    low, high = rule["array_length_bp_inclusive"]
    p_low, p_high = rule["family_period_bp_inclusive"]
    flank = rule["natural_flank_bp_each_side"]
    candidates = []
    length_eligible = []
    total = 0
    with bed_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError("BED row has fewer than four columns")
            total += 1
            contig, start_s, end_s, family = fields[:4]
            start, end = int(start_s), int(end_s)
            period = family_period[family]
            length = end - start
            if start < flank or length <= 0:
                continue
            if low <= length <= high:
                length_eligible.append({
                    "contig": contig, "start": start, "end": end,
                    "family_id": family, "length_bp": length,
                    "period_bp": period,
                    "period_eligible": p_low <= period <= p_high,
                })
            if low <= length <= high and p_low <= period <= p_high:
                candidates.append(
                    {"contig": contig, "start": start, "end": end,
                     "family_id": family, "length_bp": length,
                     "period_bp": period, "flank_bp": flank}
                )
    return candidates, total, length_eligible


def extract_contexts(fasta_path: Path, candidates: list[dict]) -> None:
    """Extract nominated [start-flank,end+flank) intervals in one FASTA pass."""
    by_contig: dict[str, list[dict]] = {}
    for candidate in candidates:
        candidate["parts"] = []
        by_contig.setdefault(candidate["contig"], []).append(candidate)
    contig = None
    offset = 0
    seen = set()
    with fasta_path.open("rt") as handle:
        for raw in handle:
            if raw.startswith(">"):
                contig = raw[1:].split()[0]
                if contig in seen:
                    raise ValueError(f"duplicate FASTA contig: {contig}")
                seen.add(contig)
                offset = 0
                continue
            if contig is None:
                raise ValueError("FASTA sequence before header")
            sequence = raw.strip().upper()
            if not sequence:
                continue
            stop = offset + len(sequence)
            for candidate in by_contig.get(contig, []):
                left = candidate["start"] - candidate["flank_bp"]
                right = candidate["end"] + candidate["flank_bp"]
                overlap_start = max(offset, left)
                overlap_end = min(stop, right)
                if overlap_end > overlap_start:
                    candidate["parts"].append(sequence[overlap_start - offset:overlap_end - offset])
            offset = stop
    for candidate in candidates:
        if candidate["contig"] not in seen:
            raise ValueError(f"missing FASTA contig: {candidate['contig']}")
        context = "".join(candidate.pop("parts"))
        expected = candidate["length_bp"] + 2 * candidate["flank_bp"]
        if len(context) != expected:
            raise ValueError(f"truncated FASTA context: {candidate['contig']}:{candidate['start']}-{candidate['end']}")
        candidate["context"] = context


def score_candidates(candidates: list[dict]) -> list[dict]:
    scored = []
    for candidate in candidates:
        flank = candidate["flank_bp"]
        array = candidate["context"][flank:flank + candidate["length_bp"]]
        period = candidate["period_bp"]
        compared = len(array) - period
        if compared <= 0:
            raise ValueError("period must be shorter than array")
        matches = sum(a == b for a, b in zip(array[:-period], array[period:]))
        scored.append({
            key: candidate[key] for key in ("contig", "start", "end", "family_id", "length_bp", "period_bp", "flank_bp")
        } | {
            "period_shift_matches": matches,
            "period_shift_compared": compared,
            "period_shift_identity": matches / compared,
            "array_sha256": sha256_text(array),
            "context_sha256": sha256_text(candidate["context"]),
        })
    return sorted(scored, key=lambda row: (-row["period_shift_identity"], -row["length_bp"], row["contig"], row["start"]))


def run(protocol_path: Path, output_dir: Path) -> dict:
    protocol = json.loads(protocol_path.read_text())
    assembly = Path(protocol["assembly_path"])
    bed = Path(protocol["prior_catalogue"]["arrays_bed"]["path"])
    families = Path(protocol["prior_catalogue"]["families_tsv"]["path"])
    expected = {
        assembly: protocol["assembly_expected_sha256"],
        bed: protocol["prior_catalogue"]["arrays_bed"]["sha256"],
        families: protocol["prior_catalogue"]["families_tsv"]["sha256"],
    }
    actual = {path: sha256(path) for path in expected}
    for path, digest in actual.items():
        if digest != expected[path]:
            raise ValueError(f"source SHA-256 mismatch: {path}")
    candidates, total, length_eligible = eligible_intervals(bed, families, protocol)
    if not candidates:
        raise ValueError("no intervals pass the frozen assembly-only rule")
    extract_contexts(assembly, candidates)
    scored = score_candidates(candidates)
    selected = scored[0]
    selected_raw = next(row for row in candidates if row["contig"] == selected["contig"] and row["start"] == selected["start"] and row["end"] == selected["end"])
    output_dir.mkdir(parents=True, exist_ok=False)
    fields = list(scored[0])
    with (output_dir / "all_candidate_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(scored)
    with (output_dir / "length_filter_trace.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(length_eligible[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(length_eligible)
    context_header = f">{selected['contig']}:{selected['start'] - selected['flank_bp']}-{selected['end'] + selected['flank_bp']} selected_assembly_only\n"
    (output_dir / "selected_context.fa").write_text(context_header + selected_raw["context"] + "\n")
    receipt = {
        "schema_version": 1,
        "selection_type": "assembly_only_before_original_read_inspection",
        "frozen_protocol_path": str(protocol_path),
        "frozen_protocol_sha256": sha256(protocol_path),
        "source_sha256": {str(path): digest for path, digest in actual.items()},
        "catalogue_rows_examined": total,
        "length_eligible_rows": len(length_eligible),
        "eligible_candidates": len(candidates),
        "selection_rank": protocol["interval_selection_before_read_mapping"]["rank"],
        "selected": selected,
        "all_candidate_scores_sha256": sha256(output_dir / "all_candidate_scores.tsv"),
        "length_filter_trace_sha256": sha256(output_dir / "length_filter_trace.tsv"),
        "selected_context_fasta_sha256": sha256(output_dir / "selected_context.fa"),
        "original_reads_inspected": False,
        "biological_truth": "unavailable",
    }
    (output_dir / "selection_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.protocol, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
