"""Audit the six selected Col-CEN reads against a complete reference PAF.

This checks mapper competition across the supplied full assembly. It cannot
establish exact donor identity, biological copy number, or exhaustive mapping.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/controlled_collapse/native_pair_audit_20260917"
FIELDS = ("read_id", "local_context", "target", "target_start0", "target_end0",
          "strand", "mapq", "identity", "left_flank_bp", "right_flank_bp",
          "reported_secondary_count", "status")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_paf_line(line: str) -> dict:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError("PAF line has fewer than 12 fields")
    query_length, query_start, query_end = map(int, fields[1:4])
    target_length, target_start, target_end = map(int, fields[6:9])
    matches, alignment_columns, mapq = map(int, fields[9:12])
    if (not 0 <= query_start < query_end <= query_length
            or not 0 <= target_start < target_end <= target_length
            or not 0 <= matches <= alignment_columns or alignment_columns < 1
            or not 0 <= mapq <= 255 or fields[4] not in {"+", "-"}):
        raise ValueError("Invalid PAF coordinates, identity, or MAPQ")
    tags = fields[12:]
    status_tags = [tag for tag in tags if tag.startswith("tp:A:")]
    if len(status_tags) != 1 or status_tags[0] not in {"tp:A:P", "tp:A:S"}:
        raise ValueError("Expected one explicit PAF primary/secondary tag")
    return {"read_id": fields[0], "query_length": query_length,
            "query_start0": query_start, "query_end0": query_end,
            "strand": fields[4], "target": fields[5],
            "target_start0": target_start, "target_end0": target_end,
            "matches": matches, "alignment_columns": alignment_columns,
            "identity": matches / alignment_columns, "mapq": mapq,
            "primary": status_tags[0] == "tp:A:P"}


def score(protocol_path: Path, paf_path: Path, outdir: Path) -> dict:
    protocol = json.loads(protocol_path.read_text())
    if protocol.get("status") != "frozen_before_full_reference_mapping":
        raise ValueError("Mapping protocol is not frozen")
    for label, path_key, hash_key in (
        ("reference gzip", "reference_compressed_path", "reference_compressed_sha256"),
        ("reference FASTA", "reference_decompressed_path", "reference_decompressed_sha256"),
        ("native reads", "reads_path", "reads_sha256"),
        ("minimap2", "minimap2_executable", "minimap2_executable_sha256"),
    ):
        if sha256_file(Path(protocol[path_key])) != protocol[hash_key]:
            raise ValueError(f"Frozen {label} hash mismatch")
    source = json.loads((SOURCE / "source_and_alignment_receipt.json").read_text())
    expected = {row["read_id"]: row for row in source["native_molecules"]}
    if len(expected) != 6 or protocol["reads_sha256"] != source["native_spanners_fastq_sha256"]:
        raise ValueError("Selected native-read set changed")
    alignments: dict[str, list[dict]] = {read_id: [] for read_id in expected}
    with paf_path.open() as stream:
        for line in stream:
            if not line.strip():
                continue
            row = parse_paf_line(line)
            if row["read_id"] not in expected:
                raise ValueError("PAF contains an undeclared read")
            alignments[row["read_id"]].append(row)
    rows = []
    for read_id, source_row in sorted(expected.items()):
        choices = alignments[read_id]
        primaries = [row for row in choices if row["primary"]]
        if len(primaries) != 1:
            raise ValueError(f"Expected exactly one primary alignment for {read_id}")
        best = primaries[0]
        contig, start, end = protocol["target_arrays_zero_based_half_open"][source_row["context"]]
        left = start - best["target_start0"]
        right = best["target_end0"] - end
        passing = (best["target"] == contig and left >= 1000 and right >= 1000
                   and best["identity"] >= 0.99 and best["mapq"] >= 20)
        rows.append({"read_id": read_id, "local_context": source_row["context"],
                     "target": best["target"], "target_start0": best["target_start0"],
                     "target_end0": best["target_end0"], "strand": best["strand"],
                     "mapq": best["mapq"], "identity": round(best["identity"], 8),
                     "left_flank_bp": left, "right_flank_bp": right,
                     "reported_secondary_count": len(choices) - 1,
                     "status": "full_reference_primary_spanner" if passing else "not_eligible"})
    if outdir.exists():
        raise FileExistsError("Preserve existing full-reference audit evidence")
    outdir.mkdir(parents=True)
    with (outdir / "per_read.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {"selected_read_count": len(rows),
               "full_reference_primary_spanner_count": sum(
                   row["status"] == "full_reference_primary_spanner" for row in rows),
               "by_local_context": {context: sum(row["local_context"] == context and
                   row["status"] == "full_reference_primary_spanner" for row in rows)
                   for context in ("C1", "C2", "C3")},
               "total_reported_secondary_count": sum(row["reported_secondary_count"] for row in rows),
               "minimum_primary_mapq": min(row["mapq"] for row in rows),
               "minimum_primary_identity": min(row["identity"] for row in rows),
               "interpretation": "mapper_dependent_full_reference_competition_for_locally_selected_reads_only",
               "donor_haplotype_pairing": "unverified",
               "native_biological_copy_truth": "unverified"}
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    receipt = {"protocol_sha256": sha256_file(protocol_path), "paf_sha256": sha256_file(paf_path),
               "per_read_sha256": sha256_file(outdir / "per_read.tsv"),
               "summary_sha256": sha256_file(outdir / "summary.json"),
               "scorer_sha256": sha256_file(Path(__file__)),
               "source_receipt_sha256": sha256_file(SOURCE / "source_and_alignment_receipt.json")}
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(score(args.protocol, args.paf, args.outdir), sort_keys=True))


if __name__ == "__main__":
    main()
