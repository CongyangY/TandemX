"""Extract independently mappable source reads from locked recovery evidence.

This is mechanical preparation only.  It neither maps reads nor evaluates a
candidate.  Output is published atomically only when every requested ID occurs
exactly once in the enrolled FASTQ and its raw SHA-256 matches the QC receipt.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from benchmarks.scripts.fastq_stream import hashed_fastq, records
from tandemx.recovery.poc import sha256


ID_FIELDS = ("read_id",)


def selected_ids(evidence: Path, locus_id: str, candidate_read_id: str) -> tuple[set[str], int]:
    """Return unique flank-supported IDs, excluding the candidate read itself."""
    with evidence.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus_id", "evidence_type", "read_id"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError("Recruitment evidence lacks required columns")
        matching = [row["read_id"] for row in reader
                    if row["locus_id"] == locus_id and row["evidence_type"] == "flank_anchored"]
    if not matching:
        raise ValueError(f"No flank_anchored reads for {locus_id}")
    if candidate_read_id not in matching:
        raise ValueError(f"Candidate read {candidate_read_id} is absent from {locus_id} flank evidence")
    ids = set(matching)
    ids.discard(candidate_read_id)
    if not ids:
        raise ValueError("Candidate exclusion left no independent read IDs")
    return ids, len(matching)


def _write_receipt(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def extract(
    evidence: Path,
    source_fastq: Path,
    source_qc: Path,
    outdir: Path,
    locus_id: str,
    candidate_read_id: str,
    command: list[str] | None = None,
) -> dict:
    """Stream an enrolled FASTQ into an independent, hash-bound FASTA subset."""
    if outdir.exists():
        raise ValueError(f"Refusing to modify existing output directory: {outdir}")
    if not evidence.is_file() or not source_fastq.is_file() or not source_qc.is_file():
        raise ValueError("Evidence, source FASTQ, and source QC receipt must be files")
    qc = json.loads(source_qc.read_text(encoding="utf-8"))
    if qc.get("complete") is not True or qc.get("fastq_records_valid") is not True or not qc.get("input_sha256"):
        raise ValueError("Source FASTQ does not have a completed valid QC receipt")
    identifiers, evidence_rows = selected_ids(evidence, locus_id, candidate_read_id)
    outdir.mkdir(parents=True)
    fasta = outdir / "reads.fa"
    ids_tsv = outdir / "IDs.tsv"
    fasta_partial = fasta.with_name(fasta.name + ".partial")
    ids_partial = ids_tsv.with_name(ids_tsv.name + ".partial")
    counts = {identifier: 0 for identifier in identifiers}
    source_records = source_bases = extracted_bases = 0
    failure: str | None = None
    try:
        with hashed_fastq(source_fastq) as (handle, source_digest), fasta_partial.open("xb") as output:
            for record in records(handle):
                source_records += 1
                source_bases += len(record.sequence)
                identifier = record.identifier.decode("ascii")
                if identifier in counts:
                    counts[identifier] += 1
                    output.write(b">" + record.identifier + b"\n" + record.sequence + b"\n")
                    extracted_bases += len(record.sequence)
        source_sha256 = source_digest.hexdigest()
        missing = sorted(identifier for identifier, count in counts.items() if count == 0)
        duplicate = sorted(identifier for identifier, count in counts.items() if count > 1)
        if source_sha256 != qc["input_sha256"]:
            raise ValueError("Source FASTQ SHA-256 differs from enrolled QC receipt")
        if missing or duplicate:
            details = []
            if missing:
                details.append(f"missing IDs: {','.join(missing[:5])}" + ("..." if len(missing) > 5 else ""))
            if duplicate:
                details.append(f"duplicate IDs: {','.join(duplicate[:5])}" + ("..." if len(duplicate) > 5 else ""))
            raise ValueError("; ".join(details))
        with ids_partial.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=ID_FIELDS, delimiter="\t")
            writer.writeheader()
            for identifier in sorted(identifiers):
                writer.writerow({"read_id": identifier})
        fasta_partial.replace(fasta)
        ids_partial.replace(ids_tsv)
        result = {
            "complete": True,
            "command": command or [],
            "scope": "mechanical_independent_read_extraction_no_mapping_or_candidate_evaluation",
            "locus_id": locus_id,
            "candidate_read_id_excluded": candidate_read_id,
            "flank_evidence_rows": evidence_rows,
            "distinct_flank_read_ids_before_exclusion": len(identifiers) + 1,
            "extracted_independent_read_ids": len(identifiers),
            "extracted_bases": extracted_bases,
            "source_fastq": str(source_fastq.resolve()),
            "source_fastq_sha256": source_sha256,
            "source_qc": str(source_qc.resolve()),
            "source_qc_sha256": sha256(source_qc),
            "recruited_reads_tsv": str(evidence.resolve()),
            "recruited_reads_tsv_sha256": sha256(evidence),
            "source_records_streamed": source_records,
            "source_bases_streamed": source_bases,
            "outputs": {"reads.fa": sha256(fasta), "IDs.tsv": sha256(ids_tsv)},
        }
    except Exception as exc:
        failure = str(exc)
        result = {
            "complete": False, "locus_id": locus_id, "candidate_read_id_excluded": candidate_read_id,
            "flank_evidence_rows": evidence_rows, "requested_independent_read_ids": len(identifiers),
            "source_fastq": str(source_fastq.resolve()), "source_qc": str(source_qc.resolve()),
            "failure": failure,
        }
        _write_receipt(outdir / "receipt.json", result)
        raise
    _write_receipt(outdir / "receipt.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recruited-reads", type=Path, required=True)
    parser.add_argument("--source-fastq", type=Path, required=True)
    parser.add_argument("--source-qc", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--locus-id", required=True)
    parser.add_argument("--exclude-candidate-read", required=True)
    args = parser.parse_args()
    import sys
    command = [sys.executable, str(Path(__file__).resolve()),
               "--recruited-reads", str(args.recruited_reads.resolve()),
               "--source-fastq", str(args.source_fastq.resolve()),
               "--source-qc", str(args.source_qc.resolve()),
               "--outdir", str(args.outdir.resolve()), "--locus-id", args.locus_id,
               "--exclude-candidate-read", args.exclude_candidate_read]
    extract(args.recruited_reads, args.source_fastq, args.source_qc, args.outdir,
            args.locus_id, args.exclude_candidate_read, command)


if __name__ == "__main__":
    main()
