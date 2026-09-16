"""Bounded bp-only editing of a receipt-verified real assembly fragment.

This DEVELOPMENT tool records exact injected base-pair deltas. It deliberately
does not infer complete monomer copy counts or biological copy number.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


LEVELS = (0, 25, 50, 75, 100)
MAX_FRAGMENT_BP = 1_000_000
MAX_FASTA_BYTES = 2_000_000
LEDGER_FIELDS = (
    "case_id", "deletion_percent", "region_id", "source_sequence_id",
    "source_region_start0", "source_region_end0", "fragment_length_bp",
    "deleted_fragment_start0", "deleted_fragment_end0",
    "deleted_genomic_start0", "deleted_genomic_end0", "injected_deleted_bp",
    "retained_bp", "edited_fragment_start0", "edited_fragment_end0",
    "source_sequence_sha256", "edited_sequence_sha256", "period_context_bp",
    "complete_unit_count_truth", "unit_count_status", "truth_scope", "read_pairing_status",
)
CHAIN_FIELDS = ("case_id", "operation", "source_start0", "source_end0", "edited_start0", "edited_end0")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_fragment(source: Path, archive_receipt: Path, archive_manifest: Path,
                  confirmation_summary: Path) -> tuple[str, str, dict, int]:
    """Check the archived extraction, archive manifest, and tandem context."""
    archive = json.loads(archive_receipt.read_text())
    if archive.get("coordinate_system") != "zero_based_half_open" or len(archive.get("extracted", [])) != 1:
        raise ValueError("unexpected archive extraction schema")
    region = archive["extracted"][0]
    if source.stat().st_size > MAX_FASTA_BYTES:
        raise ValueError("fragment FASTA exceeds bounded development input size")
    if sha256(source) != archive["output_sha256"]:
        raise ValueError("archived FASTA SHA-256 mismatch")
    with archive_manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in manifest if row.get("file") == source.name]
    if len(matches) != 1 or matches[0].get("sha256") != archive["output_sha256"]:
        raise ValueError("archive manifest FASTA checksum mismatch")
    lines = source.read_text(encoding="ascii").splitlines()
    if not lines or not lines[0].startswith(">") or any(line.startswith(">") for line in lines[1:]):
        raise ValueError("expected exactly one FASTA record")
    header = lines[0]
    match = re.fullmatch(r">(\S+) sequence_id=([^;]+);start0=(\d+);end0=(\d+)", header)
    if match is None:
        raise ValueError("unexpected archived FASTA header")
    region_id, sequence_id, start, end = match.groups()
    if (region_id != region["region_id"] or sequence_id != region["sequence_id"]
            or int(start) != region["start0"] or int(end) != region["end0"]):
        raise ValueError("FASTA header and extraction receipt disagree")
    sequence = "".join(lines[1:]).upper()
    if (not sequence or len(sequence) > MAX_FRAGMENT_BP or len(sequence) != int(end) - int(start)
            or len(sequence) != region["extracted_length_bp"] or any(base not in "ACGTN" for base in sequence)):
        raise ValueError("invalid fragment length or sequence")
    with confirmation_summary.open(newline="") as handle:
        confirmations = list(csv.DictReader(handle, delimiter="\t"))
    passing = [row for row in confirmations if row.get("verdict") == "confirmed_tandem_structure"
               and row.get("coverage_fraction") == "1.000000"
               and int(row.get("region_length_bp", -1)) == len(sequence)]
    if {row["tool"] for row in passing} != {"TRF", "TideHunter"}:
        raise ValueError("independent tandem-structure context is incomplete")
    periods = {int(row["dominant_period"]) for row in passing}
    if len(periods) != 1:
        raise ValueError("dominant-period context disagrees")
    return region_id, sequence, region, periods.pop()


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate(source: Path, archive_receipt: Path, archive_manifest: Path,
             confirmation_summary: Path, outdir: Path) -> None:
    for path in (source, archive_receipt, archive_manifest, confirmation_summary):
        if not path.is_file():
            raise ValueError(f"input file is absent: {path}")
    region_id, sequence, region, period = read_fragment(source, archive_receipt,
                                                         archive_manifest, confirmation_summary)
    outdir.mkdir(parents=True, exist_ok=False)
    cases = []
    length = len(sequence)
    for level in LEVELS:
        case_id = f"bp_deletion_{level:03d}"
        removed = length * level // 100
        retained = length - removed
        edited = sequence[:retained]
        fasta = outdir / f"{case_id}.fa"
        fasta.write_text(f">{case_id} source_region={region_id};truth=bp_delta_only\n"
                         + (edited + "\n" if edited else ""), encoding="ascii")
        row = dict(case_id=case_id, deletion_percent=level, region_id=region_id,
                   source_sequence_id=region["sequence_id"],
                   source_region_start0=region["start0"], source_region_end0=region["end0"],
                   fragment_length_bp=length, deleted_fragment_start0=retained,
                   deleted_fragment_end0=length,
                   deleted_genomic_start0=region["start0"] + retained,
                   deleted_genomic_end0=region["end0"], injected_deleted_bp=removed,
                   retained_bp=retained, edited_fragment_start0=0,
                   edited_fragment_end0=retained, source_sequence_sha256=sequence_hash(sequence),
                   edited_sequence_sha256=sequence_hash(edited), period_context_bp=period,
                   complete_unit_count_truth="unknown", unit_count_status="unknown_no_copy_boundaries",
                   truth_scope="injected_bp_delta_only", read_pairing_status="not_evaluated")
        ledger = outdir / f"{case_id}.ledger.tsv"
        write_tsv(ledger, LEDGER_FIELDS, [row])
        chain_rows = [dict(case_id=case_id, operation="match", source_start0=0,
                           source_end0=retained, edited_start0=0, edited_end0=retained)] if retained else []
        if removed:
            chain_rows.append(dict(case_id=case_id, operation="delete", source_start0=retained,
                                   source_end0=length, edited_start0=retained, edited_end0=retained))
        chain = outdir / f"{case_id}.chain.tsv"
        write_tsv(chain, CHAIN_FIELDS, chain_rows)
        cases.append(dict(case_id=case_id, fasta_sha256=sha256(fasta),
                          ledger_sha256=sha256(ledger), chain_sha256=sha256(chain),
                          empty_fasta_record=not edited))
    receipt = dict(schema_version=1, development_only=True,
                   source_archive_fasta=str(source.resolve()), source_archive_fasta_sha256=sha256(source),
                   extraction_receipt_sha256=sha256(archive_receipt),
                   archive_manifest_sha256=sha256(archive_manifest),
                   tandem_confirmation_summary_sha256=sha256(confirmation_summary),
                   generator_sha256=sha256(Path(__file__)), source_region=region,
                   source_fragment_length_bp=length, dominant_period_context_bp=period,
                   source_length_mod_period=length % period,
                   complete_unit_count_truth="unknown", truth_scope="injected_bp_delta_only",
                   read_pairing_status="not_evaluated", cases=cases)
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--archive-receipt", type=Path, required=True)
    parser.add_argument("--archive-manifest", type=Path, required=True)
    parser.add_argument("--confirmation-summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    generate(args.source, args.archive_receipt, args.archive_manifest,
             args.confirmation_summary, args.outdir)


if __name__ == "__main__":
    main()
