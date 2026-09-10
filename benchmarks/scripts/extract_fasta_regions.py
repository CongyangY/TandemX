#!/usr/bin/env python3
"""Stream a FASTA/FASTA.GZ once and extract fixed BED-style regions."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import TextIO


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def extract_regions(fasta: Path, regions: list[dict[str, object]]) -> dict[str, str]:
    by_seqid: dict[str, list[dict[str, object]]] = {}
    for region in regions:
        by_seqid.setdefault(str(region["sequence_id"]), []).append(region)
    pieces = {str(region["region_id"]): [] for region in regions}
    current_seqid: str | None = None
    sequence_offset = 0
    with open_text(fasta) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seqid = line[1:].split()[0]
                sequence_offset = 0
                continue
            if current_seqid is None:
                raise ValueError("FASTA sequence occurs before its header")
            chunk_start = sequence_offset
            chunk_end = chunk_start + len(line)
            for region in by_seqid.get(current_seqid, []):
                start0, end0 = int(region["start0"]), int(region["end0"])
                overlap_start = max(start0, chunk_start)
                overlap_end = min(end0, chunk_end)
                if overlap_start < overlap_end:
                    pieces[str(region["region_id"])].append(
                        line[overlap_start - chunk_start : overlap_end - chunk_start].upper()
                    )
            sequence_offset = chunk_end
    extracted = {region_id: "".join(parts) for region_id, parts in pieces.items()}
    for region in regions:
        region_id = str(region["region_id"])
        expected = int(region["end0"]) - int(region["start0"])
        if len(extracted[region_id]) != expected:
            raise ValueError(
                f"region {region_id} expected {expected} bp but extracted {len(extracted[region_id])} bp"
            )
    return extracted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    if not args.fasta.is_file() or not args.regions.is_file():
        parser.error("FASTA and region table must exist")
    with args.regions.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"region_id", "sequence_id", "start0", "end0"}
    if not rows or not required.issubset(rows[0]):
        parser.error("region table must contain region_id, sequence_id, start0, and end0")
    region_ids: set[str] = set()
    regions: list[dict[str, object]] = []
    for row in rows:
        region_id = row["region_id"]
        start0, end0 = int(row["start0"]), int(row["end0"])
        if region_id in region_ids or start0 < 0 or end0 <= start0:
            parser.error(f"invalid or duplicate region: {region_id}")
        region_ids.add(region_id)
        regions.append({"region_id": region_id, "sequence_id": row["sequence_id"], "start0": start0, "end0": end0})
    sequences = extract_regions(args.fasta, regions)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for region in regions:
            region_id = str(region["region_id"])
            handle.write(
                f">{region_id} sequence_id={region['sequence_id']};start0={region['start0']};end0={region['end0']}\n"
            )
            sequence = sequences[region_id]
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")
    receipt = {
        "schema_version": 1,
        "coordinate_system": "zero_based_half_open",
        "fasta": str(args.fasta.resolve()),
        "fasta_sha256": sha256(args.fasta),
        "regions": str(args.regions.resolve()),
        "regions_sha256": sha256(args.regions),
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
        "extracted": [
            {**region, "extracted_length_bp": len(sequences[str(region["region_id"])])}
            for region in regions
        ],
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
