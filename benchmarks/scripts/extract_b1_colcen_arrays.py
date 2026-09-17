"""Stream three prespecified Col-CEN v1.2 source intervals from a verified FASTA.gz."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract(source_gz: Path, protocol_path: Path, out_fasta: Path) -> dict:
    protocol = json.loads(protocol_path.read_text())
    if sha256_file(source_gz) != protocol["source_assembly_sha256"]:
        raise ValueError("Col-CEN compressed assembly SHA-256 mismatch")
    intervals = protocol["source_intervals"]
    selected: dict[str, list[str]] = {row["chromosome"]: [] for row in intervals}
    offsets = {chrom: 0 for chrom in selected}
    chrom = None
    with gzip.open(source_gz, "rt", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                if chrom not in selected:
                    continue
                offsets[chrom] = 0
                continue
            if chrom not in selected:
                continue
            sequence = line.strip().upper()
            position = offsets[chrom]
            for row in intervals:
                if row["chromosome"] != chrom:
                    continue
                start, end = row["start0"], row["end0"]
                if position < end and position + len(sequence) > start:
                    selected[chrom].append(sequence[max(0, start - position):min(len(sequence), end - position)])
            offsets[chrom] += len(sequence)
    lines = []
    catalogue = []
    for row in intervals:
        chrom = row["chromosome"]
        sequence = "".join(selected[chrom])
        period = protocol["operational_period_bp"]
        if (len(sequence) != row["end0"] - row["start0"] or len(sequence) != period * protocol["operational_tile_count"]
                or set(sequence) - set("ACGT")):
            raise ValueError(f"source interval missing, wrong length or ambiguous bases: {chrom}")
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest != row["sequence_sha256"]:
            raise ValueError(f"source interval hash mismatch: {chrom}")
        tiles = [sequence[i:i + period] for i in range(0, len(sequence), period)]
        adjacent = [sum(a == b for a, b in zip(left, right)) / period
                    for left, right in zip(tiles, tiles[1:])]
        if min(adjacent) < protocol["minimum_adjacent_tile_identity"]:
            raise ValueError(f"periodicity qualification failed: {chrom}")
        lines.extend((f">{row['array_id']}|{chrom}:{row['start0']}-{row['end0']}|period={period}", sequence))
        catalogue.append(dict(array_id=row["array_id"], chromosome=chrom,
                              start0=row["start0"], end0=row["end0"], period_bp=period,
                              operational_tile_count=len(tiles), sequence_sha256=digest,
                              adjacent_tile_identity_min=min(adjacent),
                              adjacent_tile_identity_mean=sum(adjacent) / len(adjacent),
                              biological_monomer_boundary_status="unverified_operational_phase"))
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_fasta.write_text("\n".join(lines) + "\n", encoding="ascii")
    return dict(source_assembly_sha256=protocol["source_assembly_sha256"],
                source_fasta_sha256=sha256_file(out_fasta), arrays=catalogue)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-gz", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--out-fasta", type=Path, required=True)
    parser.add_argument("--out-catalogue", type=Path, required=True)
    args = parser.parse_args()
    receipt = extract(args.source_gz, args.protocol, args.out_fasta)
    args.out_catalogue.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
