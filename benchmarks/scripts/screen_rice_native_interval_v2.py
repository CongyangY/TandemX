"""Pre-read rice v2 interval selection from frozen assembly-only v1 seeds."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from benchmarks.scripts.screen_rice_native_interval import (
    entropy,
    fasta_records,
    flank_unique_fraction,
)


def seeds_from_v1(path: Path, expected_hash: str, count: int) -> list[dict[str, str]]:
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_hash:
        raise ValueError("v1 assembly-derived seed table SHA-256 mismatch")
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 1000 or any(row["selected"] != "false" for row in rows):
        raise ValueError("v1 seed table shape or negative outcome changed")
    return sorted(rows, key=lambda r: (-float(r["flank_unique_31mer_fraction"]),
                                       -float(r["shift_identity"]), r["contig"],
                                       int(r["start_0"])))[:count]


def score_periods(array: bytes, first: int, last: int):
    bases = np.frombuffer(array, dtype=np.uint8)
    for period in range(first, last + 1):
        matches = int(np.count_nonzero(bases[:-period] == bases[period:]))
        yield period, matches, matches / (len(array) - period)


def rank_key(row: tuple) -> tuple:
    contig, start, end, period, matches, identity, h, unique, eligible = row
    return (-unique, -identity, -(end - start), contig, start, period)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.protocol.read_text())
    prior = config["prior_assembly_only_candidate_table"]
    rule = config["candidate_generation_before_any_read_alignment"]
    seeds = seeds_from_v1(Path(prior["path"]), prior["sha256"], prior["seed_count"])
    source = Path(config["assembly_path"])
    digest = hashlib.sha256()
    with source.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != config["assembly_sha256"]:
        raise ValueError("assembly SHA-256 mismatch")
    sequences = dict(fasta_records(source))
    flank = rule["natural_flank_bp_each_side"]
    first, last = rule["periods_bp_inclusive"]
    rows: list[tuple] = []
    region_rows = []
    for seed in seeds:
        contig = seed["contig"]
        center = (int(seed["start_0"]) + int(seed["end_0"])) // 2
        sequence = sequences[contig]
        for length in rule["array_lengths_bp"]:
            start = center - length // 2
            end = start + length
            if start < flank or end + flank > len(sequence):
                raise ValueError("fixed seed lacks natural flanks; no reselect permitted")
            context = sequence[start - flank:end + flank]
            array = context[flank:flank + length]
            h = entropy(array)
            unique = flank_unique_fraction(context, flank, flank + length)
            region_rows.append((contig, start, end, h, unique))
            for period, matches, identity in score_periods(array, first, last):
                eligible = (identity >= rule["eligibility_filter"]["direct_period_shift_identity_min"]
                            and h >= rule["eligibility_filter"]["array_entropy_bits_min"])
                rows.append((contig, start, end, period, matches, identity, h, unique, eligible))
    if len(rows) != rule["expected_evaluations"]:
        raise ValueError("candidate evaluation denominator mismatch")
    eligible_rows = sorted((row for row in rows if row[-1]), key=rank_key)
    chosen = eligible_rows[0] if eligible_rows else None
    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "all_candidate_scores.tsv.gz").open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="", write_through=True) as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(["contig", "start_0", "end_0", "period_bp", "matching_shift_bases", "shift_identity", "array_entropy_bits", "flank_unique_31mer_fraction", "eligible"])
                for row in rows:
                    writer.writerow([*row[:5], f"{row[5]:.12f}", f"{row[6]:.12f}", f"{row[7]:.12f}", str(row[8]).lower()])
    with (args.outdir / "region_summary.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "start_0", "end_0", "array_entropy_bits", "flank_unique_31mer_fraction"])
        for row in region_rows:
            writer.writerow([*row[:3], f"{row[3]:.12f}", f"{row[4]:.12f}"])
    if chosen is not None:
        contig, start, end, period = chosen[:4]
        context = sequences[contig][start - flank:end + flank]
        with (args.outdir / "selected_context.fa").open("w") as stream:
            stream.write(f">{contig}:{start - flank}-{end + flank} array={start}-{end} period={period}\n")
            for i in range(0, len(context), 80):
                stream.write(context[i:i + 80].decode("ascii") + "\n")
    receipt = {
        "source_sha256": digest.hexdigest(),
        "v1_seed_table_sha256": prior["sha256"], "seed_count": len(seeds),
        "region_count": len(region_rows), "period_evaluations": len(rows),
        "eligible_evaluations": len(eligible_rows),
        "selected": None if chosen is None else {
            "contig": chosen[0], "array_start_0": chosen[1], "array_end_0": chosen[2],
            "period_bp": chosen[3], "matching_shift_bases": chosen[4],
            "shift_identity": chosen[5], "array_entropy_bits": chosen[6],
            "flank_unique_31mer_fraction": chosen[7],
        },
    }
    for filename in ("all_candidate_scores.tsv.gz", "region_summary.tsv", "selected_context.fa"):
        path = args.outdir / filename
        if path.exists():
            receipt[filename + "_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (args.outdir / "selection_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
