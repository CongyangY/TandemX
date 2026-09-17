"""Select one rice tandem interval before observing native reads."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def entropy(sequence: bytes) -> float:
    counts = Counter(sequence)
    length = len(sequence)
    return -sum((n / length) * math.log2(n / length) for n in counts.values())


def flank_unique_fraction(context: bytes, array_start: int, array_end: int, k: int = 31) -> float:
    all_kmers = Counter(context[i:i + k] for i in range(len(context) - k + 1))
    positions = list(range(array_start - k + 1)) + list(range(array_end, len(context) - k + 1))
    return sum(all_kmers[context[i:i + k]] == 1 for i in positions) / len(positions)


def fasta_records(path: Path):
    with gzip.open(path, "rb") as stream:
        name = None
        parts = []
        for line in stream:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(parts).upper()
                name = line[1:].split()[0].decode("ascii")
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            yield name, b"".join(parts).upper()


def candidate_scores(sequence: bytes, period: int, length: int, stride: int, flank: int):
    bases = np.frombuffer(sequence, dtype=np.uint8)
    equal = (bases[:-period] == bases[period:]).astype(np.uint8)
    sums = np.empty(len(equal) + 1, dtype=np.int64)
    sums[0] = 0
    np.cumsum(equal, out=sums[1:])
    for start in range(flank, len(sequence) - length - flank + 1, stride):
        matches = int(sums[start + length - period] - sums[start])
        yield start, matches, matches / (length - period)


def select(candidates: list[tuple], sequences: dict[str, bytes], config: dict):
    rule = config["interval_selection_before_read_mapping"]
    flank = rule["flank_bp_each_side"]
    length = rule["array_length_bp"]
    checked = []
    for identity, contig, start, matches in sorted(candidates, key=lambda x: (-x[0], x[1], x[2])):
        context = sequences[contig][start - flank:start + length + flank]
        array = context[flank:flank + length]
        h = entropy(array)
        unique = flank_unique_fraction(context, flank, flank + length)
        passed = (identity >= rule["minimum_period_shift_identity"] and
                  h >= rule["minimum_base_shannon_entropy_bits_per_base"] and
                  unique >= rule["minimum_flank_31mer_uniqueness_fraction_within_context"])
        checked.append((contig, start, start + length, matches, identity, h, unique, passed))
        if passed:
            return checked, context
    return checked, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.protocol.read_text())
    source = Path(config["assembly_path"])
    digest = hashlib.sha256()
    with source.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != config["assembly_sha256_expected"]:
        raise ValueError("assembly SHA-256 mismatch")
    rule = config["interval_selection_before_read_mapping"]
    sequences = dict(fasta_records(source))
    candidates = []
    for contig, sequence in sequences.items():
        for start, matches, identity in candidate_scores(
            sequence, rule["period_bp"], rule["array_length_bp"],
            rule["stride_bp"], rule["flank_bp_each_side"]
        ):
            candidates.append((identity, contig, start, matches))
    checked, context = select(candidates, sequences, config)
    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "ranked_screen.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "start_0", "end_0", "matching_shift_bases", "shift_identity", "array_entropy_bits", "flank_unique_31mer_fraction", "selected"])
        for row in checked:
            writer.writerow([*row[:4], f"{row[4]:.12f}", f"{row[5]:.12f}", f"{row[6]:.12f}", str(row[7]).lower()])
    chosen = checked[-1] if context is not None else None
    if chosen is not None:
        contig, start, end = chosen[:3]
        with (args.outdir / "context.fa").open("w") as stream:
            stream.write(f">{contig}:{start - rule['flank_bp_each_side']}-{end + rule['flank_bp_each_side']} array={start}-{end} period={rule['period_bp']}\n")
            for i in range(0, len(context), 80):
                stream.write(context[i:i + 80].decode("ascii") + "\n")
    receipt = {
        "source_sha256": digest.hexdigest(), "contig_lengths": {k: len(v) for k, v in sequences.items()},
        "candidate_windows": len(candidates), "ranked_candidates_checked": len(checked),
        "selected": None if chosen is None else {"contig": chosen[0], "array_start_0": chosen[1], "array_end_0": chosen[2], "shift_identity": chosen[4], "array_entropy_bits": chosen[5], "flank_unique_31mer_fraction": chosen[6]},
    }
    (args.outdir / "selection_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
