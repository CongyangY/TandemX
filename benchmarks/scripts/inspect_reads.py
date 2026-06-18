#!/usr/bin/env python3
"""Inspect a FASTA/FASTQ file with TandemX's streaming sequence reader."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from tandemx.io.sequences import detect_sequence_format, read_sequence_records


FIELDS = [
    "input_file",
    "read_count",
    "total_bases",
    "mean_read_length",
    "n50_read_length",
    "min_read_length",
    "max_read_length",
    "format",
    "notes",
]


def calculate_n50(lengths: list[int], total_bases: int) -> int:
    cumulative = 0
    threshold = (total_bases + 1) // 2
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def inspect_reads(path: Path) -> dict[str, str | int]:
    lengths: list[int] = []
    total_bases = 0
    for record in read_sequence_records(path):
        length = len(record.sequence)
        lengths.append(length)
        total_bases += length

    read_count = len(lengths)
    return {
        "input_file": str(path),
        "read_count": read_count,
        "total_bases": total_bases,
        "mean_read_length": f"{total_bases / read_count:.3f}",
        "n50_read_length": calculate_n50(lengths, total_bases),
        "min_read_length": min(lengths),
        "max_read_length": max(lengths),
        "format": detect_sequence_format(path),
        "notes": "validated_by_tandemx_streaming_reader",
    }


def write_stats(path: Path, row: dict[str, str | int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    write_stats(args.output, inspect_reads(args.reads))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
