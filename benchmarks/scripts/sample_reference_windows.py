#!/usr/bin/env python3
"""Create deterministic nested window samples from a reference FASTA.

Every selected window is emitted as a separate FASTA record so sampling never
creates an artificial adjacency.  Selection uses complete, non-overlapping
tiles and is therefore bounded by the requested window size rather than the
reference size.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import random
import re
from typing import BinaryIO, Iterator, Protocol


SCALE_PATTERN = re.compile(r"([A-Za-z0-9][A-Za-z0-9_.-]*)=([1-9][0-9]*)")


class Digest(Protocol):
    def update(self, data: bytes) -> None: ...


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_fasta(path: Path) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def fasta_lines(
    path: Path, *, raw_digest: Digest | None = None
) -> Iterator[tuple[str, int, bytes]]:
    """Yield sequence lines with their zero-based start in each record."""
    seen: set[str] = set()
    current: str | None = None
    position = 0
    with open_fasta(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw_digest is not None:
                raw_digest.update(raw)
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                identifier = line[1:].split(maxsplit=1)[0].decode("utf-8")
                if not identifier or identifier in seen:
                    raise ValueError(
                        f"Empty or duplicate FASTA identifier at line {line_number}"
                    )
                seen.add(identifier)
                current = identifier
                position = 0
                continue
            if current is None:
                raise ValueError(f"Sequence before FASTA header at line {line_number}")
            if any(character in b" \t\r" for character in line):
                raise ValueError(f"Whitespace inside FASTA sequence at line {line_number}")
            sequence = line.upper()
            yield current, position, sequence
            position += len(sequence)
    if not seen:
        raise ValueError("FASTA contains no records")


def reference_metadata(path: Path) -> tuple[dict[str, int], str]:
    lengths: dict[str, int] = {}
    raw_digest = hashlib.sha256() if path.suffix != ".gz" else None
    for identifier, start, sequence in fasta_lines(path, raw_digest=raw_digest):
        expected = lengths.get(identifier, 0)
        if start != expected:
            raise ValueError(f"Non-contiguous FASTA coordinates for {identifier}")
        lengths[identifier] = start + len(sequence)
    input_sha256 = raw_digest.hexdigest() if raw_digest is not None else digest_file(path)
    return lengths, input_sha256


def parse_scales(values: list[str], window_size: int) -> dict[str, int]:
    if window_size <= 0:
        raise ValueError("window size must be positive")
    scales: dict[str, int] = {}
    for value in values:
        match = SCALE_PATTERN.fullmatch(value)
        if match is None:
            raise ValueError(f"Invalid scale, expected NAME=BASES: {value}")
        name, bases_text = match.groups()
        bases = int(bases_text)
        if name in scales:
            raise ValueError(f"Duplicate sample name: {name}")
        if bases % window_size:
            raise ValueError(
                f"Scale {name}={bases} is not divisible by window size {window_size}"
            )
        scales[name] = bases
    if not scales:
        raise ValueError("At least one scale is required")
    return scales


def selected_tiles(
    lengths: dict[str, int], scales: dict[str, int], window_size: int, seed: int
) -> tuple[list[dict[str, object]], dict[str, int]]:
    tiles = [
        (identifier, start, start + window_size)
        for identifier, length in lengths.items()
        for start in range(0, length - window_size + 1, window_size)
    ]
    needed = {name: bases // window_size for name, bases in scales.items()}
    maximum = max(needed.values())
    if maximum > len(tiles):
        raise ValueError(
            f"Requested {maximum} windows but reference has {len(tiles)} complete windows"
        )
    random.Random(seed).shuffle(tiles)
    selected = []
    for rank, (identifier, start, end) in enumerate(tiles[:maximum], 1):
        selected.append(
            {
                "rank": rank,
                "record_id": f"tile_{rank:06d}",
                "source_sequence": identifier,
                "source_start": start,
                "source_end": end,
                "length_bp": end - start,
                "sample_ids": [name for name, count in needed.items() if rank <= count],
            }
        )
    return selected, needed


def sample_reference(
    reference: Path,
    outdir: Path,
    scale_values: list[str],
    window_size: int,
    seed: int,
) -> dict[str, object]:
    reference = reference.resolve()
    if not reference.is_file():
        raise ValueError(f"Reference FASTA is missing: {reference}")
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    scales = parse_scales(scale_values, window_size)
    lengths, reference_sha256 = reference_metadata(reference)
    selected, needed = selected_tiles(lengths, scales, window_size, seed)
    selected_by_key = {
        (str(row["source_sequence"]), int(row["source_start"])): row
        for row in selected
    }

    outdir.mkdir(parents=True)
    output_paths = {name: outdir / f"{name}.fa" for name in scales}
    handles = {name: path.open("wb") for name, path in output_paths.items()}
    opened: set[tuple[str, str, int]] = set()
    written_bases = {name: 0 for name in scales}
    written_records = {name: 0 for name in scales}
    try:
        for identifier, line_start, sequence in fasta_lines(reference):
            offset = 0
            while offset < len(sequence):
                absolute = line_start + offset
                tile_start = (absolute // window_size) * window_size
                tile_end = tile_start + window_size
                take = min(len(sequence) - offset, tile_end - absolute)
                row = selected_by_key.get((identifier, tile_start))
                if row is not None:
                    for sample_id in row["sample_ids"]:
                        key = (sample_id, identifier, tile_start)
                        if key not in opened:
                            if absolute != tile_start:
                                raise ValueError(
                                    f"Selected tile did not begin at its source start: {key}"
                                )
                            header = f">{row['record_id']}\n"
                            handles[sample_id].write(header.encode("utf-8"))
                            opened.add(key)
                            written_records[sample_id] += 1
                        handles[sample_id].write(sequence[offset : offset + take] + b"\n")
                        written_bases[sample_id] += take
                offset += take
    finally:
        for handle in handles.values():
            handle.close()

    sample_rows: dict[str, dict[str, object]] = {}
    for name, expected_bases in scales.items():
        if written_records[name] != needed[name] or written_bases[name] != expected_bases:
            raise ValueError(
                f"Incomplete sample {name}: records={written_records[name]}, "
                f"bases={written_bases[name]}"
            )
        path = output_paths[name]
        sample_rows[name] = {
            "path": path.name,
            "requested_bases": expected_bases,
            "window_count": written_records[name],
            "total_bases": written_bases[name],
            "bytes": path.stat().st_size,
            "sha256": digest_file(path),
        }

    window_path = outdir / "windows.tsv"
    with window_path.open("w", newline="", encoding="utf-8") as handle:
        fields = (
            "rank",
            "record_id",
            "source_sequence",
            "source_start",
            "source_end",
            "length_bp",
            "sample_ids",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in sorted(
            selected,
            key=lambda item: (str(item["source_sequence"]), int(item["source_start"])),
        ):
            writer.writerow({**row, "sample_ids": ",".join(row["sample_ids"])})

    receipt: dict[str, object] = {
        "schema_version": 1,
        "complete": True,
        "method": "seeded_permutation_of_complete_nonoverlapping_reference_windows",
        "warning": "separate_window_records_preserve_source_coordinates_but_not_whole_genome_context",
        "reference": {
            "path": str(reference),
            "bytes": reference.stat().st_size,
            "sha256": reference_sha256,
            "sequence_count": len(lengths),
            "total_bases": sum(lengths.values()),
            "sequence_lengths": lengths,
        },
        "seed": seed,
        "window_size": window_size,
        "eligible_window_count": sum(length // window_size for length in lengths.values()),
        "samples": sample_rows,
        "windows": {
            "path": window_path.name,
            "rows": len(selected),
            "bytes": window_path.stat().st_size,
            "sha256": digest_file(window_path),
        },
    }
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--scale",
        action="append",
        required=True,
        help="Nested sample as NAME=BASES; may be supplied more than once",
    )
    parser.add_argument("--window-size", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=8101)
    args = parser.parse_args()
    try:
        sample_reference(
            args.reference, args.outdir, args.scale, args.window_size, args.seed
        )
    except (OSError, ValueError, UnicodeDecodeError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
