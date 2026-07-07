"""Streaming FASTA/FASTQ sequence readers for TandemX."""

from __future__ import annotations

import csv
import gzip
import io
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence, TextIO

from tandemx.discover.rust_backend import (
    RustBackendUnavailable,
    count_sequence_paths_stats as rust_count_sequence_paths_stats,
    rust_backend_available,
)


@dataclass(frozen=True)
class SequenceRecord:
    """A normalized FASTA or FASTQ sequence record."""

    id: str
    sequence: str
    quality: str | None = None
    description: str = ""


@dataclass(frozen=True)
class SequenceStats:
    record_count: int
    total_bases: int
    max_read_length: int


class SequenceFormatError(ValueError):
    """Raised when a sequence file is empty, malformed, or internally inconsistent."""


FASTA_SUFFIXES = (".fa", ".fasta")
FASTQ_SUFFIXES = (".fq", ".fastq")
VALID_BASES = frozenset("ACGTN")


def read_sequence_records(path: Path) -> Iterator[SequenceRecord]:
    """Read FASTA/FASTQ, optionally gzip-compressed, as a streaming iterator."""
    kind = detect_sequence_format(path)
    with open_text(path) as handle:
        if kind == "fasta":
            yield from read_fasta_records(handle, path)
        else:
            yield from read_fastq_records(handle, path)


def read_sequence_records_many(paths: Path | Sequence[Path]) -> Iterator[SequenceRecord]:
    """Read one or more FASTA/FASTQ inputs as a single streaming record iterator."""
    sequence_paths = normalize_sequence_paths(paths)
    seen_ids: set[str] = set()
    for path in sequence_paths:
        for record in read_sequence_records(path):
            if record.id in seen_ids:
                raise SequenceFormatError(
                    f"Duplicate sequence id across input read files: {record.id}"
                )
            seen_ids.add(record.id)
            yield record


def normalize_sequence_paths(paths: Path | Sequence[Path]) -> tuple[Path, ...]:
    if isinstance(paths, Path):
        return (paths,)
    normalized = tuple(paths)
    if not normalized:
        raise SequenceFormatError("At least one sequence input file is required")
    return normalized


def count_sequence_records(path: Path) -> SequenceStats:
    record_count = 0
    total_bases = 0
    max_read_length = 0
    for record in read_sequence_records(path):
        record_count += 1
        read_length = len(record.sequence)
        total_bases += read_length
        max_read_length = max(max_read_length, read_length)
    return SequenceStats(
        record_count=record_count,
        total_bases=total_bases,
        max_read_length=max_read_length,
    )


def count_sequence_records_many(
    paths: Path | Sequence[Path],
    *,
    threads: int = 1,
) -> SequenceStats:
    sequence_paths = normalize_sequence_paths(paths)
    seqkit_stats = count_sequence_records_with_seqkit(sequence_paths, threads=threads)
    if seqkit_stats is not None:
        return seqkit_stats
    rust_stats = count_sequence_records_with_rust(sequence_paths, threads=threads)
    if rust_stats is not None:
        return rust_stats
    workers = max(1, min(threads, len(sequence_paths)))
    if workers == 1:
        stats = [count_sequence_records(path) for path in sequence_paths]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            stats = list(executor.map(count_sequence_records, sequence_paths))
    return SequenceStats(
        record_count=sum(item.record_count for item in stats),
        total_bases=sum(item.total_bases for item in stats),
        max_read_length=max((item.max_read_length for item in stats), default=0),
    )


def count_sequence_records_with_seqkit(
    paths: Sequence[Path],
    *,
    threads: int,
) -> SequenceStats | None:
    seqkit = shutil.which("seqkit")
    if seqkit is None:
        return None
    command = [
        seqkit,
        "stats",
        "-T",
        "-j",
        str(max(1, threads)),
        *[str(path) for path in paths],
    ]
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    try:
        return parse_seqkit_stats_table(result.stdout)
    except ValueError:
        return None


def count_sequence_records_with_rust(
    paths: Sequence[Path],
    *,
    threads: int,
) -> SequenceStats | None:
    if not rust_backend_available():
        return None
    try:
        stats = rust_count_sequence_paths_stats(paths, threads=max(1, threads))
    except RustBackendUnavailable:
        return None
    except AttributeError:
        return None
    except ValueError:
        return None
    return SequenceStats(
        record_count=stats.record_count,
        total_bases=stats.total_bases,
        max_read_length=stats.max_read_length,
    )


def parse_seqkit_stats_table(text: str) -> SequenceStats:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    required = {"num_seqs", "sum_len", "max_len"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError("seqkit stats output is missing required fields")
    rows = [row for row in reader if any((value or "").strip() for value in row.values())]
    if not rows:
        raise ValueError("seqkit stats output contains no data rows")
    return SequenceStats(
        record_count=sum(parse_seqkit_int(row["num_seqs"]) for row in rows),
        total_bases=sum(parse_seqkit_int(row["sum_len"]) for row in rows),
        max_read_length=max(parse_seqkit_int(row["max_len"]) for row in rows),
    )


def parse_seqkit_int(value: str) -> int:
    return int(value.replace(",", "").strip())


def detect_sequence_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    suffix = Path(name).suffix
    if suffix in FASTA_SUFFIXES:
        return "fasta"
    if suffix in FASTQ_SUFFIXES:
        return "fastq"
    raise SequenceFormatError(
        f"Unsupported sequence file extension for {path}. "
        "Expected .fa, .fasta, .fq, .fastq, or .gz-compressed variants."
    )


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def record_id(header: str) -> str:
    identifier = header.split()[0].split(";")[0]
    if not identifier:
        raise SequenceFormatError("Sequence record has an empty identifier")
    return identifier


def validate_sequence(sequence: str, path: Path, line_number: int) -> str:
    sequence = sequence.strip().upper()
    if not sequence:
        raise SequenceFormatError(f"Empty sequence line in {path} at line {line_number}")
    invalid = sorted(set(sequence).difference(VALID_BASES))
    if invalid:
        joined = "".join(invalid)
        raise SequenceFormatError(f"Invalid base(s) '{joined}' in {path} at line {line_number}")
    return sequence


def read_fasta_records(handle: TextIO, path: Path) -> Iterator[SequenceRecord]:
    seen_ids: set[str] = set()
    current_header: str | None = None
    parts: list[str] = []
    yielded = False
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                yield make_fasta_record(current_header, parts, seen_ids, path)
                yielded = True
            current_header = line[1:].strip()
            if not current_header:
                raise SequenceFormatError(f"Empty FASTA header in {path} at line {line_number}")
            parts = []
            continue
        if current_header is None:
            raise SequenceFormatError(f"Invalid FASTA in {path}: sequence before header at line {line_number}")
        parts.append(validate_sequence(line, path, line_number))
    if current_header is not None:
        yield make_fasta_record(current_header, parts, seen_ids, path)
        yielded = True
    if not yielded:
        raise SequenceFormatError(f"Sequence file is empty or contains no records: {path}")


def make_fasta_record(
    header: str,
    parts: list[str],
    seen_ids: set[str],
    path: Path,
) -> SequenceRecord:
    if len(parts) == 1:
        sequence = parts[0]
    else:
        sequence = "".join(parts)
    if not sequence:
        raise SequenceFormatError(f"Invalid FASTA in {path}: empty sequence for record {header}")
    identifier = record_id(header)
    if identifier in seen_ids:
        raise SequenceFormatError(f"Duplicate sequence id in {path}: {identifier}")
    seen_ids.add(identifier)
    return SequenceRecord(id=identifier, sequence=sequence, quality=None, description=header)


def read_fastq_records(handle: TextIO, path: Path) -> Iterator[SequenceRecord]:
    seen_ids: set[str] = set()
    yielded = False
    line_number = 0
    while True:
        header = handle.readline()
        if not header:
            break
        line_number += 1
        header = header.rstrip("\n\r")
        if not header:
            continue
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality:
            raise SequenceFormatError(f"Truncated FASTQ record in {path} starting at line {line_number}")
        sequence_line = line_number + 1
        line_number += 3
        header_text = header[1:].strip() if header.startswith("@") else ""
        if not header.startswith("@") or not header_text:
            raise SequenceFormatError(f"Invalid FASTQ header in {path} at line {line_number - 3}")
        if not plus.startswith("+"):
            raise SequenceFormatError(f"Invalid FASTQ separator in {path} at line {line_number - 1}")
        sequence_text = validate_sequence(sequence, path, sequence_line)
        quality_text = quality.strip()
        if len(sequence_text) != len(quality_text):
            raise SequenceFormatError(
                f"FASTQ sequence and quality lengths differ for {record_id(header_text)} in {path}"
            )
        identifier = record_id(header_text)
        if identifier in seen_ids:
            raise SequenceFormatError(f"Duplicate sequence id in {path}: {identifier}")
        seen_ids.add(identifier)
        yielded = True
        yield SequenceRecord(
            id=identifier,
            sequence=sequence_text,
            quality=quality_text,
            description=header_text,
        )
    if not yielded:
        raise SequenceFormatError(f"Sequence file is empty or contains no records: {path}")
