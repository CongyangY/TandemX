"""Streaming FASTA/FASTQ sequence readers for TandemX."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


@dataclass(frozen=True)
class SequenceRecord:
    """A normalized FASTA or FASTQ sequence record."""

    id: str
    sequence: str
    quality: str | None = None
    description: str = ""


class SequenceFormatError(ValueError):
    """Raised when a sequence file is empty, malformed, or internally inconsistent."""


FASTA_SUFFIXES = (".fa", ".fasta")
FASTQ_SUFFIXES = (".fq", ".fastq")


def read_sequence_records(path: Path) -> Iterator[SequenceRecord]:
    """Read FASTA/FASTQ, optionally gzip-compressed, as a streaming iterator."""
    kind = detect_sequence_format(path)
    with open_text(path) as handle:
        if kind == "fasta":
            yield from read_fasta_records(handle, path)
        else:
            yield from read_fastq_records(handle, path)


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
    invalid = sorted(set(sequence).difference("ACGTN"))
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
