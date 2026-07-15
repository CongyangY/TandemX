"""Streaming FASTA/FASTQ sequence readers for TandemX."""

from __future__ import annotations

import csv
import gzip
import io
import sqlite3
import shutil
import subprocess
import tempfile
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


@dataclass(frozen=True)
class FastaChunk:
    """One bounded sequence chunk with its 0-based offset in a FASTA record."""

    id: str
    description: str
    start: int
    sequence: str


class SequenceFormatError(ValueError):
    """Raised when a sequence file is empty, malformed, or internally inconsistent."""


class DuplicateIdTracker:
    """Exact duplicate detector that spills identifiers to SQLite at a fixed limit."""

    def __init__(self, memory_limit: int = 100_000) -> None:
        if memory_limit < 1:
            raise ValueError("memory_limit must be positive")
        self.memory_limit = memory_limit
        self._memory: set[str] = set()
        self._database: sqlite3.Connection | None = None
        self._database_path: Path | None = None

    @property
    def spilled(self) -> bool:
        return self._database is not None

    def add(self, identifier: str) -> bool:
        """Add an identifier and return False when it was already present."""
        if self._database is None and len(self._memory) < self.memory_limit:
            if identifier in self._memory:
                return False
            self._memory.add(identifier)
            return True
        if self._database is None:
            self._spill()
        assert self._database is not None
        try:
            self._database.execute("INSERT INTO sequence_ids(identifier) VALUES (?)", (identifier,))
        except sqlite3.IntegrityError:
            return False
        return True

    def _spill(self) -> None:
        temporary = tempfile.NamedTemporaryFile(prefix="tandemx-ids-", suffix=".sqlite3", delete=False)
        temporary.close()
        self._database_path = Path(temporary.name)
        self._database = sqlite3.connect(self._database_path)
        self._database.execute("PRAGMA journal_mode=OFF")
        self._database.execute("PRAGMA synchronous=OFF")
        self._database.execute("CREATE TABLE sequence_ids(identifier TEXT PRIMARY KEY) WITHOUT ROWID")
        self._database.executemany(
            "INSERT INTO sequence_ids(identifier) VALUES (?)",
            ((identifier,) for identifier in self._memory),
        )
        self._memory.clear()

    def close(self) -> None:
        if self._database is not None:
            self._database.close()
            self._database = None
        if self._database_path is not None:
            self._database_path.unlink(missing_ok=True)
            self._database_path = None
        self._memory.clear()

    def __enter__(self) -> DuplicateIdTracker:
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()


FASTA_SUFFIXES = (".fa", ".fasta")
FASTQ_SUFFIXES = (".fq", ".fastq")
VALID_BASES = frozenset("ACGTN")


def read_sequence_records(path: Path) -> Iterator[SequenceRecord]:
    """Read FASTA/FASTQ, optionally gzip-compressed, as a streaming iterator."""
    with DuplicateIdTracker() as seen_ids:
        yield from _read_sequence_records(path, seen_ids)


def _read_sequence_records(path: Path, seen_ids: DuplicateIdTracker) -> Iterator[SequenceRecord]:
    kind = detect_sequence_format(path)
    with open_text(path) as handle:
        if kind == "fasta":
            yield from read_fasta_records(handle, path, seen_ids=seen_ids)
        else:
            yield from read_fastq_records(handle, path, seen_ids=seen_ids)


def read_sequence_records_many(
    paths: Path | Sequence[Path],
    *,
    check_duplicate_ids_across_files: bool = True,
) -> Iterator[SequenceRecord]:
    """Read one or more FASTA/FASTQ inputs as a single streaming record iterator."""
    sequence_paths = normalize_sequence_paths(paths)
    if not check_duplicate_ids_across_files:
        for path in sequence_paths:
            yield from read_sequence_records(path)
        return
    with DuplicateIdTracker() as seen_ids:
        for path in sequence_paths:
            try:
                yield from _read_sequence_records(path, seen_ids)
            except SequenceFormatError as exc:
                if "Duplicate sequence id in" in str(exc):
                    identifier = str(exc).rsplit(": ", 1)[-1]
                    raise SequenceFormatError(
                        f"Duplicate sequence id across input read files: {identifier}"
                    ) from exc
                raise


def normalize_sequence_paths(paths: Path | Sequence[Path]) -> tuple[Path, ...]:
    if isinstance(paths, Path):
        return (paths,)
    normalized = tuple(paths)
    if not normalized:
        raise SequenceFormatError("At least one sequence input file is required")
    return normalized


def read_fasta_chunks(path: Path, *, chunk_bases: int = 1_000_000) -> Iterator[FastaChunk]:
    """Stream bounded chunks without materializing an entire assembly contig."""
    if chunk_bases <= 0:
        raise ValueError("chunk_bases must be positive")
    if detect_sequence_format(path) != "fasta":
        raise SequenceFormatError(f"Assembly input must be FASTA: {path}")
    with DuplicateIdTracker() as seen_ids, open_text(path) as handle:
        current_id: str | None = None
        current_description = ""
        parts: list[str] = []
        buffered_bases = 0
        offset = 0
        record_has_sequence = False
        yielded = False
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    if not record_has_sequence:
                        raise SequenceFormatError(
                            f"Invalid FASTA in {path}: empty sequence for record {current_description}"
                        )
                    if parts:
                        sequence = "".join(parts)
                        yield FastaChunk(current_id, current_description, offset, sequence)
                        yielded = True
                current_description = line[1:].strip()
                if not current_description:
                    raise SequenceFormatError(f"Empty FASTA header in {path} at line {line_number}")
                current_id = record_id(current_description)
                if not seen_ids.add(current_id):
                    raise SequenceFormatError(f"Duplicate sequence id in {path}: {current_id}")
                parts = []
                buffered_bases = 0
                offset = 0
                record_has_sequence = False
                continue
            if current_id is None:
                raise SequenceFormatError(
                    f"Invalid FASTA in {path}: sequence before header at line {line_number}"
                )
            validated = validate_sequence(line, path, line_number)
            parts.append(validated)
            buffered_bases += len(validated)
            record_has_sequence = True
            if buffered_bases >= chunk_bases:
                combined = "".join(parts)
                while len(combined) >= chunk_bases:
                    sequence = combined[:chunk_bases]
                    yield FastaChunk(current_id, current_description, offset, sequence)
                    yielded = True
                    offset += len(sequence)
                    combined = combined[chunk_bases:]
                parts = [combined] if combined else []
                buffered_bases = len(combined)
        if current_id is not None:
            if not record_has_sequence:
                raise SequenceFormatError(
                    f"Invalid FASTA in {path}: empty sequence for record {current_description}"
                )
            if parts:
                sequence = "".join(parts)
                yield FastaChunk(current_id, current_description, offset, sequence)
                yielded = True
        if not yielded:
            raise SequenceFormatError(f"Sequence file is empty or contains no records: {path}")


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


def read_fasta_records(
    handle: TextIO,
    path: Path,
    *,
    seen_ids: DuplicateIdTracker | None = None,
) -> Iterator[SequenceRecord]:
    if seen_ids is None:
        with DuplicateIdTracker() as tracker:
            yield from _read_fasta_records(handle, path, tracker)
        return
    yield from _read_fasta_records(handle, path, seen_ids)


def _read_fasta_records(
    handle: TextIO,
    path: Path,
    seen_ids: DuplicateIdTracker,
) -> Iterator[SequenceRecord]:
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
    seen_ids: DuplicateIdTracker,
    path: Path,
) -> SequenceRecord:
    if len(parts) == 1:
        sequence = parts[0]
    else:
        sequence = "".join(parts)
    if not sequence:
        raise SequenceFormatError(f"Invalid FASTA in {path}: empty sequence for record {header}")
    identifier = record_id(header)
    if not seen_ids.add(identifier):
        raise SequenceFormatError(f"Duplicate sequence id in {path}: {identifier}")
    return SequenceRecord(id=identifier, sequence=sequence, quality=None, description=header)


def read_fastq_records(
    handle: TextIO,
    path: Path,
    *,
    seen_ids: DuplicateIdTracker | None = None,
) -> Iterator[SequenceRecord]:
    if seen_ids is None:
        with DuplicateIdTracker() as tracker:
            yield from _read_fastq_records(handle, path, tracker)
        return
    yield from _read_fastq_records(handle, path, seen_ids)


def _read_fastq_records(
    handle: TextIO,
    path: Path,
    seen_ids: DuplicateIdTracker,
) -> Iterator[SequenceRecord]:
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
        if not seen_ids.add(identifier):
            raise SequenceFormatError(f"Duplicate sequence id in {path}: {identifier}")
        yielded = True
        yield SequenceRecord(
            id=identifier,
            sequence=sequence_text,
            quality=quality_text,
            description=header_text,
        )
    if not yielded:
        raise SequenceFormatError(f"Sequence file is empty or contains no records: {path}")
