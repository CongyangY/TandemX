"""Bounded four-line FASTQ parsing with an optional exact raw-file digest."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import gzip
import hashlib
import io
from pathlib import Path
from typing import BinaryIO, Iterator


class DigestReader(io.RawIOBase):
    """Hash raw bytes as consumed, including compressed headers and trailers."""

    def __init__(self, raw: BinaryIO):
        self.raw = raw
        self.digest = hashlib.sha256()

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        count = self.raw.readinto(buffer)
        if count:
            self.digest.update(memoryview(buffer)[:count])
        return count


@contextmanager
def hashed_fastq(path: Path):
    """The digest is final only after records have been consumed through EOF."""
    with path.open('rb') as raw:
        stream = DigestReader(raw)
        with io.BufferedReader(stream) as buffered:
            if path.suffix == '.gz':
                with gzip.GzipFile(fileobj=buffered, mode='rb') as handle:
                    yield handle, stream.digest
            else:
                yield buffered, stream.digest


@dataclass(frozen=True)
class FastqRecord:
    header: bytes
    sequence: bytes
    plus: bytes
    quality: bytes
    identifier: bytes

    def as_bytes(self) -> bytes:
        return b'\n'.join((self.header, self.sequence, self.plus, self.quality)) + b'\n'


def records(handle: BinaryIO, line_limit: int = 5_000_002) -> Iterator[FastqRecord]:
    """Validate individual records; duplicate identifiers require a separate index."""
    count = 0
    while header := handle.readline(line_limit):
        count += 1
        lines = [header] + [handle.readline(line_limit) for _ in range(3)]
        if any(not line or len(line) >= line_limit for line in lines):
            raise ValueError(f'Incomplete or oversized FASTQ record {count}')
        header, seq, plus, quality = [line.rstrip(b'\r\n') for line in lines]
        if (not header.startswith(b'@') or not plus.startswith(b'+') or not seq
                or len(seq) != len(quality) or set(seq.upper()) - set(b'ACGTN')
                or min(quality) < 33 or max(quality) > 126):
            raise ValueError(f'Invalid FASTQ record {count}')
        identifiers = header[1:].split()
        if not identifiers:
            raise ValueError('Empty FASTQ identifier')
        identifier = identifiers[0]
        if plus[1:] and (not plus[1:].split() or plus[1:].split()[0] != identifier):
            raise ValueError('FASTQ plus identifier does not match header')
        yield FastqRecord(header, seq, plus, quality, identifier)
