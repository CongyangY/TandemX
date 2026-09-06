"""Small, strict records shared by benchmark generation and evaluation."""

from __future__ import annotations

import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


@dataclass(frozen=True)
class ArrayRecord:
    """One array in 0-based half-open read coordinates."""

    read_id: str
    start: int
    end: int
    period: int
    sequence: str = ""
    family_id: str = ""

    def __post_init__(self) -> None:
        if not self.read_id or self.start < 0 or self.end <= self.start or self.period < 1:
            raise ValueError(f"Invalid array record: {self}")
        if self.sequence and set(self.sequence.upper()) - set("ACGTN"):
            raise ValueError("Array consensus must contain only ACGTN")


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_table(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None or isinstance(v, float) and not math.isfinite(v) else v
                          for k, v in row.items()} for row in rows)


def iter_table(path: Path, required: set[str] | None = None) -> Iterator[dict[str, str]]:
    """Validate and yield one TSV row at a time, including late malformed rows."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Duplicate TSV field: {path}")
        if required and required - set(reader.fieldnames):
            raise ValueError(f"Missing required TSV fields in {path}: {sorted(required - set(reader.fieldnames))}")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"Malformed TSV row: {path}")
            yield row


def read_table(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    return list(iter_table(path, required))
