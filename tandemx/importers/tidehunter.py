"""Strict TideHunter ``-f 2`` import with source-read validation."""
from __future__ import annotations

import csv
from dataclasses import dataclass, replace
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import sqlite3
import tempfile
from typing import Iterator, Sequence
from urllib.parse import quote

from tandemx.discover.clustering import cluster_monomers, write_membership
from tandemx.discover.family_audit import write_family_audit
from tandemx.discover.mvp import (
    CandidateRepeat,
    is_low_complexity,
    wrap_sequence,
    write_candidate_reads,
    write_families,
    write_monomers,
)
from tandemx.discover.status import file_sha256, write_discovery_summary
from tandemx.io.sequences import read_sequence_records_many


IMPORT_AUDIT_FIELDS = [
    "candidate_id",
    "read_id",
    "native_read_id",
    "tidehunter_repeat_id",
    "native_copy_number",
    "native_read_length",
    "native_start",
    "native_end",
    "native_consensus_length",
    "native_average_match_percent",
    "native_full_length_copy_number",
    "native_subunit_starts",
    "normalized_start",
    "normalized_end",
    "status",
    "warning",
]


@dataclass(frozen=True)
class TideHunterRecord:
    read_id: str
    repeat_id: str
    copy_number: float
    read_length: int
    start: int
    end: int
    period: int
    average_match_percent: float
    full_length_copy_number: int
    subunit_starts: tuple[int, ...]
    consensus: str


def iter_tidehunter_f2(path: Path) -> Iterator[TideHunterRecord]:
    """Yield strict native records; coordinates remain 1-based inclusive here."""
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 11:
                raise ValueError(
                    f"Expected TideHunter -f 2 output with 11 fields: {path}:{line_number}"
                )
            try:
                copy_number = float(fields[2])
                read_length = int(fields[3])
                start = int(fields[4])
                end = int(fields[5])
                period = int(fields[6])
                average_match = float(fields[7])
                full_length = int(fields[8])
                subunit_starts = tuple(int(value) for value in fields[9].split(","))
            except ValueError as error:
                raise ValueError(
                    f"Non-numeric TideHunter -f 2 field: {path}:{line_number}"
                ) from error
            consensus = fields[10].upper()
            if (
                not fields[0]
                or not fields[1]
                or not math.isfinite(copy_number)
                or copy_number <= 0
                or read_length < 1
                or not 1 <= start <= end <= read_length
                or period < 1
                or not math.isfinite(average_match)
                or not 0 <= average_match <= 100
                or full_length < 0
                or not subunit_starts
                or any(value < start or value > end for value in subunit_starts)
            ):
                raise ValueError(f"Invalid TideHunter -f 2 values: {path}:{line_number}")
            if any(a >= b for a, b in zip(subunit_starts, subunit_starts[1:])):
                raise ValueError(
                    f"TideHunter subunit starts are not strictly increasing: {path}:{line_number}"
                )
            if not consensus or set(consensus) - set("ACGTN"):
                raise ValueError(
                    f"TideHunter consensus must contain nonempty ACGTN: {path}:{line_number}"
                )
            if len(consensus) != period:
                raise ValueError(
                    f"TideHunter consensus length differs from field 7: {path}:{line_number}"
                )
            yield TideHunterRecord(
                fields[0],
                fields[1],
                copy_number,
                read_length,
                start,
                end,
                period,
                average_match,
                full_length,
                subunit_starts,
                consensus,
            )


def _read_index(
    reads: Sequence[Path], database: Path
) -> tuple[sqlite3.Connection, int, int, str]:
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute(
        "CREATE TABLE reads (read_id TEXT PRIMARY KEY, native_id TEXT UNIQUE NOT NULL, "
        "length INTEGER NOT NULL) WITHOUT ROWID"
    )
    digest = hashlib.sha256()
    read_count = total_bases = 0
    try:
        for record in read_sequence_records_many(reads):
            identifier = record.id.encode("utf-8")
            sequence = record.sequence.encode("ascii")
            native_id = record.description.split()[0]
            digest.update(len(identifier).to_bytes(8, "big"))
            digest.update(identifier)
            digest.update(len(sequence).to_bytes(8, "big"))
            digest.update(sequence)
            try:
                connection.execute(
                    "INSERT INTO reads(read_id, native_id, length) VALUES (?, ?, ?)",
                    (record.id, native_id, len(record.sequence)),
                )
            except sqlite3.IntegrityError as error:  # pragma: no cover; reader rejects first
                raise ValueError(f"Duplicate read identifier: {record.id}") from error
            read_count += 1
            total_bases += len(record.sequence)
            if read_count % 10_000 == 0:
                connection.commit()
        connection.commit()
        if read_count < 1 or total_bases < 1:
            raise ValueError("TideHunter import requires nonempty source reads")
        return connection, read_count, total_bases, digest.hexdigest()
    except Exception:
        connection.close()
        raise


def _candidate(
    record: TideHunterRecord, ordinal: int, normalized_read_id: str
) -> CandidateRepeat:
    low_complexity = is_low_complexity(record.consensus)
    warnings = [
        "external_detector=tidehunter",
        "native_format=f2",
        "source_coordinates=1_based_inclusive",
        "uncalibrated_confidence",
    ]
    if low_complexity:
        warnings.append("low_complexity_candidate")
    if "N" in record.consensus:
        warnings.append("ambiguous_consensus_bases")
    return CandidateRepeat(
        read_id=normalized_read_id,
        candidate_id=f"TXC{ordinal:09d}",
        sequence=record.consensus,
        read_start=record.start - 1,
        read_end=record.end,
        strand=".",
        period_bp=record.period,
        repeat_span_bp=record.end - record.start + 1,
        unit_count=record.copy_number,
        score=record.average_match_percent / 100,
        low_complexity_flag=low_complexity,
        confidence=(
            "high"
            if record.average_match_percent >= 90 and not low_complexity and "N" not in record.consensus
            else "medium"
        ),
        warning=";".join(warnings),
    )


def import_tidehunter_catalog(
    native: Path,
    reads: Sequence[Path],
    outdir: Path,
    *,
    min_support_reads: int = 2,
    cluster_identity: float = 0.95,
    backend: str = "python",
    family_audit: str = "related",
    kmer_size: int = 11,
    logger: logging.Logger | None = None,
) -> dict[str, object]:
    """Build a TandemX catalogue while retaining external-detector provenance."""
    if min_support_reads < 1:
        raise ValueError("--min-support-reads must be positive")
    if not 0 < cluster_identity <= 1:
        raise ValueError("--cluster-identity must be in (0,1]")
    if backend not in {"python", "rust"}:
        raise ValueError("--backend must resolve to python or rust")
    if family_audit not in {"full", "related"}:
        raise ValueError("--family-audit must be full or related")
    if kmer_size < 1:
        raise ValueError("--kmer-size must be positive")
    protected = (
        "candidate_reads.tsv",
        "candidate_monomers.fa",
        "monomer_membership.tsv",
        "families.tsv",
        "monomers.fa",
        "family_similarity.tsv",
        "family_audit_summary.json",
        "discovery_summary.json",
        "tidehunter_import.tsv",
        "import_summary.json",
    )
    existing = [name for name in protected if (outdir / name).exists()]
    if existing:
        raise ValueError(f"Refusing to overwrite existing import outputs: {','.join(existing)}")
    logger = logger or logging.getLogger("tandemx.import.tidehunter")
    native_sha256 = file_sha256(native)
    descriptor, database_name = tempfile.mkstemp(
        prefix=".tidehunter-read-index-", suffix=".sqlite3", dir=outdir
    )
    os.close(descriptor)
    database = Path(database_name)
    connection: sqlite3.Connection | None = None
    try:
        connection, read_count, total_bases, semantic_sha256 = _read_index(reads, database)
        candidates: list[CandidateRepeat] = []
        audit_rows: list[dict[str, object]] = []
        for ordinal, record in enumerate(iter_tidehunter_f2(native), 1):
            source = connection.execute(
                "SELECT read_id, length FROM reads WHERE native_id=?", (record.read_id,)
            ).fetchone()
            if source is None:
                raise ValueError(
                    f"TideHunter output contains read absent from source inputs: {record.read_id}"
                )
            if int(source[1]) != record.read_length:
                raise ValueError(
                    f"TideHunter read length differs from source input: {record.read_id}"
                )
            candidate = _candidate(record, ordinal, str(source[0]))
            candidates.append(candidate)
            audit_rows.append(
                {
                    "candidate_id": candidate.candidate_id,
                    "read_id": candidate.read_id,
                    "native_read_id": record.read_id,
                    "tidehunter_repeat_id": record.repeat_id,
                    "native_copy_number": f"{record.copy_number:g}",
                    "native_read_length": record.read_length,
                    "native_start": record.start,
                    "native_end": record.end,
                    "native_consensus_length": record.period,
                    "native_average_match_percent": f"{record.average_match_percent:g}",
                    "native_full_length_copy_number": record.full_length_copy_number,
                    "native_subunit_starts": ",".join(map(str, record.subunit_starts)),
                    "normalized_start": candidate.read_start,
                    "normalized_end": candidate.read_end,
                    "status": "imported",
                    "warning": "external_detector=tidehunter;coordinates_normalized_to_0_based_half_open",
                }
            )
        if file_sha256(native) != native_sha256:
            raise ValueError("TideHunter native output changed during import")
    finally:
        if connection is not None:
            connection.close()
        database.unlink(missing_ok=True)

    write_candidate_reads(outdir / "candidate_reads.tsv", candidates)
    with (outdir / "candidate_monomers.fa").open("w", encoding="utf-8") as handle:
        for candidate in candidates:
            handle.write(
                f">candidate_id={candidate.candidate_id};"
                f"read_id={quote(candidate.read_id, safe='')};"
                f"length_bp={len(candidate.sequence)}\n{wrap_sequence(candidate.sequence)}\n"
            )
    with (outdir / "tidehunter_import.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=IMPORT_AUDIT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    families, membership = cluster_monomers(
        candidates, min_support_reads, cluster_identity, backend
    )
    families = [
        replace(
            family,
            warning=";".join(
                [part for part in (family.warning, "external_detector=tidehunter") if part]
            ),
        )
        for family in families
    ]
    write_membership(outdir / "monomer_membership.tsv", membership)
    families, _redundant = write_family_audit(
        outdir / "family_similarity.tsv",
        families,
        k=kmer_size,
        backend=backend,
        mode=family_audit,
        logger=logger,
    )
    write_monomers(outdir / "monomers.fa", families)
    write_families(outdir / "families.tsv", families)
    write_discovery_summary(
        outdir, read_count, total_bases, len(candidates), len(families)
    )
    outputs = [
        "candidate_reads.tsv",
        "candidate_monomers.fa",
        "monomer_membership.tsv",
        "families.tsv",
        "monomers.fa",
        "family_similarity.tsv",
        "family_audit_summary.json",
        "discovery_summary.json",
        "tidehunter_import.tsv",
    ]
    summary: dict[str, object] = {
        "schema_version": 1,
        "complete": True,
        "source_detector": "TideHunter",
        "native_format": "-f 2",
        "native_input": str(native.resolve()),
        "native_input_sha256": native_sha256,
        "reads": [
            {"path": str(path.resolve()), "bytes": path.stat().st_size} for path in reads
        ],
        "reads_semantic_sha256": semantic_sha256,
        "read_count": read_count,
        "total_bases": total_bases,
        "candidate_count": len(candidates),
        "family_count": len(families),
        "minimum_support_reads": min_support_reads,
        "cluster_identity": cluster_identity,
        "backend": backend,
        "family_audit": family_audit,
        "kmer_size": kmer_size,
        "output_sha256": {name: file_sha256(outdir / name) for name in outputs},
        "warning": (
            "external_detector_calls_not_tandemx_native_discovery;"
            "reads_semantic_sha256_hashes_ordered_ids_and_sequences_not_quality"
        ),
    }
    (outdir / "import_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary
