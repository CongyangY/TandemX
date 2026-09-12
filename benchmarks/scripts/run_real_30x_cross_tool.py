"""Chunked, benchmark-only real-read discovery comparison.

This runner is deliberately separate from TandemX production commands.  It
materializes a single verified FASTA into deterministic record-bound chunks and
runs TandemX, TRF and TideHunter on every one of those exact chunks. It retains
every execution receipt and produces only a preliminary exact-normalization
screen. A permissive family-matching screen remains mandatory before any
exclusive-candidate conclusion.

It does not establish biological novelty or discovery accuracy: an absent call
from a comparator is a parameter- and implementation-specific observation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import sys
from typing import Iterator

from benchmarks.challenge.adapters import build_command, iter_arrays
from benchmarks.challenge.run import run_process, source_manifest
from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.score_known_repeat_exclusion import circular_glocal_identity
from tandemx.discover.mvp import orient_monomer
from tandemx.discover.status import has_verified_empty_catalog
from tandemx.io.sequences import read_sequence_records


TOOLS = ("tandemx", "trf", "tidehunter")
COMPARATOR_TOOLS = tuple(tool for tool in TOOLS if tool != "tandemx")
DIRECT_IDENTITY = 0.90
DIRECT_SHORTER_TO_LONGER_RATIO_MIN = 0.90
MULTIPLE_IDENTITY = 0.80
MULTIPLE_RELATIVE_ERROR_MAX = 0.05
SEED_LENGTH = 4
MAX_CANDIDATES_PER_QUERY = 10_000
FROZEN_MIN_PERIOD = 30
FROZEN_MAX_PERIOD = 1000
FROZEN_MIN_REPEAT_SPAN = 100
FROZEN_TOP_PERIODS = 5
FROZEN_MIN_SUPPORT_READS = 5
FORMAL_PARTITION_COUNT = 3
FORMAL_PARTITION_TARGET_COVERAGE = 10.0
FAMILY_FIELDS = [
    "tool", "family_key", "canonical_unit", "primitive_period_bp",
    "source_record_count", "chunk_count", "source_kind",
]
EXACT_SCREEN_FIELDS = FAMILY_FIELDS + ["matched_other_tools", "screen_status", "warning"]
RECURRENCE_FIELDS = FAMILY_FIELDS + ["cross_partition_recurrence_status"]
EXACT_COMPARATOR_SUPPORT_FIELDS = [
    "tandemx_family_key", "tandemx_canonical_unit", "tandemx_length_bp",
    "tandemx_source_record_count", "tandemx_chunk_count", "cross_partition_recurrence_status",
    "comparator_tool", "comparator_family_key", "comparator_source_record_count",
    "comparator_chunk_count", "exact_comparator_support_status", "warning",
]
PERMISSIVE_FIELDS = [
    "tandemx_family_key", "tandemx_canonical_unit", "tandemx_length_bp",
    "matched_comparator", "matched_family_key", "matched_length_bp", "glocal_identity",
    "best_orientation", "relation", "shorter_to_longer_length_ratio", "nearest_integer_multiple",
    "multiple_relative_error", "candidate_count", "candidate_index_status",
    "permissive_match_status", "warning",
]


def primitive_unit(sequence: str) -> str:
    """Return the exact shortest repeated unit using a linear KMP prefix scan."""
    sequence = sequence.upper()
    if not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("Consensus must be nonempty ACGTN")
    prefix = [0] * len(sequence)
    matched = 0
    for index in range(1, len(sequence)):
        while matched and sequence[index] != sequence[matched]:
            matched = prefix[matched - 1]
        if sequence[index] == sequence[matched]:
            matched += 1
        prefix[index] = matched
    period = len(sequence) - prefix[-1]
    return sequence[:period] if len(sequence) % period == 0 else sequence


def normalized_family(sequence: str) -> tuple[str, str, int]:
    """Normalize exact period multiples, rotations and reverse complements."""
    unit = primitive_unit(sequence)
    canonical = orient_monomer(unit)
    return hashlib.sha256(canonical.encode("ascii")).hexdigest(), canonical, len(canonical)


def build_frozen_tandemx_command(reads: Path, outdir: Path, *, threads: int,
                                 executable: str = "python") -> list[str]:
    """Return the frozen YSD56 production discovery command for one partition."""
    return [
        executable, "-m", "tandemx.cli", "discover", "--reads", str(reads), "--outdir", str(outdir),
        "--discovery-method", "cascade", "--clustering-method", "sequence", "--cluster-identity", ".95",
        "--family-audit", "related", "--min-period", str(FROZEN_MIN_PERIOD),
        "--max-period", str(FROZEN_MAX_PERIOD), "--top-periods", str(FROZEN_TOP_PERIODS),
        "--min-support-reads", str(FROZEN_MIN_SUPPORT_READS), "--min-repeat-span", str(FROZEN_MIN_REPEAT_SPAN),
        "--kmer-backend", "rust", "--threads", str(threads), "--no-progress",
    ]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _circular_words(sequence: str, k: int = SEED_LENGTH) -> set[str]:
    if len(sequence) < k:
        return {sequence}
    words = set()
    for oriented in (sequence, _reverse_complement(sequence)):
        circular = oriented + oriented[: k - 1]
        words.update(circular[index:index + k] for index in range(len(oriented)))
    return words


def classify_permissive_match(identity: float, shorter_to_longer_length_ratio: float) -> tuple[str, int, float]:
    """Classify one glocal pair under the frozen cross-tool family rule."""
    if not 0 <= identity <= 1 or not 0 < shorter_to_longer_length_ratio <= 1:
        raise ValueError("identity must be in [0,1] and shorter_to_longer_length_ratio must be in (0,1]")
    longer_to_shorter = 1 / shorter_to_longer_length_ratio
    multiple = max(1, math.floor(longer_to_shorter + .5))
    error = abs(longer_to_shorter - multiple) / multiple
    if (identity + 1e-12 >= DIRECT_IDENTITY
            and shorter_to_longer_length_ratio + 1e-12 >= DIRECT_SHORTER_TO_LONGER_RATIO_MIN):
        return "direct_glocal", multiple, error
    if (multiple >= 2 and identity + 1e-12 >= MULTIPLE_IDENTITY
            and error <= MULTIPLE_RELATIVE_ERROR_MAX + 1e-12):
        return "integer_period_multiple", multiple, error
    return "no_qualifying_match", multiple, error


def _family_id_from_header(description: str, ordinal: int) -> str:
    for field in description.split(";"):
        if field.startswith("family_id=") and field[10:]:
            return field[10:]
    return f"unparsed_tandemx_family_{ordinal:09d}"


def _append_family(connection: sqlite3.Connection, *, tool: str, sequence: str,
                   source_kind: str, chunk_id: str) -> None:
    key, canonical, period = normalized_family(sequence)
    connection.execute(
        "INSERT INTO families(tool, family_key, canonical_unit, primitive_period_bp, "
        "source_record_count, chunk_count, source_kind) VALUES (?, ?, ?, ?, 1, 1, ?) "
        "ON CONFLICT(tool, family_key) DO UPDATE SET "
        "source_record_count=source_record_count+1, "
        "chunk_count=chunk_count+(SELECT CASE WHEN EXISTS(SELECT 1 FROM family_chunks "
        "WHERE tool=excluded.tool AND family_key=excluded.family_key AND chunk_id=?) THEN 0 ELSE 1 END)",
        (tool, key, canonical, period, source_kind, chunk_id),
    )
    connection.execute(
        "INSERT OR IGNORE INTO family_chunks(tool, family_key, chunk_id) VALUES (?, ?, ?)",
        (tool, key, chunk_id),
    )


def _iter_tandemx_families(path: Path) -> Iterator[str]:
    for ordinal, record in enumerate(read_sequence_records(path), 1):
        _family_id_from_header(record.description or record.id, ordinal)
        yield record.sequence


def create_chunks(input_fasta: Path, outdir: Path, database: Path, chunk_bases: int,
                  *, label: str = "chunk", required_count: int | None = None) -> dict:
    """Build common chunks and a bounded SQLite read index without loading reads."""
    if chunk_bases < 1:
        raise ValueError("--chunk-bases must be positive")
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite existing chunk directory: {outdir}")
    outdir.mkdir(parents=True)
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("CREATE TABLE reads(read_id TEXT PRIMARY KEY, length INTEGER NOT NULL, chunk_id TEXT NOT NULL) WITHOUT ROWID")
    chunks: list[dict] = []
    chunk_handle = None
    chunk_path: Path | None = None
    chunk_reads = chunk_total = read_count = total_bases = 0
    digest = hashlib.sha256()

    def close_chunk() -> None:
        nonlocal chunk_handle, chunk_path, chunk_reads, chunk_total
        if chunk_handle is None or chunk_path is None:
            return
        chunk_handle.close()
        chunks.append({"chunk_id": chunk_path.stem, "path": str(chunk_path), "sha256": digest_file(chunk_path),
                       "read_count": chunk_reads, "total_bases": chunk_total})
        chunk_handle = None
        chunk_path = None
        chunk_reads = chunk_total = 0

    try:
        for record in read_sequence_records(input_fasta):
            sequence = record.sequence
            if (chunk_handle is not None and chunk_total and chunk_total + len(sequence) > chunk_bases
                    and (required_count is None or len(chunks) < required_count - 1)):
                close_chunk()
            if chunk_handle is None:
                chunk_path = outdir / f"{label}_{len(chunks) + 1:05d}.fa"
                chunk_handle = chunk_path.open("x", encoding="ascii")
            assert chunk_path is not None
            try:
                connection.execute("INSERT INTO reads VALUES (?, ?, ?)", (record.id, len(sequence), chunk_path.stem))
            except sqlite3.IntegrityError as exc:
                raise ValueError(f"Duplicate read ID: {record.id}") from exc
            chunk_handle.write(f">{record.id}\n{sequence}\n")
            digest.update(record.id.encode("utf-8")); digest.update(b"\0")
            digest.update(sequence.encode("ascii")); digest.update(b"\n")
            chunk_reads += 1; chunk_total += len(sequence); read_count += 1; total_bases += len(sequence)
            if read_count % 10_000 == 0:
                connection.commit()
        close_chunk()
        if not chunks:
            raise ValueError("Input FASTA contains no reads")
        if required_count is not None and len(chunks) != required_count:
            raise ValueError(f"Require exactly {required_count} disjoint source partitions; observed {len(chunks)}")
        connection.commit()
    except Exception:
        if chunk_handle is not None:
            chunk_handle.close()
        raise
    finally:
        connection.close()
    result = {"input_fasta": str(input_fasta.resolve()), "input_fasta_sha256": digest_file(input_fasta),
              "semantic_sha256": digest.hexdigest(), "read_count": read_count, "total_bases": total_bases,
              "chunk_bases_target": chunk_bases, "chunks": chunks}
    if label != "chunk":
        result["partition_target_bases"] = result.pop("chunk_bases_target")
        result["partitions"] = result.pop("chunks")
    return result


def _normalize_chunk(connection: sqlite3.Connection, tool: str, native: Path, chunk_id: str,
                     min_period: int, max_period: int, min_span: int) -> dict:
    calls = 0
    if tool == "tandemx":
        monomers = native.parent / "monomers.fa"
        if not monomers.is_file():
            raise ValueError("TandemX completed without monomers.fa")
        if monomers.stat().st_size == 0:
            if not has_verified_empty_catalog(native.parent):
                raise ValueError("Empty TandemX monomer catalogue lacks a verified no_families receipt")
            return {"normalized_calls": None, "family_source": "catalogue_family"}
        for sequence in _iter_tandemx_families(monomers):
            _append_family(connection, tool=tool, sequence=sequence, source_kind="catalogue_family", chunk_id=chunk_id)
        return {"normalized_calls": None, "family_source": "catalogue_family"}
    for array in iter_arrays(tool, native, min_period, max_period, min_span):
        row = connection.execute("SELECT length, chunk_id FROM reads WHERE read_id=?", (array.read_id,)).fetchone()
        if row is None or row[1] != chunk_id or array.end > row[0]:
            raise ValueError(f"Native prediction has unknown/cross-chunk/out-of-bounds read: {array.read_id}")
        if not array.sequence:
            raise ValueError(f"{tool} native output lacks a consensus sequence")
        _append_family(connection, tool=tool, sequence=array.sequence, source_kind="array_consensus", chunk_id=chunk_id)
        calls += 1
    return {"normalized_calls": calls, "family_source": "array_consensus"}


def write_exact_match_screen(database: Path, output: Path, complete: bool,
                             *, comparator_tools: tuple[str, ...] = COMPARATOR_TOOLS,
                             unavailable_comparator_tools: tuple[str, ...] = ()) -> list[dict]:
    """Write an exact-only preliminary screen, never an exclusivity conclusion."""
    connection = sqlite3.connect(database)
    try:
        rows: list[dict] = []
        for row in connection.execute("SELECT tool, family_key, canonical_unit, primitive_period_bp, source_record_count, chunk_count, source_kind FROM families WHERE tool='tandemx' ORDER BY family_key"):
            base = dict(zip(FAMILY_FIELDS, row, strict=True))
            matched = [name for name in comparator_tools if connection.execute(
                "SELECT 1 FROM families WHERE tool=? AND family_key=?", (name, base["family_key"])).fetchone()]
            if not complete:
                status, warning = "unresolved", "one_or_more_screened_tool_chunks_failed_or_timed_out"
            elif matched:
                status, warning = "matched_by_exact_normalization_preliminary", "exact_period_multiple_rotation_and_RC_match;permissive_screen_still_pending"
            else:
                status, warning = "requires_permissive_family_match_screen", "exact_unmatched_requires_nonexact_family_comparison_before_interpretation"
            if unavailable_comparator_tools:
                warning = f"{warning};unavailable_comparator_tools={','.join(unavailable_comparator_tools)}"
            rows.append({**base, "matched_other_tools": ",".join(matched) or "none", "screen_status": status, "warning": warning})
        write_table(output, rows, EXACT_SCREEN_FIELDS)
        return rows
    finally:
        connection.close()


def write_normalized_families(database: Path, output: Path) -> None:
    """Publish the complete comparator-normalized family table for audit."""
    connection = sqlite3.connect(database)
    try:
        rows = [dict(zip(FAMILY_FIELDS, row, strict=True)) for row in connection.execute(
            "SELECT tool, family_key, canonical_unit, primitive_period_bp, source_record_count, chunk_count, source_kind "
            "FROM families ORDER BY tool, family_key"
        )]
        write_table(output, rows, FAMILY_FIELDS)
    finally:
        connection.close()


def write_cross_partition_recurrence(database: Path, output: Path, partition_count: int) -> list[dict]:
    """Write TandemX family recurrence across the fixed source partitions."""
    if partition_count < 1:
        raise ValueError("partition_count must be positive")
    connection = sqlite3.connect(database)
    try:
        rows = []
        for row in connection.execute(
            "SELECT tool, family_key, canonical_unit, primitive_period_bp, source_record_count, chunk_count, source_kind "
            "FROM families WHERE tool='tandemx' ORDER BY family_key"
        ):
            base = dict(zip(FAMILY_FIELDS, row, strict=True))
            status = ("recurred_all_formal_partitions" if int(base["chunk_count"]) == partition_count
                      else "not_recurred_all_formal_partitions")
            rows.append({**base, "cross_partition_recurrence_status": status})
        write_table(output, rows, RECURRENCE_FIELDS)
        return rows
    finally:
        connection.close()


def write_exact_comparator_support(database: Path, output: Path, partition_count: int,
                                   comparator_tool: str) -> list[dict]:
    """Join exact family support to its per-tool partition recurrence, without a uniqueness claim."""
    if comparator_tool not in COMPARATOR_TOOLS:
        raise ValueError(f"Unknown comparator tool: {comparator_tool}")
    if partition_count < 1:
        raise ValueError("partition_count must be positive")
    connection = sqlite3.connect(database)
    try:
        rows = []
        query = (
            "SELECT t.family_key, t.canonical_unit, t.primitive_period_bp, t.source_record_count, t.chunk_count, "
            "c.family_key, c.source_record_count, c.chunk_count "
            "FROM families t LEFT JOIN families c ON c.tool=? AND c.family_key=t.family_key "
            "WHERE t.tool='tandemx' ORDER BY t.family_key"
        )
        for row in connection.execute(query, (comparator_tool,)):
            (key, sequence, length, source_count, chunk_count,
             comparator_key, comparator_source_count, comparator_chunk_count) = row
            recurrence = ("recurred_all_formal_partitions" if int(chunk_count) == partition_count
                          else "not_recurred_all_formal_partitions")
            if comparator_key is None:
                support, warning = "exact_unmatched_requires_permissive_screen", "exact_unmatched_is_not_absent_comparator_support"
            else:
                support, warning = "exact_comparator_support", "exact_period_multiple_rotation_and_RC_match;permissive_screen_still_required_for_unmatched_families"
            rows.append({
                "tandemx_family_key": str(key), "tandemx_canonical_unit": str(sequence),
                "tandemx_length_bp": str(length), "tandemx_source_record_count": str(source_count),
                "tandemx_chunk_count": str(chunk_count), "cross_partition_recurrence_status": recurrence,
                "comparator_tool": comparator_tool, "comparator_family_key": str(comparator_key or "NA"),
                "comparator_source_record_count": str(comparator_source_count or "NA"),
                "comparator_chunk_count": str(comparator_chunk_count or "NA"),
                "exact_comparator_support_status": support, "warning": warning,
            })
        write_table(output, rows, EXACT_COMPARATOR_SUPPORT_FIELDS)
        return rows
    finally:
        connection.close()


def _build_permissive_candidate_index(connection: sqlite3.Connection) -> set[int]:
    """Index comparator representatives by length and rotation/strand-aware seeds."""
    connection.execute(
        "CREATE TABLE permissive_family_seeds(tool TEXT NOT NULL, seed_length INTEGER NOT NULL, "
        "seed TEXT NOT NULL, family_key TEXT NOT NULL, PRIMARY KEY(tool, seed_length, seed, family_key)) WITHOUT ROWID"
    )
    lengths: set[int] = set()
    for tool, key, sequence, length in connection.execute(
        "SELECT tool, family_key, canonical_unit, primitive_period_bp FROM families "
        "WHERE tool IN ('trf', 'tidehunter')"
    ):
        lengths.add(int(length))
        size = min(SEED_LENGTH, len(sequence))
        connection.executemany(
            "INSERT OR IGNORE INTO permissive_family_seeds VALUES (?, ?, ?, ?)",
            ((tool, size, word, key) for word in _circular_words(sequence, size)),
        )
    connection.execute("CREATE INDEX permissive_seed_lookup ON permissive_family_seeds(tool, seed_length, seed, family_key)")
    connection.commit()
    return lengths


def _candidate_lengths(query_length: int, comparator_lengths: set[int]) -> tuple[int, ...]:
    return tuple(
        length for length in sorted(comparator_lengths)
        if classify_permissive_match(1.0, min(query_length, length) / max(query_length, length))[0] != "no_qualifying_match"
    )


def _candidate_rows(connection: sqlite3.Connection, tool: str, query_sequence: str,
                    comparator_lengths: set[int], maximum_candidates: int) -> tuple[list[tuple[str, str, int]], bool]:
    """Return bounded seed-and-length candidates; saturation is unresolved, never negative evidence."""
    allowed_lengths = _candidate_lengths(len(query_sequence), comparator_lengths)
    if not allowed_lengths:
        return [], False
    seed_size = min(SEED_LENGTH, len(query_sequence))
    seeds = sorted(_circular_words(query_sequence, seed_size))
    candidates: dict[str, tuple[str, str, int]] = {}
    for start in range(0, len(allowed_lengths), 500):
        lengths = allowed_lengths[start:start + 500]
        placeholders = ",".join("?" for _ in seeds)
        length_placeholders = ",".join("?" for _ in lengths)
        remaining = maximum_candidates - len(candidates) + 1
        query = (
            "SELECT f.family_key, f.canonical_unit, f.primitive_period_bp "
            "FROM permissive_family_seeds s JOIN families f "
            "ON f.tool=s.tool AND f.family_key=s.family_key "
            f"WHERE s.tool=? AND s.seed_length=? AND s.seed IN ({placeholders}) "
            f"AND f.primitive_period_bp IN ({length_placeholders}) "
            "GROUP BY f.family_key, f.canonical_unit, f.primitive_period_bp LIMIT ?"
        )
        for key, sequence, length in connection.execute(query, (tool, seed_size, *seeds, *lengths, remaining)):
            candidates[str(key)] = (str(key), str(sequence), int(length))
        if len(candidates) >= maximum_candidates:
            return list(candidates.values())[:maximum_candidates], True
    return list(candidates.values()), False


def _seed_length_mismatch_possible(query_sequence: str, comparator_lengths: set[int]) -> bool:
    """A short primitive unit cannot use the 4-mer index of a longer candidate."""
    query_seed_length = min(SEED_LENGTH, len(query_sequence))
    return any(
        min(SEED_LENGTH, length) != query_seed_length
        for length in _candidate_lengths(len(query_sequence), comparator_lengths)
    )


def screen_permissive_matches(database: Path, output: Path, complete: bool) -> tuple[list[dict], dict]:
    """Run the frozen permissive screen with bounded candidate generation.

    The 4-mer seed index is only a candidate generator.  Exact circular glocal
    alignment supplies every accepted result.  A candidate-limit saturation is
    explicitly unresolved so it cannot become a negative comparator result.
    """
    connection = sqlite3.connect(database)
    try:
        queries = list(connection.execute(
            "SELECT family_key, canonical_unit, primitive_period_bp FROM families WHERE tool='tandemx' ORDER BY family_key"
        ))
        rows: list[dict] = []
        saturated_queries = seed_length_mismatch_queries = 0
        comparator_lengths: set[int] = set()
        if complete:
            comparator_lengths = _build_permissive_candidate_index(connection)
        for key, sequence, length in queries:
            base = {
                "tandemx_family_key": str(key), "tandemx_canonical_unit": str(sequence),
                "tandemx_length_bp": str(length), "matched_comparator": "NA", "matched_family_key": "NA",
                "matched_length_bp": "NA", "glocal_identity": "NA", "best_orientation": "NA",
                "relation": "NA", "shorter_to_longer_length_ratio": "NA", "nearest_integer_multiple": "NA",
                "multiple_relative_error": "NA", "candidate_count": "0", "candidate_index_status": "not_run",
            }
            if not complete:
                rows.append({**base, "permissive_match_status": "technical_unresolved",
                             "warning": "one_or_more_tool_chunks_failed_or_timed_out"})
                continue
            matches: list[tuple[int, float, float, str, str, int, str, int, float]] = []
            candidate_count = 0
            saturated = False
            seed_mismatch = _seed_length_mismatch_possible(str(sequence), comparator_lengths)
            for tool in ("trf", "tidehunter"):
                remaining = MAX_CANDIDATES_PER_QUERY - candidate_count
                if remaining <= 0:
                    saturated = True
                    break
                candidates, tool_saturated = _candidate_rows(connection, tool, str(sequence), comparator_lengths, remaining)
                candidate_count += len(candidates)
                saturated = saturated or tool_saturated
                for target_key, target_sequence, target_length in candidates:
                    identity, _distance, orientation = circular_glocal_identity(str(sequence), target_sequence)
                    shorter_to_longer = min(int(length), target_length) / max(int(length), target_length)
                    relation, multiple, multiple_error = classify_permissive_match(identity, shorter_to_longer)
                    if relation != "no_qualifying_match":
                        matches.append((2 if relation == "direct_glocal" else 1, identity, -multiple_error,
                                        tool, target_key, target_length, orientation, multiple, multiple_error))
            if matches:
                priority, identity, _negative_error, tool, target_key, target_length, orientation, multiple, multiple_error = max(
                    matches, key=lambda value: (value[0], value[1], value[2], value[3], value[4])
                )
                relation = "direct_glocal" if priority == 2 else "integer_period_multiple"
                shorter_to_longer = min(int(length), target_length) / max(int(length), target_length)
                status, warning = "matched_by_permissive_screen", "comparator_recurrence_under_frozen_glocal_rule"
                index_status = "saturated_but_qualifying_match_found" if saturated else "complete"
            elif saturated:
                saturated_queries += 1
                status, warning, index_status = "permissive_screen_unresolved_candidate_limit", "candidate_index_limit_prevents_negative_match_claim", "saturated"
                tool = target_key = orientation = "NA"; target_length = multiple = "NA"; identity = multiple_error = shorter_to_longer = None; relation = "NA"
            elif seed_mismatch:
                seed_length_mismatch_queries += 1
                status, warning, index_status = "permissive_screen_unresolved_seed_length_mismatch", "short_primitive_seed_length_mismatch_prevents_negative_match_claim", "seed_length_mismatch"
                tool = target_key = orientation = "NA"; target_length = multiple = "NA"; identity = multiple_error = shorter_to_longer = None; relation = "NA"
            else:
                status, warning, index_status = "eligible_for_recurrence_and_biological_triage", "no_qualifying_comparator_representative_under_frozen_screen", "complete"
                tool = target_key = orientation = "NA"; target_length = multiple = "NA"; identity = multiple_error = shorter_to_longer = None; relation = "no_qualifying_match"
            rows.append({**base, "matched_comparator": tool, "matched_family_key": target_key,
                         "matched_length_bp": str(target_length), "glocal_identity": "NA" if identity is None else f"{identity:.6f}",
                         "best_orientation": orientation, "relation": relation,
                         "shorter_to_longer_length_ratio": "NA" if shorter_to_longer is None else f"{shorter_to_longer:.6f}",
                         "nearest_integer_multiple": str(multiple),
                         "multiple_relative_error": "NA" if multiple_error is None else f"{multiple_error:.6f}",
                         "candidate_count": str(candidate_count), "candidate_index_status": index_status,
                         "permissive_match_status": status, "warning": warning})
        write_table(output, rows, PERMISSIVE_FIELDS)
        receipt = {
            "schema_version": 1, "screen_complete": bool(complete and saturated_queries == 0 and seed_length_mismatch_queries == 0),
            "query_family_count": len(queries), "candidate_index_saturated_queries": saturated_queries,
            "seed_length_mismatch_queries": seed_length_mismatch_queries,
            "candidate_seed_length_bp": SEED_LENGTH, "maximum_candidates_per_query": MAX_CANDIDATES_PER_QUERY,
            "direct_glocal_identity_minimum": DIRECT_IDENTITY,
            "direct_shorter_to_longer_length_ratio_minimum": DIRECT_SHORTER_TO_LONGER_RATIO_MIN,
            "period_multiple_identity_minimum": MULTIPLE_IDENTITY,
            "period_multiple_relative_error_maximum": MULTIPLE_RELATIVE_ERROR_MAX,
            "output_sha256": digest_file(output),
            "interpretation_boundary": "cross-tool recurrence screen; eligible rows still require recurrence and biological triage and are neither novel nor tool-exclusive",
        }
        return rows, receipt
    finally:
        connection.close()


def run(input_fasta: Path, outdir: Path, trf: Path, tidehunter: Path, timeout: float,
        chunk_bases: int | None, threads: int, *, assembly_span_bp: int,
        min_period: int = FROZEN_MIN_PERIOD, max_period: int = FROZEN_MAX_PERIOD,
        min_span: int = FROZEN_MIN_REPEAT_SPAN, plan_only: bool = False) -> dict:
    if timeout <= 0 or not math.isfinite(timeout) or not 1 <= threads <= 64 or assembly_span_bp < 1:
        raise ValueError("timeout must be finite and positive; threads must be in [1,64]")
    if (min_period, max_period, min_span) != (FROZEN_MIN_PERIOD, FROZEN_MAX_PERIOD, FROZEN_MIN_REPEAT_SPAN):
        raise ValueError("The formal 30x comparison requires the frozen 30--1000 bp and 100-bp span scope")
    if outdir.exists():
        raise FileExistsError(f"Choose a new output directory: {outdir}")
    if not input_fasta.is_file():
        raise ValueError(f"Missing input FASTA: {input_fasta}")
    executables = {"tandemx": Path(sys.executable), "trf": trf, "tidehunter": tidehunter}
    if not plan_only and any(not path.is_file() or not os.access(path, os.X_OK) for path in executables.values()):
        raise ValueError("Missing executable")
    outdir.mkdir(parents=True)
    database = outdir / "normalization.sqlite"
    partition_target_bases = round(assembly_span_bp * FORMAL_PARTITION_TARGET_COVERAGE)
    if chunk_bases is not None and chunk_bases != partition_target_bases:
        raise ValueError("--chunk-bases is deprecated for formal runs and must equal 10x --assembly-span-bp")
    manifest = create_chunks(input_fasta, outdir / "partitions", database, partition_target_bases,
                             label="partition", required_count=FORMAL_PARTITION_COUNT)
    partition_depths = []
    for partition in manifest["partitions"]:
        partition["partition_id"] = partition.pop("chunk_id")
        observed_depth = partition["total_bases"] / assembly_span_bp
        partition["observed_depth"] = observed_depth
        partition_depths.append(observed_depth)
    if any(abs(depth - FORMAL_PARTITION_TARGET_COVERAGE) > 1.0 for depth in partition_depths):
        raise ValueError("Each formal source partition must be within 1x of the 10x target depth")
    root = Path(__file__).resolve().parents[2]
    manifest.update(source_manifest(root), schema_version=1, tools={name: str(path.resolve()) for name, path in executables.items()},
                    tool_sha256={name: digest_file(path) for name, path in executables.items()},
                    runner_sha256=digest_file(Path(__file__)),
                    scope={"min_period": min_period, "max_period": max_period, "min_span": min_span},
                    formal_partitions={
                        "partition_count": FORMAL_PARTITION_COUNT,
                        "target_coverage_per_partition": FORMAL_PARTITION_TARGET_COVERAGE,
                        "assembly_span_bp": assembly_span_bp,
                        "target_bases_per_partition": partition_target_bases,
                        "observed_partition_depths": partition_depths,
                        "source_order": "preserved_and_disjoint",
                    },
                    timeout_per_chunk_seconds=timeout, threads=threads,
                    resource_scope="direct_child_wait4; subprocess descendants may be excluded; sum CPU and max RSS are descriptive",
                    exact_match_screen_rule="exact primitive period multiple plus cyclic rotation and reverse-complement normalization;preliminary_only",
                    permissive_family_match_rule={
                        "direct_glocal_identity_minimum": DIRECT_IDENTITY,
                        "direct_shorter_to_longer_length_ratio_minimum": DIRECT_SHORTER_TO_LONGER_RATIO_MIN,
                        "integer_period_multiple_identity_minimum": MULTIPLE_IDENTITY,
                        "integer_period_multiple_relative_error_maximum": MULTIPLE_RELATIVE_ERROR_MAX,
                        "candidate_seed_length_bp": SEED_LENGTH,
                        "maximum_candidates_per_query": MAX_CANDIDATES_PER_QUERY,
                        "glocal_alignment_helper_sha256": digest_file(
                            Path(circular_glocal_identity.__code__.co_filename)
                        ),
                    },
                    accuracy="not_assessed_without_independent_real-read truth",
                    storage_warning=("none" if str(outdir.resolve()).startswith("/Volumes/T7/Codex/TandemX/")
                                     else "large_real_read_outputs_should_be_placed_under_/Volumes/T7/Codex/TandemX"))
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if plan_only:
        return {"status": "planned", "partition_count": len(manifest["partitions"])}
    connection = sqlite3.connect(database)
    connection.execute("CREATE TABLE families(tool TEXT NOT NULL, family_key TEXT NOT NULL, canonical_unit TEXT NOT NULL, primitive_period_bp INTEGER NOT NULL, source_record_count INTEGER NOT NULL, chunk_count INTEGER NOT NULL, source_kind TEXT NOT NULL, PRIMARY KEY(tool, family_key)) WITHOUT ROWID")
    connection.execute("CREATE TABLE family_chunks(tool TEXT NOT NULL, family_key TEXT NOT NULL, chunk_id TEXT NOT NULL, PRIMARY KEY(tool, family_key, chunk_id)) WITHOUT ROWID")
    rows: list[dict] = []
    try:
        for partition in manifest["partitions"]:
            chunk_path = Path(partition["path"])
            for tool in TOOLS:
                folder = outdir / "runs" / tool / partition["partition_id"]
                folder.mkdir(parents=True)
                command, native = build_command(tool, str(executables[tool]), chunk_path, folder, min_period, max_period, min_span, threads)
                if tool == "tandemx":
                    command = build_frozen_tandemx_command(chunk_path, native.parent, threads=threads,
                                                           executable=sys.executable)
                measured = run_process(command, native if tool == "trf" else folder / "stdout.log", folder / "stderr.log", timeout,
                                       working_directory=root)
                state = "completed" if measured["exit_code"] == 0 and not measured["timed_out"] else "failed_or_timed_out"
                detail = {"tool": tool, "partition_id": partition["partition_id"], "command": command, "native_output": str(native), "status": state, **measured}
                if state == "completed":
                    try:
                        detail.update(_normalize_chunk(connection, tool, native, partition["partition_id"], min_period, max_period, min_span), normalization="ok")
                        connection.commit()
                    except Exception as exc:
                        detail.update(status="normalization_failed", normalization="failed", failure=str(exc))
                        connection.rollback()
                (folder / "execution.json").write_text(json.dumps(detail, indent=2) + "\n")
                rows.append(detail)
    finally:
        connection.close()
    summary_fields = ["tool", "partition_id", "status", "exit_code", "timed_out", "runtime_seconds", "cpu_user_seconds", "cpu_system_seconds", "peak_rss_mib", "normalization", "normalized_calls", "family_source", "failure"]
    write_table(outdir / "execution_summary.tsv", [{field: row.get(field) for field in summary_fields} for row in rows], summary_fields)
    complete = len(rows) == len(manifest["partitions"]) * len(TOOLS) and all(row["status"] == "completed" and row.get("normalization") == "ok" for row in rows)
    write_normalized_families(database, outdir / "all_tool_normalized_families.tsv")
    exact_screen = write_exact_match_screen(database, outdir / "tandemx_exact_match_screen.tsv", complete)
    permissive_rows, permissive_receipt = screen_permissive_matches(
        database, outdir / "tandemx_permissive_family_match.tsv", complete
    )
    (outdir / "permissive_family_match_receipt.json").write_text(json.dumps(permissive_receipt, indent=2) + "\n")
    permissive_complete = bool(permissive_receipt["screen_complete"])
    summary = {
        "status": ("completed_permissive_screen" if permissive_complete else
                   "completed_with_unresolved_permissive_screen" if complete else
                   "completed_with_technical_unresolved_comparator"),
        "partition_count": len(manifest["partitions"]),
        "execution_count": len(rows),
        "exact_match_screen_status": "complete_preliminary_only" if complete else "unresolved",
        "permissive_family_match_status": "complete" if permissive_complete else "technical_or_candidate_index_unresolved",
        "tandemx_family_rows": len(exact_screen),
        "eligible_for_recurrence_and_biological_triage_rows": sum(
            row["permissive_match_status"] == "eligible_for_recurrence_and_biological_triage"
            for row in permissive_rows
        ),
    }
    (outdir / "completion.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--trf", type=Path, required=True)
    parser.add_argument("--tidehunter", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=7200)
    parser.add_argument("--assembly-span-bp", type=int, required=True,
                        help="Verified assembly span; formal input is split into three disjoint ~10x partitions.")
    parser.add_argument("--chunk-bases", type=int,
                        help="Deprecated formal guard; when supplied it must equal 10x --assembly-span-bp.")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-period", type=int, default=30)
    parser.add_argument("--max-period", type=int, default=1000)
    parser.add_argument("--min-span", type=int, default=100)
    parser.add_argument("--plan-only", action="store_true", help="Materialize and hash common chunks, but do not execute tools")
    args = parser.parse_args()
    result = run(args.input_fasta, args.outdir, args.trf, args.tidehunter, args.timeout, args.chunk_bases, args.threads,
                 assembly_span_bp=args.assembly_span_bp, min_period=args.min_period, max_period=args.max_period,
                 min_span=args.min_span, plan_only=args.plan_only)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
