#!/usr/bin/env python3
"""Audit the frozen YSD56 TandemX--TideHunter permissive family rule.

This is a benchmark-only postprocessor.  It never invokes discovery or
abundance estimation.  Given a completed real-30x run directory, it verifies
the execution receipt, makes a small 4-mer-posting preflight measurement, and
then scores every length-eligible TandemX--TideHunter representative pair with
the frozen edlib circular-glocal rule.  The 4-mer measurement is diagnostic
only: final matching retrieves exact length-eligible pairs from a local SQLite
index in bounded batches and does not use a seed index,
so posting saturation cannot create a negative TideHunter result.

If a remount has removed the source run, the script deliberately writes only a
``blocked_missing_post_remount_source_outputs`` receipt with no support counts.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile
from typing import Any, Iterator

from benchmarks.scripts.run_real_30x_cross_tool import (
    MAX_CANDIDATES_PER_QUERY,
    SEED_LENGTH,
    _circular_words,
    classify_permissive_match,
)
from benchmarks.scripts.score_known_repeat_exclusion import circular_glocal_identity


SOURCE_FILES = ("normalization.sqlite", "execution_summary.tsv", "run_manifest.json")
FORMAL_PARTITION_COUNT = 3
DEFAULT_MAX_CANDIDATE_PAIRS = 5_000_000
MATCH_FIELDS = (
    "tandemx_family_key", "tandemx_length_bp", "cross_partition_recurrence_status", "exact_tidehunter_support",
    "tidehunter_family_key", "tidehunter_length_bp", "glocal_identity",
    "tidehunter_chunk_count", "tidehunter_recurrence_status", "best_orientation", "relation", "shorter_to_longer_length_ratio",
    "nearest_integer_multiple", "multiple_relative_error", "candidate_pair_count",
    "postprocess_status", "interpretation_boundary",
)
LENGTH_BREAKDOWN_FIELDS = (
    "tandemx_length_bp", "recurrent_exact_unmatched_query_count", "direct_length_eligible_pair_count",
    "integer_multiple_length_eligible_pair_count", "length_eligible_pair_union_count",
)
ABUNDANCE_BREAKDOWN_FIELDS = (
    "source_record_count_quantile", "query_count", "direct_length_eligible_pair_count",
    "integer_multiple_length_eligible_pair_count", "length_eligible_pair_union_count",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _valid_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(letter in "0123456789abcdef" for letter in value.lower())


def run_manifest_contract(manifest: Path) -> tuple[bool, str, tuple[str, ...]]:
    """Validate the portable frozen-input facts required for this replay."""
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        return False, f"run_manifest_unreadable:{type(error).__name__}", ()
    if not isinstance(data, dict):
        return False, "run_manifest_not_json_object", ()
    if not _valid_sha256(data.get("input_fasta_sha256")) or not _valid_sha256(data.get("semantic_sha256")):
        return False, "run_manifest_missing_valid_input_hashes", ()
    partitions = data.get("partitions")
    if not isinstance(partitions, list) or len(partitions) != FORMAL_PARTITION_COUNT:
        return False, "run_manifest_missing_three_formal_partitions", ()
    identifiers: list[str] = []
    for partition in partitions:
        if not isinstance(partition, dict):
            return False, "run_manifest_invalid_partition_entry", ()
        identifier = partition.get("partition_id")
        if not isinstance(identifier, str) or not identifier or not _valid_sha256(partition.get("sha256")):
            return False, "run_manifest_invalid_partition_identity_or_hash", ()
        if not isinstance(partition.get("read_count"), int) or partition["read_count"] <= 0:
            return False, "run_manifest_invalid_partition_read_count", ()
        if not isinstance(partition.get("total_bases"), int) or partition["total_bases"] <= 0:
            return False, "run_manifest_invalid_partition_base_count", ()
        identifiers.append(identifier)
    if len(set(identifiers)) != FORMAL_PARTITION_COUNT:
        return False, "run_manifest_duplicate_partition_ids", ()
    return True, "validated_frozen_input_and_three_partitions", tuple(sorted(identifiers))


def tidehunter_execution_complete(summary: Path, expected_partitions: tuple[str, ...]) -> tuple[bool, str]:
    try:
        with summary.open(encoding="utf-8", newline="") as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("tool") == "tidehunter"]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        return False, f"execution_summary_unreadable:{type(error).__name__}"
    if len(rows) != len(expected_partitions) or {row.get("partition_id") for row in rows} != set(expected_partitions):
        return False, "tidehunter_execution_rows_do_not_match_frozen_partitions"
    invalid = [row for row in rows if row.get("status") != "completed" or row.get("normalization") != "ok"]
    if invalid:
        return False, "tidehunter_execution_or_normalization_incomplete"
    return True, "all_frozen_tidehunter_partitions_completed_and_normalized"


def family_histogram(database: Path, tool: str, *, recurrent_exact_unmatched_only: bool = False) -> dict[int, int]:
    uri = f"file:{database.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    try:
        if recurrent_exact_unmatched_only:
            statement = (
                "SELECT t.primitive_period_bp, COUNT(*) FROM families t "
                "WHERE t.tool='tandemx' AND t.chunk_count=? AND NOT EXISTS "
                "(SELECT 1 FROM families c WHERE c.tool='tidehunter' AND c.family_key=t.family_key) "
                "GROUP BY t.primitive_period_bp"
            )
            values: tuple[object, ...] = (FORMAL_PARTITION_COUNT,)
        else:
            statement = "SELECT primitive_period_bp, COUNT(*) FROM families WHERE tool=? GROUP BY primitive_period_bp"
            values = (tool,)
        return {int(length): int(count) for length, count in connection.execute(statement, values)}
    finally:
        connection.close()


def length_eligible(first_length: int, second_length: int) -> bool:
    ratio = min(first_length, second_length) / max(first_length, second_length)
    return classify_permissive_match(1.0, ratio)[0] != "no_qualifying_match"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def sample_families(database: Path, tool: str, limit: int, *, recurrent_exact_unmatched_only: bool = False) -> list[dict[str, object]]:
    uri = f"file:{database.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    try:
        if recurrent_exact_unmatched_only:
            statement = (
                "SELECT t.family_key, t.canonical_unit, t.primitive_period_bp FROM families t "
                "WHERE t.tool='tandemx' AND t.chunk_count=? AND NOT EXISTS "
                "(SELECT 1 FROM families c WHERE c.tool='tidehunter' AND c.family_key=t.family_key) "
                "ORDER BY t.family_key LIMIT ?"
            )
            values: tuple[object, ...] = (FORMAL_PARTITION_COUNT, limit)
        else:
            statement = "SELECT family_key, canonical_unit, primitive_period_bp FROM families WHERE tool=? ORDER BY family_key LIMIT ?"
            values = (tool, limit)
        return [{"family_key": str(key), "sequence": str(sequence), "length": int(length)}
                for key, sequence, length in connection.execute(statement, values)]
    finally:
        connection.close()


def seed_preflight(
    queries: list[dict[str, object]], comparators: list[dict[str, object]], *, query_limit: int, comparator_limit: int
) -> dict[str, object]:
    """Measure a bounded 4-mer posting sample without using it for matching."""
    sampled_queries = queries[:query_limit]
    sampled_comparators = comparators[:comparator_limit]
    postings: dict[str, set[str]] = {}
    for family in sampled_comparators:
        for seed in _circular_words(str(family["sequence"]), SEED_LENGTH):
            postings.setdefault(seed, set()).add(str(family["family_key"]))
    union_sizes = []
    for family in sampled_queries:
        candidates: set[str] = set()
        for seed in _circular_words(str(family["sequence"]), SEED_LENGTH):
            candidates.update(postings.get(seed, set()))
        union_sizes.append(len(candidates))
    maximum_posting = max((len(value) for value in postings.values()), default=0)
    maximum_union = max(union_sizes, default=0)
    return {
        "method": "sampled_rotation_and_strand_aware_4mer_postings_diagnostic_only",
        "query_sample_count": len(sampled_queries),
        "comparator_sample_count": len(sampled_comparators),
        "seed_length_bp": SEED_LENGTH,
        "maximum_single_seed_posting": maximum_posting,
        "maximum_union_candidate_count": maximum_union,
        "historical_candidate_limit": MAX_CANDIDATES_PER_QUERY,
        "would_saturate_historical_4mer_index": maximum_union >= MAX_CANDIDATES_PER_QUERY,
        "interpretation": "preflight only; final pair selection is exhaustive among frozen length-eligible representatives",
    }


def eligible_lengths(query_length: int, comparator_lengths: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(length for length in comparator_lengths if length_eligible(query_length, length))


def candidate_pair_estimate(query_histogram: dict[int, int], comparator_histogram: dict[int, int]) -> int:
    """Count every length-eligible pair before any sequence alignment."""
    comparator_lengths = tuple(sorted(comparator_histogram))
    return sum(
        query_count * sum(comparator_histogram[length] for length in eligible_lengths(query_length, comparator_lengths))
        for query_length, query_count in query_histogram.items()
    )


def length_pair_breakdown(query_length: int, comparator_histogram: dict[int, int]) -> tuple[int, int, int]:
    """Exact pre-alignment candidate count by the two frozen length branches."""
    direct = multiple = 0
    for comparator_length, count in comparator_histogram.items():
        ratio = min(query_length, comparator_length) / max(query_length, comparator_length)
        relation, _nearest, _error = classify_permissive_match(1.0, ratio)
        if relation == "direct_glocal":
            direct += count
        elif relation == "integer_period_multiple":
            multiple += count
    return direct, multiple, direct + multiple


def write_pair_breakdowns(database: Path, outdir: Path, query_histogram: dict[int, int],
                          comparator_histogram: dict[int, int]) -> dict[str, str]:
    """Write exact, alignment-free length and source-count distribution receipts."""
    by_length = outdir / "candidate_pair_breakdown_by_query_length.tsv"
    length_rows: list[dict[str, object]] = []
    cache: dict[int, tuple[int, int, int]] = {}
    for length, query_count in sorted(query_histogram.items()):
        cache[length] = length_pair_breakdown(length, comparator_histogram)
        direct, multiple, union = cache[length]
        length_rows.append({"tandemx_length_bp": length, "recurrent_exact_unmatched_query_count": query_count,
                            "direct_length_eligible_pair_count": direct * query_count,
                            "integer_multiple_length_eligible_pair_count": multiple * query_count,
                            "length_eligible_pair_union_count": union * query_count})
    with by_length.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LENGTH_BREAKDOWN_FIELDS, delimiter="\t")
        writer.writeheader(); writer.writerows(length_rows)
    uri = f"file:{database.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    try:
        profiles = [(str(key), int(length), int(source_count)) for key, length, source_count in connection.execute(
            "SELECT t.family_key, t.primitive_period_bp, t.source_record_count FROM families t "
            "WHERE t.tool='tandemx' AND t.chunk_count=? AND NOT EXISTS "
            "(SELECT 1 FROM families c WHERE c.tool='tidehunter' AND c.family_key=t.family_key) ORDER BY t.source_record_count, t.family_key",
            (FORMAL_PARTITION_COUNT,),
        )]
    finally:
        connection.close()
    bins: dict[str, list[int]] = {f"Q{index}": [0, 0, 0, 0] for index in range(1, 5)}
    for index, (_key, length, _source_count) in enumerate(profiles):
        label = f"Q{min(4, index * 4 // max(1, len(profiles)) + 1)}"
        direct, multiple, union = cache[length]
        bins[label][0] += 1; bins[label][1] += direct; bins[label][2] += multiple; bins[label][3] += union
    by_abundance = outdir / "candidate_pair_breakdown_by_source_record_count_quantile.tsv"
    with by_abundance.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ABUNDANCE_BREAKDOWN_FIELDS, delimiter="\t")
        writer.writeheader()
        for label, (query_count, direct, multiple, union) in bins.items():
            writer.writerow({"source_record_count_quantile": label, "query_count": query_count,
                             "direct_length_eligible_pair_count": direct,
                             "integer_multiple_length_eligible_pair_count": multiple,
                             "length_eligible_pair_union_count": union})
    return {by_length.name: sha256(by_length), by_abundance.name: sha256(by_abundance)}


def qgram_necessary_condition(first: str, second: str, identity_threshold: float, *, q: int = 4) -> bool:
    """A lossless necessary q-gram condition for the frozen edlib HW score.

    If edit distance is at most ``floor((1-threshold)*len(shorter))``, each edit
    disrupts at most ``q`` q-grams of the complete shorter query.  Therefore the
    q-gram multiset intersection with the doubled longer target must retain the
    stated lower bound in at least one query orientation.  This helper is only
    used for a diagnostic preflight; final support calls still align every pair.
    """
    shorter, longer = (first, second) if len(first) <= len(second) else (second, first)
    if len(shorter) < q:
        return True
    maximum_edits = int((1.0 - identity_threshold) * len(shorter) + 1e-12)
    required_shared = max(0, len(shorter) - q + 1 - q * maximum_edits)
    target_counts = Counter((longer + longer)[index:index + q] for index in range(len(longer) * 2 - q + 1))
    for oriented in (shorter, reverse_complement(shorter)):
        query_counts = Counter(oriented[index:index + q] for index in range(len(oriented) - q + 1))
        shared = sum(min(count, target_counts.get(word, 0)) for word, count in query_counts.items())
        if shared >= required_shared:
            return True
    return False


def qgram_sample_preflight(queries: list[dict[str, object]], comparators: list[dict[str, object]]) -> dict[str, object]:
    """Measure the q-gram necessary-condition retention without changing calls."""
    total = retained = 0
    for query in queries:
        for comparator in comparators:
            ratio = min(int(query["length"]), int(comparator["length"])) / max(int(query["length"]), int(comparator["length"]))
            relation, _multiple, _error = classify_permissive_match(1.0, ratio)
            if relation == "no_qualifying_match":
                continue
            total += 1
            threshold = 0.90 if relation == "direct_glocal" else 0.80
            retained += qgram_necessary_condition(str(query["sequence"]), str(comparator["sequence"]), threshold)
    return {"q": 4, "length_eligible_sample_pairs": total, "qgram_necessary_condition_retained_pairs": retained,
            "qgram_necessary_condition_removed_pairs": total - retained,
            "interpretation": "diagnostic necessary-condition preflight only; not used to call unmatched comparator support"}


def candidate_store_required_bytes(comparator_sequence_bytes: int) -> int:
    """Conservative local temporary-index budget; source SQLite remains read-only."""
    return max(256 << 20, comparator_sequence_bytes * 3)


def build_candidate_store(source: Path, directory: Path) -> Path:
    """Copy only TideHunter representatives into a local, indexed temporary store."""
    target = directory / "tidehunter_length_index.sqlite"
    source_uri = f"file:{source.resolve()}?mode=ro"
    reader = sqlite3.connect(source_uri, uri=True)
    writer = sqlite3.connect(target)
    try:
        writer.execute("CREATE TABLE families(length INTEGER NOT NULL, family_key TEXT PRIMARY KEY, sequence TEXT NOT NULL, chunk_count INTEGER NOT NULL) WITHOUT ROWID")
        writer.execute("BEGIN")
        cursor = reader.execute("SELECT family_key, canonical_unit, primitive_period_bp, chunk_count FROM families WHERE tool='tidehunter' ORDER BY family_key")
        for ordinal, (key, sequence, length, chunk_count) in enumerate(cursor, 1):
            writer.execute("INSERT INTO families VALUES (?, ?, ?, ?)", (int(length), str(key), str(sequence), int(chunk_count)))
            if ordinal % 10_000 == 0:
                writer.commit()
                writer.execute("BEGIN")
        writer.commit()
        writer.execute("CREATE INDEX families_by_length ON families(length, family_key)")
        writer.commit()
    finally:
        writer.close()
        reader.close()
    return target


def recurrent_exact_unmatched_count(database: Path) -> tuple[int, int]:
    uri = f"file:{database.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    try:
        recurrent = int(connection.execute(
            "SELECT COUNT(*) FROM families WHERE tool='tandemx' AND chunk_count=?", (FORMAL_PARTITION_COUNT,)
        ).fetchone()[0])
        exact = int(connection.execute(
            "SELECT COUNT(*) FROM families t WHERE t.tool='tandemx' AND t.chunk_count=? AND EXISTS "
            "(SELECT 1 FROM families c WHERE c.tool='tidehunter' AND c.family_key=t.family_key)",
            (FORMAL_PARTITION_COUNT,),
        ).fetchone()[0])
        return recurrent, exact
    finally:
        connection.close()


def comparator_sequence_bytes(database: Path) -> int:
    uri = f"file:{database.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    try:
        return int(connection.execute(
            "SELECT COALESCE(SUM(LENGTH(canonical_unit)), 0) FROM families WHERE tool='tidehunter'"
        ).fetchone()[0])
    finally:
        connection.close()


def iter_recurrent_exact_unmatched_queries(database: Path) -> Iterator[tuple[str, str, int]]:
    uri = f"file:{database.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    try:
        for key, sequence, length in connection.execute(
            "SELECT t.family_key, t.canonical_unit, t.primitive_period_bp FROM families t "
            "WHERE t.tool='tandemx' AND t.chunk_count=? AND NOT EXISTS "
            "(SELECT 1 FROM families c WHERE c.tool='tidehunter' AND c.family_key=t.family_key) ORDER BY t.family_key",
            (FORMAL_PARTITION_COUNT,),
        ):
            yield str(key), str(sequence), int(length)
    finally:
        connection.close()


def best_permissive_match(query_sequence: str, query_length: int, candidates: sqlite3.Connection,
                          comparator_lengths: tuple[int, ...], *, block_size: int) -> tuple[dict[str, object] | None, int]:
    """Score all length-eligible stored representatives; no seed index can remove a pair."""
    lengths = eligible_lengths(query_length, comparator_lengths)
    if not lengths:
        return None, 0
    placeholders = ",".join("?" for _ in lengths)
    candidates_seen = 0
    best: tuple[int, float, float, str, str, int, int, str, int, float] | None = None
    cursor = candidates.execute(
        f"SELECT family_key, sequence, length, chunk_count FROM families WHERE length IN ({placeholders}) ORDER BY family_key", lengths
    )
    while batch := cursor.fetchmany(block_size):
        for target_key, target_sequence, target_length, target_chunk_count in batch:
            candidates_seen += 1
            identity, _distance, orientation = circular_glocal_identity(query_sequence, str(target_sequence))
            ratio = min(query_length, int(target_length)) / max(query_length, int(target_length))
            relation, multiple, error = classify_permissive_match(identity, ratio)
            if relation == "no_qualifying_match":
                continue
            proposal = (2 if relation == "direct_glocal" else 1, identity, -error, str(target_key),
                        str(target_sequence), int(target_length), int(target_chunk_count), orientation, multiple, error)
            if best is None or proposal[:4] > best[:4]:
                best = proposal
    if best is None:
        return None, candidates_seen
    priority, identity, _negative_error, target_key, target_sequence, target_length, target_chunk_count, orientation, multiple, error = best
    return {"target": {"family_key": target_key, "sequence": target_sequence, "length": target_length, "chunk_count": target_chunk_count},
            "identity": identity, "orientation": orientation,
            "relation": "direct_glocal" if priority == 2 else "integer_period_multiple",
            "multiple": multiple, "error": error,
            "ratio": min(query_length, target_length) / max(query_length, target_length)}, candidates_seen


def postprocess(run_dir: Path, outdir: Path, *, block_size: int = 500, preflight_query_limit: int = 50,
                preflight_comparator_limit: int = 2000,
                max_candidate_pairs: int = DEFAULT_MAX_CANDIDATE_PAIRS) -> dict[str, object]:
    """Write a reproducible receipt and, when inputs exist, exact/permissive support tables."""
    outdir.mkdir(parents=True, exist_ok=False)
    missing = [str(run_dir / name) for name in SOURCE_FILES if not (run_dir / name).is_file()]
    receipt: dict[str, Any] = {
        "schema_version": 1,
        "run_dir": str(run_dir.resolve()),
        "required_source_files": list(SOURCE_FILES),
        "matcher": "complete_shorter_vs_doubled_longer_edlib_HW_both_orientations",
        "frozen_rule": {
            "direct_glocal_identity_minimum": 0.90,
            "direct_shorter_to_longer_length_ratio_minimum": 0.90,
            "integer_period_multiple_identity_minimum": 0.80,
            "integer_period_multiple_relative_error_maximum": 0.05,
        },
        "final_candidate_method": "length_eligible_pairwise_edlib_in_bounded_comparator_blocks; no_4mer_index_for_final_matching",
        "script_sha256": sha256(Path(__file__)),
    }
    if missing:
        receipt.update(
            status="blocked_missing_post_remount_source_outputs", missing_source_files=missing,
            exact_tidehunter_supported_family_count="NA", permissive_tidehunter_supported_family_count="NA",
            unresolved_family_count="NA",
            interpretation_boundary="No family support number was inferred: the post-remount source run outputs are absent.",
        )
        (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    if (run_dir / "normalization.sqlite").stat().st_size == 0:
        receipt.update(
            status="blocked_zero_byte_normalization_database", exact_tidehunter_supported_family_count="NA",
            permissive_tidehunter_supported_family_count="NA", unresolved_family_count="NA",
            interpretation_boundary="A zero-byte normalized family database cannot establish comparator support or absence.",
        )
        (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    manifest_valid, manifest_state, expected_partitions = run_manifest_contract(run_dir / "run_manifest.json")
    receipt["run_manifest_state"] = manifest_state
    if not manifest_valid:
        receipt.update(
            status="blocked_invalid_frozen_run_manifest", exact_tidehunter_supported_family_count="NA",
            permissive_tidehunter_supported_family_count="NA", unresolved_family_count="NA",
            interpretation_boundary="The frozen input and three-partition contract is not verifiable; no comparator-support number was inferred.",
        )
        (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    complete, execution_state = tidehunter_execution_complete(run_dir / "execution_summary.tsv", expected_partitions)
    receipt["source_sha256"] = {name: sha256(run_dir / name) for name in SOURCE_FILES}
    receipt["tidehunter_execution_state"] = execution_state
    if not complete:
        receipt.update(
            status="blocked_incomplete_tidehunter_execution", exact_tidehunter_supported_family_count="NA",
            permissive_tidehunter_supported_family_count="NA", unresolved_family_count="NA",
            interpretation_boundary="Incomplete TideHunter execution is technical unresolved, never absent comparator support.",
        )
        (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    database = run_dir / "normalization.sqlite"
    recurrent_count, exact_count = recurrent_exact_unmatched_count(database)
    query_histogram = family_histogram(database, "tandemx", recurrent_exact_unmatched_only=True)
    comparator_histogram = family_histogram(database, "tidehunter")
    query_count = sum(query_histogram.values())
    comparator_count = sum(comparator_histogram.values())
    candidate_pairs = candidate_pair_estimate(query_histogram, comparator_histogram)
    queries = sample_families(database, "tandemx", preflight_query_limit, recurrent_exact_unmatched_only=True)
    comparators = sample_families(database, "tidehunter", preflight_comparator_limit)
    preflight = seed_preflight(queries, comparators, query_limit=preflight_query_limit,
                               comparator_limit=preflight_comparator_limit)
    qgram_preflight = qgram_sample_preflight(queries, comparators)
    preflight.update({
        "screen_population": "cross_partition_recurrent_TandemX_families_exact_unmatched_against_TideHunter",
        "formal_partition_count": FORMAL_PARTITION_COUNT,
        "recurrent_tandemx_family_count": recurrent_count,
        "exact_tidehunter_supported_recurrent_family_count": exact_count,
        "recurrent_exact_unmatched_tandemx_family_count": query_count,
        "tidehunter_family_count": comparator_count,
        "length_eligible_pair_alignment_estimate": candidate_pairs,
        "maximum_permitted_pair_alignments": max_candidate_pairs,
        "tidehunter_recurrence_requirement": "not_required_by_frozen_support_screen; chunk_count is reported for every matched TideHunter representative",
    })
    (outdir / "candidate_index_preflight.json").write_text(json.dumps(preflight, indent=2) + "\n", encoding="utf-8")
    (outdir / "qgram_necessary_condition_preflight.json").write_text(json.dumps(qgram_preflight, indent=2) + "\n", encoding="utf-8")
    breakdown_sha256 = write_pair_breakdowns(database, outdir, query_histogram, comparator_histogram)
    table = outdir / "ysd56_tandemx_tidehunter_permissive_matches.tsv"
    if candidate_pairs > max_candidate_pairs:
        with table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MATCH_FIELDS, delimiter="\t")
            writer.writeheader()
            for query_key, _sequence, query_length in iter_recurrent_exact_unmatched_queries(database):
                writer.writerow({
                    "tandemx_family_key": query_key, "tandemx_length_bp": str(query_length),
                    "cross_partition_recurrence_status": "recurred_all_formal_partitions", "exact_tidehunter_support": "false",
                    "tidehunter_family_key": "NA", "tidehunter_length_bp": "NA", "glocal_identity": "NA",
                    "tidehunter_chunk_count": "NA", "tidehunter_recurrence_status": "NA",
                    "best_orientation": "NA", "relation": "NA", "shorter_to_longer_length_ratio": "NA",
                    "nearest_integer_multiple": "NA", "multiple_relative_error": "NA", "candidate_pair_count": "NA",
                    "postprocess_status": "technical_unresolved_candidate_pair_budget",
                    "interpretation_boundary": "The predeclared exhaustive alignment budget was exceeded. This row is unresolved, not absent TideHunter support or TandemX-only evidence.",
                })
        receipt.update(
            status="blocked_candidate_pair_estimate_exceeds_gate",
            recurrent_tandemx_family_count=recurrent_count,
            exact_tidehunter_supported_family_count=exact_count,
            recurrent_exact_unmatched_tandemx_family_count=query_count,
            permissive_only_tidehunter_supported_family_count="NA",
            permissive_tidehunter_supported_family_count="NA",
            unmatched_no_tidehunter_support_not_tandemx_only_count="NA",
            unresolved_family_count=query_count,
            preflight_file="candidate_index_preflight.json",
            qgram_preflight_file="qgram_necessary_condition_preflight.json",
            output_sha256={table.name: sha256(table), "candidate_index_preflight.json": sha256(outdir / "candidate_index_preflight.json"),
                           "qgram_necessary_condition_preflight.json": sha256(outdir / "qgram_necessary_condition_preflight.json"), **breakdown_sha256},
            interpretation_boundary="Exact canonical-key support is reported for recurrent TandemX families. The nonexact recurrent population exceeded the exhaustive alignment budget and remains unresolved; no negative comparator conclusion was produced.",
        )
        (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    required_bytes = candidate_store_required_bytes(comparator_sequence_bytes(database))
    available_bytes = shutil.disk_usage(outdir).free
    receipt["candidate_store_required_bytes"] = required_bytes
    receipt["candidate_store_available_bytes"] = available_bytes
    if available_bytes < required_bytes:
        receipt.update(
            status="blocked_candidate_store_storage_budget", recurrent_tandemx_family_count=recurrent_count,
            exact_tidehunter_supported_family_count=exact_count, recurrent_exact_unmatched_tandemx_family_count=query_count,
            permissive_only_tidehunter_supported_family_count="NA", permissive_tidehunter_supported_family_count="NA",
            unmatched_no_tidehunter_support_not_tandemx_only_count="NA", unresolved_family_count=query_count,
            preflight_file="candidate_index_preflight.json", qgram_preflight_file="qgram_necessary_condition_preflight.json",
            interpretation_boundary="The local temporary length index lacks storage budget. Nonexact recurrent families remain unresolved, not negative comparator evidence.",
        )
        (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    permissive_only_count = unmatched_count = 0
    comparator_lengths = tuple(sorted(comparator_histogram))
    with tempfile.TemporaryDirectory(prefix="tidehunter_length_index_", dir=outdir) as temporary:
        store_path = build_candidate_store(database, Path(temporary))
        store = sqlite3.connect(store_path)
        try:
            with table.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=MATCH_FIELDS, delimiter="\t")
                writer.writeheader()
                for query_key, sequence, query_length in iter_recurrent_exact_unmatched_queries(database):
                    base = {
                        "tandemx_family_key": query_key, "tandemx_length_bp": str(query_length),
                        "cross_partition_recurrence_status": "recurred_all_formal_partitions", "exact_tidehunter_support": "false",
                        "interpretation_boundary": "TideHunter-only comparator postprocess. An unmatched row is not TandemX-only because other required comparators and biological triage remain separate.",
                    }
                    match, candidate_count = best_permissive_match(sequence, query_length, store, comparator_lengths,
                                                                    block_size=block_size)
                    if match is None:
                        unmatched_count += 1
                        writer.writerow({**base, "tidehunter_family_key": "NA", "tidehunter_length_bp": "NA", "glocal_identity": "NA",
                                         "tidehunter_chunk_count": "NA", "tidehunter_recurrence_status": "NA", "best_orientation": "NA",
                                         "relation": "no_qualifying_match", "shorter_to_longer_length_ratio": "NA",
                                         "nearest_integer_multiple": "NA", "multiple_relative_error": "NA",
                                         "candidate_pair_count": str(candidate_count),
                                         "postprocess_status": "unmatched_no_tidehunter_support_not_tandemx_only"})
                        continue
                    target = match["target"]
                    assert isinstance(target, dict)
                    permissive_only_count += 1
                    target_chunk_count = int(target["chunk_count"])
                    writer.writerow({**base, "tidehunter_family_key": str(target["family_key"]),
                                     "tidehunter_length_bp": str(target["length"]), "glocal_identity": f"{float(match['identity']):.6f}",
                                     "tidehunter_chunk_count": str(target_chunk_count),
                                     "tidehunter_recurrence_status": ("recurred_all_formal_partitions" if target_chunk_count == FORMAL_PARTITION_COUNT else "not_required_by_frozen_support_screen"),
                                     "best_orientation": str(match["orientation"]), "relation": str(match["relation"]),
                                     "shorter_to_longer_length_ratio": f"{float(match['ratio']):.6f}",
                                     "nearest_integer_multiple": str(match["multiple"]), "multiple_relative_error": f"{float(match['error']):.6f}",
                                     "candidate_pair_count": str(candidate_count), "postprocess_status": "permissive_tidehunter_support"})
        finally:
            store.close()
    receipt.update(
        status="completed_tidehunter_only_postprocess_not_full_comparator_screen",
        recurrent_tandemx_family_count=recurrent_count, exact_tidehunter_supported_family_count=exact_count,
        recurrent_exact_unmatched_tandemx_family_count=query_count, tidehunter_family_count=comparator_count,
        permissive_only_tidehunter_supported_family_count=permissive_only_count,
        permissive_tidehunter_supported_family_count=exact_count + permissive_only_count,
        unmatched_no_tidehunter_support_not_tandemx_only_count=unmatched_count,
        unresolved_family_count=0,
        preflight_file="candidate_index_preflight.json", qgram_preflight_file="qgram_necessary_condition_preflight.json",
        output_sha256={table.name: sha256(table), "candidate_index_preflight.json": sha256(outdir / "candidate_index_preflight.json"),
                       "qgram_necessary_condition_preflight.json": sha256(outdir / "qgram_necessary_condition_preflight.json"), **breakdown_sha256},
        edlib_version=version("edlib"),
        interpretation_boundary="Complete TideHunter postprocess only. It does not establish TandemX-only recovery, novelty, biological recurrence or independent validation.",
    )
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--block-size", type=int, default=500)
    parser.add_argument("--preflight-query-limit", type=int, default=50)
    parser.add_argument("--preflight-comparator-limit", type=int, default=2000)
    parser.add_argument("--max-candidate-pairs", type=int, default=DEFAULT_MAX_CANDIDATE_PAIRS,
                        help="Predeclared maximum exhaustive length-eligible edlib alignments; excess is unresolved.")
    args = parser.parse_args()
    if args.block_size < 1 or args.preflight_query_limit < 1 or args.preflight_comparator_limit < 1 or args.max_candidate_pairs < 1:
        parser.error("block size, preflight limits and candidate-pair budget must be positive")
    receipt = postprocess(args.run_dir, args.outdir, block_size=args.block_size,
                          preflight_query_limit=args.preflight_query_limit,
                          preflight_comparator_limit=args.preflight_comparator_limit,
                          max_candidate_pairs=args.max_candidate_pairs)
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
