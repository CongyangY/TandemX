import csv
import json
import sqlite3

from benchmarks.scripts.postprocess_ysd56_tidehunter_permissive import (
    candidate_pair_estimate,
    qgram_necessary_condition,
    postprocess,
)
from benchmarks.scripts.run_real_30x_cross_tool import _append_family


def _run_dir(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    database = run_dir / "normalization.sqlite"
    connection = sqlite3.connect(database)
    connection.execute("CREATE TABLE families(tool TEXT NOT NULL, family_key TEXT NOT NULL, canonical_unit TEXT NOT NULL, primitive_period_bp INTEGER NOT NULL, source_record_count INTEGER NOT NULL, chunk_count INTEGER NOT NULL, source_kind TEXT NOT NULL, PRIMARY KEY(tool, family_key)) WITHOUT ROWID")
    connection.execute("CREATE TABLE family_chunks(tool TEXT NOT NULL, family_key TEXT NOT NULL, chunk_id TEXT NOT NULL, PRIMARY KEY(tool, family_key, chunk_id)) WITHOUT ROWID")
    query = "ACGTTGCACTGATCGTAGCTAACGTTGCAA"
    direct = "ACGTTGCACTGATCGTAGCTAACGTTGCGC"
    related = "TTTCCCAAAGGGTTTCCCAAAGGGTTTCCCA"
    unrelated = "GATCTAGCGTACCGATGCTAGCTAGCGTAC"
    for chunk in ("p1", "p2", "p3"):
        _append_family(connection, tool="tandemx", sequence=query, source_kind="catalogue_family", chunk_id=chunk)
        _append_family(connection, tool="tandemx", sequence=related, source_kind="catalogue_family", chunk_id=chunk)
        _append_family(connection, tool="tandemx", sequence=unrelated, source_kind="catalogue_family", chunk_id=chunk)
    _append_family(connection, tool="tidehunter", sequence=query, source_kind="array_consensus", chunk_id="p1")
    _append_family(connection, tool="tidehunter", sequence=direct, source_kind="array_consensus", chunk_id="p1")
    _append_family(connection, tool="tidehunter", sequence=related + "GCGCATATGCGCATATGCGCATATGCGCAT", source_kind="array_consensus", chunk_id="p1")
    connection.commit(); connection.close()
    (run_dir / "run_manifest.json").write_text(json.dumps({
        "schema_version": 1, "input_fasta_sha256": "a" * 64, "semantic_sha256": "b" * 64,
        "partitions": [
            {"partition_id": partition, "sha256": "c" * 64, "read_count": 10, "total_bases": 100}
            for partition in ("partition_00001", "partition_00002", "partition_00003")
        ],
    }) + "\n")
    with (run_dir / "execution_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("tool", "partition_id", "status", "normalization"), delimiter="\t")
        writer.writeheader()
        for partition in ("partition_00001", "partition_00002", "partition_00003"):
            writer.writerow({"tool": "tidehunter", "partition_id": partition, "status": "completed", "normalization": "ok"})
    return run_dir


def test_postprocess_reports_exact_and_permissive_support_without_seed_negative_calls(tmp_path):
    run_dir = _run_dir(tmp_path)
    receipt = postprocess(run_dir, tmp_path / "audit", block_size=1, preflight_query_limit=3, preflight_comparator_limit=3)
    assert receipt["status"] == "completed_tidehunter_only_postprocess_not_full_comparator_screen"
    assert receipt["exact_tidehunter_supported_family_count"] == 1
    assert receipt["recurrent_tandemx_family_count"] == 3
    assert receipt["recurrent_exact_unmatched_tandemx_family_count"] == 2
    assert receipt["permissive_only_tidehunter_supported_family_count"] == 1
    assert receipt["permissive_tidehunter_supported_family_count"] == 2
    assert receipt["unresolved_family_count"] == 0
    with (tmp_path / "audit" / "ysd56_tandemx_tidehunter_permissive_matches.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert any(row["postprocess_status"] == "unmatched_no_tidehunter_support_not_tandemx_only" for row in rows)
    preflight = json.loads((tmp_path / "audit" / "candidate_index_preflight.json").read_text())
    assert preflight["interpretation"].startswith("preflight only")
    assert preflight["screen_population"].startswith("cross_partition_recurrent")


def test_pair_budget_is_unresolved_not_negative_and_nonrecurrent_rows_are_excluded(tmp_path):
    run_dir = _run_dir(tmp_path)
    receipt = postprocess(run_dir, tmp_path / "budget", max_candidate_pairs=1)
    assert receipt["status"] == "blocked_candidate_pair_estimate_exceeds_gate"
    assert receipt["permissive_tidehunter_supported_family_count"] == "NA"
    assert receipt["unresolved_family_count"] == 2
    with (tmp_path / "budget" / "ysd56_tandemx_tidehunter_permissive_matches.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert {row["postprocess_status"] for row in rows} == {"technical_unresolved_candidate_pair_budget"}


def test_length_pair_estimate_counts_direct_multiple_and_union_without_overlap():
    query_histogram = {40: 2}
    comparator_histogram = {40: 3, 43: 5, 80: 7, 82: 11}
    # 40/40 and 40/43 are direct; 40/80 and 40/82 meet the two-period branch.
    assert candidate_pair_estimate(query_histogram, comparator_histogram) == 2 * (3 + 5 + 7 + 11)


def test_lossless_qgram_necessary_condition_keeps_direct_and_period_multiple_qualifiers():
    direct = "ACGTTGCA" * 5
    direct_mutated = direct[:4] + "A" + direct[5:20] + "T" + direct[21:]
    assert qgram_necessary_condition(direct, direct_mutated, 0.90)
    short = "ACGTTGCA" * 5
    assert qgram_necessary_condition(short, short * 2, 0.80)


def test_postprocess_writes_explicit_remount_block_without_invented_counts(tmp_path):
    receipt = postprocess(tmp_path / "missing_run", tmp_path / "blocked")
    assert receipt["status"] == "blocked_missing_post_remount_source_outputs"
    assert receipt["exact_tidehunter_supported_family_count"] == "NA"
    assert receipt["permissive_tidehunter_supported_family_count"] == "NA"
    saved = json.loads((tmp_path / "blocked" / "receipt.json").read_text())
    assert saved["unresolved_family_count"] == "NA"
