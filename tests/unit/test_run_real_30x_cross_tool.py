import json
import sqlite3

from benchmarks.scripts.run_real_30x_cross_tool import (
    _append_family,
    _build_permissive_candidate_index,
    _candidate_rows,
    build_frozen_tandemx_command,
    create_chunks,
    classify_permissive_match,
    normalized_family,
    primitive_unit,
    screen_permissive_matches,
    write_exact_comparator_support,
    write_cross_partition_recurrence,
    write_normalized_families,
    write_exact_match_screen,
    run,
)


def _family_database(tmp_path):
    database = tmp_path / "normalization.sqlite"
    connection = sqlite3.connect(database)
    connection.execute("CREATE TABLE families(tool TEXT NOT NULL, family_key TEXT NOT NULL, canonical_unit TEXT NOT NULL, primitive_period_bp INTEGER NOT NULL, source_record_count INTEGER NOT NULL, chunk_count INTEGER NOT NULL, source_kind TEXT NOT NULL, PRIMARY KEY(tool, family_key)) WITHOUT ROWID")
    connection.execute("CREATE TABLE family_chunks(tool TEXT NOT NULL, family_key TEXT NOT NULL, chunk_id TEXT NOT NULL, PRIMARY KEY(tool, family_key, chunk_id)) WITHOUT ROWID")
    return database, connection


def test_family_normalization_accepts_rotation_reverse_complement_and_period_multiple():
    assert primitive_unit("ACGACG") == "ACG"
    assert normalized_family("ACGACG") == normalized_family("CGACGA")
    assert normalized_family("ACG") == normalized_family("CGT")


def test_failed_tool_makes_all_tandemx_exact_screen_rows_unresolved(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="c1")
    _append_family(connection, tool="trf", sequence="CGACGA", source_kind="array_consensus", chunk_id="c1")
    connection.commit(); connection.close()
    rows = write_exact_match_screen(database, tmp_path / "exact_screen.tsv", complete=False)
    assert rows[0]["matched_other_tools"] == "trf"
    assert rows[0]["screen_status"] == "unresolved"
    assert "failed_or_timed_out" in rows[0]["warning"]
    write_normalized_families(database, tmp_path / "all.tsv")
    assert "canonical_unit" in (tmp_path / "all.tsv").read_text()


def test_exact_unmatched_family_requires_permissive_screen_never_unique(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="c1")
    connection.commit(); connection.close()
    rows = write_exact_match_screen(database, tmp_path / "exact_screen.tsv", complete=True)
    assert rows[0]["screen_status"] == "requires_permissive_family_match_screen"
    assert "unique" not in rows[0]["screen_status"]
    assert "nonexact_family_comparison" in rows[0]["warning"]


def test_exact_screen_can_report_one_completed_comparator_without_erasing_another_timeout(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="c1")
    _append_family(connection, tool="tidehunter", sequence="CGACGA", source_kind="array_consensus", chunk_id="c1")
    connection.commit(); connection.close()
    rows = write_exact_match_screen(
        database, tmp_path / "exact_screen.tsv", True,
        comparator_tools=("tidehunter",), unavailable_comparator_tools=("trf",),
    )
    assert rows[0]["matched_other_tools"] == "tidehunter"
    assert rows[0]["screen_status"] == "matched_by_exact_normalization_preliminary"
    assert "unavailable_comparator_tools=trf" in rows[0]["warning"]


def test_cross_partition_recurrence_is_separate_from_comparator_support(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="p1")
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="p2")
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="p3")
    connection.commit(); connection.close()
    rows = write_cross_partition_recurrence(database, tmp_path / "recurrence.tsv", 3)
    assert rows[0]["cross_partition_recurrence_status"] == "recurred_all_formal_partitions"


def test_exact_comparator_support_keeps_exact_unmatched_families_unresolved(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="p1")
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="p2")
    _append_family(connection, tool="tandemx", sequence="ACGACG", source_kind="catalogue_family", chunk_id="p3")
    connection.commit(); connection.close()
    rows = write_exact_comparator_support(database, tmp_path / "support.tsv", 3, "tidehunter")
    assert rows[0]["cross_partition_recurrence_status"] == "recurred_all_formal_partitions"
    assert rows[0]["exact_comparator_support_status"] == "exact_unmatched_requires_permissive_screen"
    assert rows[0]["comparator_chunk_count"] == "NA"


def test_permissive_screen_accepts_rotation_reverse_complement_direct_match(tmp_path):
    database, connection = _family_database(tmp_path)
    query = "ACGTTGCACTGATCGTAGCTAACGTTGCAA"
    reverse_complement_rotation = "TTGCAACGTTAGCTACGATCAGTGCAACGT"
    _append_family(connection, tool="tandemx", sequence=query, source_kind="catalogue_family", chunk_id="c1")
    _append_family(connection, tool="trf", sequence=reverse_complement_rotation, source_kind="array_consensus", chunk_id="c1")
    connection.commit(); connection.close()
    rows, receipt = screen_permissive_matches(database, tmp_path / "permissive.tsv", complete=True)
    assert rows[0]["permissive_match_status"] == "matched_by_permissive_screen"
    assert rows[0]["relation"] == "direct_glocal"
    assert float(rows[0]["glocal_identity"]) >= .80
    assert receipt["candidate_index_saturated_queries"] == 0


def test_permissive_screen_accepts_integer_period_multiple(tmp_path):
    database, connection = _family_database(tmp_path)
    query = "ACGTTGCACTGATCGTAGCTAACGTTGCAA"
    related_longer = query + "TTTTGCGCGATATATCGCGCTTTTGGGCCC"
    _append_family(connection, tool="tandemx", sequence=query, source_kind="catalogue_family", chunk_id="c1")
    _append_family(connection, tool="tidehunter", sequence=related_longer, source_kind="array_consensus", chunk_id="c1")
    connection.commit(); connection.close()
    rows, _ = screen_permissive_matches(database, tmp_path / "permissive.tsv", complete=True)
    assert rows[0]["relation"] == "integer_period_multiple"
    assert rows[0]["nearest_integer_multiple"] == "2"
    assert float(rows[0]["glocal_identity"]) >= .80


def test_permissive_threshold_boundaries_and_no_match_state():
    state, multiple, error = classify_permissive_match(.90, .90)
    assert (state, multiple) == ("direct_glocal", 1)
    assert abs(error - (1 / .90 - 1)) < 1e-12
    assert classify_permissive_match(.90, min(90, 100) / max(90, 100))[0] == "direct_glocal"
    assert classify_permissive_match(.90, min(100, 90) / max(100, 90))[0] == "direct_glocal"
    assert classify_permissive_match(.90, .899999)[0] == "no_qualifying_match"
    state, multiple, error = classify_permissive_match(.80, 1 / 2.10)
    assert (state, multiple) == ("integer_period_multiple", 2)
    assert abs(error - .05) < 1e-12
    assert classify_permissive_match(.80, 1 / 2.100001)[0] == "no_qualifying_match"


def test_eighty_percent_length_30_candidate_retains_a_four_mer_seed(tmp_path):
    database, connection = _family_database(tmp_path)
    query = "ACGTTGCACTGATCGTAGCTAACGTTGCAA"
    mutated = "TTTTTT" + query[6:]
    related_longer = mutated + "GCGCATATGCGCATATGCGCATATGCGCAT"
    _append_family(connection, tool="tidehunter", sequence=related_longer, source_kind="array_consensus", chunk_id="c1")
    connection.commit()
    lengths = _build_permissive_candidate_index(connection)
    candidates, saturated = _candidate_rows(connection, "tidehunter", query, lengths, 10_000)
    connection.close()
    assert not saturated
    assert len(candidates) == 1


def test_candidate_index_uses_only_each_family_own_seed_length(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="trf", sequence="ACGTTGCACTGATCGTAGCTAACGTTGCAA", source_kind="array_consensus", chunk_id="c1")
    connection.commit()
    _build_permissive_candidate_index(connection)
    seed_lengths = [row[0] for row in connection.execute("SELECT DISTINCT seed_length FROM permissive_family_seeds")]
    connection.close()
    assert seed_lengths == [4]


def test_frozen_tandemx_command_matches_ysd56_production_parameters(tmp_path):
    command = build_frozen_tandemx_command(tmp_path / "partition.fa", tmp_path / "discover", threads=4)
    assert command == [
        "python", "-m", "tandemx.cli", "discover", "--reads", str(tmp_path / "partition.fa"),
        "--outdir", str(tmp_path / "discover"), "--discovery-method", "cascade",
        "--clustering-method", "sequence", "--cluster-identity", ".95", "--family-audit", "related",
        "--min-period", "30", "--max-period", "1000", "--top-periods", "5",
        "--min-support-reads", "5", "--min-repeat-span", "100", "--kmer-backend", "rust",
        "--threads", "4", "--no-progress",
    ]


def test_permissive_no_match_is_eligible_only_after_complete_screen(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGTTGCACTGATCGTAGCTAACGTTGCAA", source_kind="catalogue_family", chunk_id="c1")
    _append_family(connection, tool="trf", sequence="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGG", source_kind="array_consensus", chunk_id="c1")
    connection.commit(); connection.close()
    rows, _ = screen_permissive_matches(database, tmp_path / "permissive.tsv", complete=True)
    assert rows[0]["permissive_match_status"] == "eligible_for_recurrence_and_biological_triage"
    assert "novel" not in " ".join(rows[0].values()).lower()


def test_permissive_tool_failure_is_technical_unresolved(tmp_path):
    database, connection = _family_database(tmp_path)
    _append_family(connection, tool="tandemx", sequence="ACGTTGCACTGATCGTAGCTAACGTTGCAA", source_kind="catalogue_family", chunk_id="c1")
    connection.commit(); connection.close()
    rows, _ = screen_permissive_matches(database, tmp_path / "permissive.tsv", complete=False)
    assert rows[0]["permissive_match_status"] == "technical_unresolved"


def test_chunk_plan_keeps_each_read_intact_and_hashes_the_common_input(tmp_path):
    reads = tmp_path / "reads.fa"
    reads.write_text(">r1 description\nACGT\n>r2\nAACCGG\n")
    result = create_chunks(reads, tmp_path / "chunks", tmp_path / "reads.sqlite", 6)
    assert result["read_count"] == 2 and result["total_bases"] == 10
    assert [chunk["read_count"] for chunk in result["chunks"]] == [1, 1]
    assert (tmp_path / "chunks" / "chunk_00001.fa").read_text() == ">r1\nACGT\n"


def test_runner_retains_complete_execution_receipts_with_empty_tandemx_catalogue(tmp_path):
    reads = tmp_path / "reads.fa"
    unit = "ACGTTGCACTGATCGTAGCTAACGTTGCAA"
    reads.write_text("".join(f">r{index}\n{unit * 4}\n" for index in range(1, 4)))
    trf = tmp_path / "trf"
    trf.write_text("#!/bin/sh\nid=$(sed -n '1s/>//p' \"$1\")\nprintf '@%s\\n1 120 30 0 0 0 0 0 0 0 0 0 0 ACGTTGCACTGATCGTAGCTAACGTTGCAA\\n' \"$id\"\n")
    tidehunter = tmp_path / "tidehunter"
    tidehunter.write_text("#!/bin/sh\nwhile [ \"$1\" != '-o' ]; do shift; done\nshift\nout=$1\nshift\nfor arg in \"$@\"; do reads=$arg; done\nid=$(sed -n '1s/>//p' \"$reads\")\nprintf '%s R1 4 120 1 120 30 99 4 1,31,61,91 ACGTTGCACTGATCGTAGCTAACGTTGCAA\\n' \"$id\" > \"$out\"\n")
    trf.chmod(0o755); tidehunter.chmod(0o755)
    result = run(reads, tmp_path / "out", trf, tidehunter, 60, 120, 1, assembly_span_bp=12)
    assert result["status"] == "completed_permissive_screen"
    assert result["exact_match_screen_status"] == "complete_preliminary_only"
    assert result["permissive_family_match_status"] == "complete"
    manifest = json.loads((tmp_path / "out" / "run_manifest.json").read_text())
    assert [partition["observed_depth"] for partition in manifest["partitions"]] == [10.0, 10.0, 10.0]
    command = json.loads((tmp_path / "out" / "runs" / "tandemx" / "partition_00001" / "execution.json").read_text())["command"]
    assert command[3:] == build_frozen_tandemx_command(
        tmp_path / "out" / "partitions" / "partition_00001.fa",
        tmp_path / "out" / "runs" / "tandemx" / "partition_00001" / "discover",
        threads=1,
    )[3:]
    summary = (tmp_path / "out" / "execution_summary.tsv").read_text()
    assert summary.count("\tcompleted\t0\tFalse") == 9
