"""Benchmark truth, normalization and metrics must not favor the target tool."""

import json
import math
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from benchmarks.challenge.adapters import build_command, parse_tidehunter, parse_trf, parse_ultra, read_fasta
from benchmarks.challenge.evaluate import circular_identity, maximum_matching, score_arrays, score_families, wilson_interval
from benchmarks.challenge.run import run_process, summarize, source_manifest
from benchmarks.challenge.schema import ArrayRecord, read_table
from benchmarks.challenge.simulate import Scenario, generate_dataset, mutate, reverse_complement


def test_challenge_reproducible_indels_and_reverse_complements(tmp_path: Path) -> None:
    config = Scenario("test", period=37, copies=6, read_count=20, read_length=1000,
                      insertion_rate=0.02, deletion_rate=0.01)
    a = generate_dataset(config, 11, tmp_path / "a")
    b = generate_dataset(config, 11, tmp_path / "b")
    assert a == b
    fasta = read_fasta(tmp_path / "a" / "reads.fa")
    reads = read_table(tmp_path / "a" / "truth_reads.tsv")
    assert {row["strand"] for row in reads} == {"+", "-"}
    assert all(len(fasta[row["read_id"]]) == int(row["length_bp"]) for row in reads)
    arrays = read_table(tmp_path / "a" / "truth_arrays.tsv")
    assert len(arrays) == 14
    assert all(0 <= int(r["start"]) < int(r["end"]) <= len(fasta[r["read_id"]]) for r in arrays)
    with pytest.raises(FileExistsError):
        generate_dataset(config, 11, tmp_path / "a")


def test_clean_truth_is_observable_in_read_coordinates(tmp_path: Path) -> None:
    config = Scenario("clean", period=37, copies=4, read_count=12, read_length=1000)
    generate_dataset(config, 23, tmp_path / "data")
    fasta = read_fasta(tmp_path / "data" / "reads.fa")
    for row in read_table(tmp_path / "data" / "truth_arrays.tsv"):
        sequence = fasta[row["read_id"]][int(row["start"]):int(row["end"])]
        assert len(sequence) == 4 * 37
        assert sequence[:37] * 4 == sequence
        assert circular_identity(sequence[:37], row["sequence"]) == 1


@pytest.mark.parametrize("changes", [{"period": 0}, {"read_count": 0}, {"substitution_rate": -1},
                                     {"copies": 999}, {"deletion_rate": 1}, {"name": "../oops"}])
def test_invalid_simulations_fail(changes: dict) -> None:
    with pytest.raises(ValueError):
        replace(Scenario("valid"), **changes).validate()


def test_mutation_extremes() -> None:
    import random
    assert mutate("ACGT", random.Random(1), deletion=1) == ""
    assert len(mutate("ACGT", random.Random(1), ins=1)) == 8
    assert all(a != b for a, b in zip("ACGT", mutate("ACGT", random.Random(1), sub=1)))


def test_maximum_matching_does_not_use_greedy_first_hit() -> None:
    assert maximum_matching([[0, 1], [0]]) == {0: 1, 1: 0}


def test_array_metrics_penalize_duplicates_wrong_periods_and_false_positives() -> None:
    truth = [ArrayRecord("p", 100, 600, 100), ArrayRecord("p", 800, 1300, 100)]
    predicted = [truth[0], truth[0], ArrayRecord("p", 800, 1300, 200), ArrayRecord("n", 0, 500, 100)]
    metrics, details = score_arrays(predicted, truth, {"p": 1500, "n": 1500})
    assert metrics["array_recall"] == 0.5
    assert metrics["array_precision"] == 0.25
    assert metrics["read_detection_recall"] == 1
    assert metrics["read_detection_precision"] == 0.5
    assert metrics["negative_read_call_rate"] == 1
    assert sum(row["status"] == "matched" for row in details) == 1


def test_bad_intervals_and_unknown_reads_raise() -> None:
    with pytest.raises(ValueError):
        ArrayRecord("r", 5, 5, 100)
    with pytest.raises(ValueError):
        score_arrays([ArrayRecord("unknown", 0, 50, 10)], [], {"r": 100})
    with pytest.raises(ValueError):
        score_arrays([ArrayRecord("r", 0, 150, 10)], [], {"r": 100})


def test_absent_denominators_are_na_not_perfect_or_zero_recall() -> None:
    metrics, _ = score_arrays([], [], {"negative": 100})
    assert math.isnan(metrics["array_recall"])
    assert math.isnan(metrics["array_precision"])
    assert metrics["negative_read_call_rate"] == 0
    lower, upper = wilson_interval(0, 100)
    assert lower == 0
    assert 0.036 < upper < 0.038


def test_independent_family_matching_and_orientation() -> None:
    seq = "ACCCGTATTGACG"
    assert circular_identity(seq, reverse_complement(seq[3:] + seq[:3])) == 1
    assert circular_identity(seq, "T" * len(seq)) < 0.5
    assert circular_identity(seq, seq + "A") == 0
    # A single consensus cannot count as two distinct recovered families.
    metrics, _ = score_families([seq, seq[3:] + seq[:3], reverse_complement(seq)], {"a": seq, "b": seq})
    assert metrics["sequence_family_recall"] == 0.5


def test_strict_external_parsers_and_one_based_conversion(tmp_path: Path) -> None:
    tide = tmp_path / "tide.tsv"
    tide.write_text("r rep0 4 1000 101 600 5 99 0 101,106 ACGTA\n")
    trf = tmp_path / "trf.txt"
    trf.write_text("@r\n101 600 5 100 5 99 0 100 25 25 25 25 2.0 ACGTA ACGTA\n")
    assert parse_tidehunter(tide) == parse_trf(trf) == [ArrayRecord("r", 100, 600, 5, "ACGTA")]
    tide.write_text("broken row\n")
    with pytest.raises(ValueError):
        parse_tidehunter(tide)
    with pytest.raises(FileNotFoundError):
        parse_trf(tmp_path / "missing")


def test_commands_never_receive_truth(tmp_path: Path) -> None:
    for tool in ("tandemx", "trf", "tidehunter", "ultra"):
        command, _ = build_command(tool, tool, tmp_path / "reads.fa", tmp_path, 30, 1000, 100)
        assert not any("truth" in value for value in command)
        if tool == "tandemx":
            assert command[command.index("--min-repeat-span") + 1] == "100"
            assert command[command.index("--threads") + 1] == "1"

    threaded = {
        tool: build_command(tool, tool, tmp_path / "reads.fa", tmp_path, 30, 1000, 100, 4)[0]
        for tool in ("tandemx", "tidehunter", "ultra")
    }
    assert threaded["tandemx"][threaded["tandemx"].index("--threads") + 1] == "4"
    assert threaded["tidehunter"][threaded["tidehunter"].index("-t") + 1] == "4"
    assert threaded["ultra"][threaded["ultra"].index("-t") + 1] == "4"
    with pytest.raises(ValueError, match="Threads"):
        build_command("tandemx", "tandemx", tmp_path / "reads.fa", tmp_path, 30, 1000, 100, 0)


def test_ultra_half_open_coordinates_and_unknown_consensus(tmp_path):
    path = tmp_path / "ultra.tsv"
    path.write_text("SeqID\tStart\tEnd\tPeriod\tScore\tConsensus\n"
                    "r\t100\t600\t5\t12.5\tACG*A\n"
                    "s\t0\t100\t10\t10\t.\n")
    assert parse_ultra(path) == [ArrayRecord("r", 100, 600, 5, "ACGNA"), ArrayRecord("s", 0, 100, 10)]
    path.write_text("bad\tformat\n")
    with pytest.raises(ValueError):
        parse_ultra(path)


def test_source_snapshot_preserves_dirty_or_unversioned_code(tmp_path):
    root = tmp_path / "source"
    (root / "tandemx").mkdir(parents=True)
    (root / "tandemx" / "__init__.py").write_text("")
    source = root / "tandemx" / "sample.py"
    source.write_text("value = 42\n")
    snapshot = tmp_path / "snapshot"
    manifest = source_manifest(root, snapshot)
    source.write_text("value = 43\n")
    assert (snapshot / "tandemx" / "sample.py").read_text() == "value = 42\n"
    assert manifest["file_hashes"]["tandemx/sample.py"]
    assert manifest["git_head"] is None
    assert manifest["revision_warning"] == "not_a_git_checkout_use_source_digest"
    run = run_process([sys.executable, "-c", "from tandemx.sample import value; print(value)"],
                      tmp_path / "out", tmp_path / "err", 5,
                      {**os.environ, "PYTHONPATH": str(snapshot)}, snapshot)
    assert run["exit_code"] == 0, (tmp_path / "err").read_text()
    assert (tmp_path / "out").read_text().strip() == "42"


def test_source_manifest_marks_tracked_source_changes_against_git_head(tmp_path):
    root = tmp_path / "repo"
    source = root / "tandemx" / "sample.py"
    source.parent.mkdir(parents=True)
    source.write_text("value = 1\n")
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    subprocess.run(["git", "add", "tandemx/sample.py"], cwd=root, check=True)
    subprocess.run([
        "git", "-c", "user.name=TandemX test", "-c", "user.email=test@example.invalid",
        "commit", "-q", "-m", "fixture"
    ], cwd=root, check=True)
    clean = source_manifest(root)
    assert clean["git_head"]
    assert clean["revision_warning"] is None
    source.write_text("value = 2\n")
    dirty = source_manifest(root)
    assert dirty["git_head"] == clean["git_head"]
    assert dirty["revision_warning"] == "worktree_differs_from_git_head_use_source_digest"
    source.write_text("value = 1\n")
    (root / ".gitignore").write_text("*.so\n")
    native = root / "tandemx" / "generated.so"
    native.write_bytes(b"native build")
    generated = source_manifest(root)
    assert generated["revision_warning"] == "worktree_differs_from_git_head_use_source_digest"
    assert generated["git_untracked_snapshot_files"] == ["tandemx/generated.so"]


def test_process_failure_and_timeout_are_observable(tmp_path: Path) -> None:
    result = run_process([sys.executable, "-c", "raise SystemExit(3)"], tmp_path / "out", tmp_path / "err", 5)
    assert result["exit_code"] == 3
    result = run_process([sys.executable, "-c", "import time; time.sleep(5)"], tmp_path / "out", tmp_path / "err", 0.1)
    assert result["timed_out"] and result["exit_code"] != 0


def test_failed_runs_never_enter_accuracy_summary() -> None:
    failed = dict(scenario="a", dataset_id="a_s1", tool="t", seed=1, split="development", status="failed")
    row = summarize([failed])[0]
    assert row["array_recall"] == "NA"
    assert row["median_runtime_seconds"] == "NA"
    assert row["successful_runs"] == 0
