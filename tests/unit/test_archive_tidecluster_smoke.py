from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_tidecluster_smoke import archive, validate


IMAGE = "sha256:" + "a" * 64


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _time_text(rss: int = 100) -> str:
    return (
        "\tUser time (seconds): 1.0\n\tSystem time (seconds): 0.5\n"
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:01.50\n"
        f"\tMaximum resident set size (kbytes): {rss}\n\tExit status: 0\n"
    )


def _fixture(tmp_path: Path) -> tuple[Path, Path]:
    source = tmp_path / "run"
    container = tmp_path / "container"
    (source / "tc_consensus").mkdir(parents=True)
    (source / "normalized").mkdir()
    (source / "profile").mkdir()
    container.mkdir()
    for name in ("tc_tidehunter.gff3", "tc_clustering.gff3"):
        (source / name).write_text("chr1\tx\ttandem_repeat\t1\t10\t.\t.\t.\tName=x\n" * 3)
    (source / "tc_cmd_args.json").write_text("{}\n")
    (source / "tc_consensus/consensus_sequences_all.fasta").write_text(">x\nACGT\n")
    (source / "tidehunter.gnu_time.txt").write_text(_time_text())
    (source / "clustering.gnu_time.txt").write_text(_time_text(200))
    (source / "profile/receipt.json").write_text(json.dumps({
        "complete": True, "requested_stage_count": 2, "completed_stage_count": 2,
        "rss_scope": "aggregate_live_process_tree",
    }))
    _write_tsv(source / "profile/stages.tsv", [
        {"stage": "tidehunter", "exit_code": 0},
        {"stage": "clustering", "exit_code": 0},
    ])
    for name in ("tidehunter.stdout.log", "tidehunter.stderr.log", "clustering.stdout.log", "clustering.stderr.log"):
        (source / "profile" / name).write_text("")
    metrics = {
        "truth_array_count": 3, "predicted_array_count": 3, "matched_array_count": 3,
        "array_recall": 1.0, "array_precision": 1.0, "matched_period_mae_bp": 0.0,
        "truth_family_count": 3, "recovered_family_count": 3, "sequence_family_recall": 1.0,
        "base_union_recall": 1.0, "base_union_precision": 0.99,
        "warning": "known_planted_assembly_truth;smoke_test_not_publication_scale",
    }
    (source / "normalized/metrics.json").write_text(json.dumps(metrics))
    for name in ("normalized_arrays.tsv", "matches.tsv", "family_recovery.tsv"):
        (source / "normalized" / name).write_text("x\n")
    external = tmp_path / "input.fa"
    external.write_text(">x\nACGT\n")
    outputs = {name: digest_file(source / "normalized" / name) for name in (
        "normalized_arrays.tsv", "matches.tsv", "family_recovery.tsv", "metrics.json"
    )}
    receipt = {"complete": True, "inputs": {str(external): {"sha256": digest_file(external), "bytes": external.stat().st_size}}, "outputs": outputs}
    (source / "normalized/receipt.json").write_text(json.dumps(receipt))
    (container / "Dockerfile").write_text("FROM example\n")
    (container / "README.md").write_text("container\n")
    return source, container


def test_archive_tidecluster_smoke_checks_and_copies_compact_evidence(tmp_path: Path) -> None:
    source, container = _fixture(tmp_path)
    result = archive(source, container, IMAGE, tmp_path / "archive")
    assert result["complete"] is True
    assert len(result["files"]) == 20
    summary = json.loads((tmp_path / "archive/summary.json").read_text())
    assert summary["internal_gnu_time"]["clustering"]["maximum_rss_kb"] == 200


def test_validate_tidecluster_smoke_rejects_changed_evaluation_input(tmp_path: Path) -> None:
    source, container = _fixture(tmp_path)
    external = tmp_path / "input.fa"
    external.write_text(">x\nTGCA\n")
    with pytest.raises(ValueError, match="evaluation input changed"):
        validate(source, container, IMAGE)
