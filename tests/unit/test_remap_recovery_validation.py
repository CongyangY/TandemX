from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts import remap_recovery_validation as remap_validation


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def _inputs(tmp_path: Path, *, include_candidate: bool = False) -> tuple[Path, Path, Path]:
    proxy, reads = tmp_path / "proxy", tmp_path / "reads"
    spans = _write(proxy / "comparison_spans.fa", ">L0001_old\n" + "ACGT" * 25 + "\n")
    _write(proxy / "proxy_completion.json", json.dumps({
        "status": "posthoc_proxy_complete", "outputs": {"comparison_spans.fa": _sha256(spans)},
    }))
    read_records = (
        ">independent-1\n" + "ACGT" * 25
        + "\n>independent-2\n" + "TGCA" * 25
        + "\n>low-coverage\n" + "ACGT" * 25
        + "\n>low-identity\n" + "TGCA" * 25 + "\n"
    )
    if include_candidate:
        read_records += ">candidate-source\n" + "ACGT" * 25 + "\n"
    fasta = _write(reads / "reads.fa", read_records)
    _write(reads / "receipt.json", json.dumps({
        "complete": True,
        "candidate_read_id_excluded": "candidate-source",
        "extracted_independent_read_ids": 5 if include_candidate else 4,
        "outputs": {"reads.fa": _sha256(fasta)},
    }))
    mapper = _write(tmp_path / "minimap2", "pinned mapper\n")
    return proxy, reads, mapper


def test_remap_counts_distinct_reads_only_above_coverage_and_identity_thresholds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    proxy, reads, mapper = _inputs(tmp_path)

    def mock_mapping(_command: list[str], output: Path) -> None:
        # query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, block, mapq
        output.write_text(
            "L0001_old\t100\t0\t100\t+\tindependent-1\t100\t0\t100\t100\t100\t60\n"
            "L0001_old\t100\t0\t100\t+\tindependent-1\t100\t0\t100\t100\t100\t60\n"
            "L0001_old\t100\t0\t80\t+\tlow-coverage\t100\t0\t80\t80\t80\t60\n"
            "L0001_old\t100\t0\t100\t+\tlow-identity\t100\t0\t100\t97\t100\t60\n"
            "L0001_old\t100\t0\t90\t+\tindependent-2\t100\t0\t90\t90\t90\t60\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(remap_validation, "run_mapping", mock_mapping)
    outdir = tmp_path / "remapped"
    remap_validation.remap(proxy, reads, mapper, outdir)

    row = _row(outdir / "read_remapping.tsv")
    assert row["sequence_id"] == "L0001_old"
    assert row["independently_supporting_reads"] == "2"
    assert row["supplied_independent_reads"] == "4"


def test_remap_rejects_candidate_source_read_before_mapping(tmp_path: Path) -> None:
    proxy, reads, mapper = _inputs(tmp_path, include_candidate=True)
    outdir = tmp_path / "source-included"

    with pytest.raises(ValueError, match="Independent read membership is invalid"):
        remap_validation.remap(proxy, reads, mapper, outdir)
    assert not outdir.exists()


def test_remap_rejects_paf_target_outside_enrolled_reads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    proxy, reads, mapper = _inputs(tmp_path)

    def mock_mapping(_command: list[str], output: Path) -> None:
        output.write_text(
            "L0001_old\t100\t0\t100\t+\tunknown-read\t100\t0\t100\t100\t100\t60\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(remap_validation, "run_mapping", mock_mapping)
    outdir = tmp_path / "unknown-target"
    with pytest.raises(ValueError, match="outside the enrolled subset"):
        remap_validation.remap(proxy, reads, mapper, outdir)
    assert not (outdir / "read_remapping.tsv").exists()


@pytest.mark.parametrize("kind", ("incomplete_receipt", "proxy_hash", "read_hash"))
def test_remap_rejects_bad_receipts_or_locked_hashes_before_mapping(tmp_path: Path, kind: str) -> None:
    proxy, reads, mapper = _inputs(tmp_path)
    if kind == "incomplete_receipt":
        receipt = json.loads((reads / "receipt.json").read_text())
        receipt["complete"] = False
        _write(reads / "receipt.json", json.dumps(receipt))
        message = "Completed proxy and independent-read receipts required"
    elif kind == "proxy_hash":
        with (proxy / "comparison_spans.fa").open("a", encoding="utf-8") as handle:
            handle.write("A\n")
        message = "Remapping input receipt/hash mismatch"
    else:
        with (reads / "reads.fa").open("a", encoding="utf-8") as handle:
            handle.write("A\n")
        message = "Remapping input receipt/hash mismatch"

    outdir = tmp_path / kind
    with pytest.raises(ValueError, match=message):
        remap_validation.remap(proxy, reads, mapper, outdir)
    assert not outdir.exists()
