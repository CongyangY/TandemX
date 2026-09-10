from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts import validate_recovery_proxy as proxy


F1 = "ACGTTGCACTGATCGTACGATGCTAGTCGATCGTACGTTG"
F2 = "TGACCTAGGCTAACGATTCGAGCTTACCGGATCAGTTCGA"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name, chunks = line[1:], []
        else:
            chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _locked_recovery(tmp_path: Path) -> tuple[Path, Path, Path, str]:
    """Build one locked locus with an old, a read candidate, and two monomers."""
    recovery = tmp_path / "recovery"
    recovery.mkdir()
    old = tmp_path / "old.fa"
    newer = tmp_path / "newer.fa"
    catalog = tmp_path / "monomers.fa"
    logical_span = F1 * 5 + F2 * 5
    _write(catalog, f">family_id=F1;monomer_id=M1\n{F1}\n>family_id=F2;monomer_id=M2\n{F2}\n")
    _write(old, f">chrOld\n{'N' * 100}{logical_span}{'N' * 500}\n")
    # The proxy interval is mapped on the reverse strand.  Extraction must restore
    # the sequence orientation before applying the frozen localizer.
    _write(newer, f">chrNew\n{'N' * 400}{_reverse_complement(logical_span)}{'N' * 500}\n")
    candidates = _write(
        recovery / "recovery_candidates.tsv",
        "locus_id\tfamily_id\tcandidate_read_id\tchromosome\told_between_anchor_start\told_between_anchor_end\n"
        "L0001\tF1\tread-1\tchrOld\t100\t500\n",
    )
    targets = _write(recovery / "recruitment_targets.fa", ">left_anchor\n" + "A" * 100 + "\n>right_anchor\n" + "C" * 100 + "\n")
    _write(recovery / "recovered_sequences.fasta", f">L0001;source=read-1\n{logical_span}\n")
    _write(recovery / "recovery_validation.tsv", "locus_id\trecovery_status\nL0001\tpartially_resolved\n")
    _write(recovery / "recovery_loci.bed", "chrOld\t100\t500\tL0001\n")
    _write(recovery / "recruited_reads.tsv", "locus_id\tread_id\nL0001\tread-1\n")
    mapper = _write(tmp_path / "mock-minimap2", "pinned mapper bytes\n")
    state = {
        "config": {"minimap2": str(mapper), "inputs": {"old_assembly": str(old), "catalog": str(catalog)}},
        "anchors": [
            {"anchor_id": "left_anchor", "locus_id": "L0001", "side": "left", "offset": 0, "eligible": True},
            {"anchor_id": "right_anchor", "locus_id": "L0001", "side": "right", "offset": 0, "eligible": True},
        ],
        "targets_sha256": _sha256(targets),
    }
    _write(recovery / "prepared.json", json.dumps(state))
    manifest = {"inputs": {
        key: {"path": str(path), "sha256": _sha256(path)}
        for key, path in (("old_assembly", old), ("catalog", catalog))
    }, "minimap2_sha256": _sha256(mapper)}
    _write(recovery / "input_manifest.json", json.dumps(manifest))
    _write(recovery / "prepare_lock.json", json.dumps({
        "prepared_sha256": _sha256(recovery / "prepared.json"),
        "input_manifest_sha256": _sha256(recovery / "input_manifest.json"),
    }))
    locked = ("recovery_candidates.tsv", "recovery_validation.tsv", "recovery_loci.bed", "recruited_reads.tsv", "recovered_sequences.fasta")
    _write(recovery / "candidate_lock.json", json.dumps({
        "status": "candidate_generation_complete",
        "outputs": {name: _sha256(recovery / name) for name in locked},
    }))
    return recovery, old, newer, logical_span


def test_validate_requires_complete_unchanged_lock_before_reading_newer(tmp_path: Path) -> None:
    recovery = tmp_path / "recovery"
    recovery.mkdir()
    lock = recovery / "candidate_lock.json"
    _write(lock, json.dumps({"status": "incomplete", "outputs": {}}))
    outdir = tmp_path / "out-incomplete"

    with pytest.raises(ValueError, match="Completed candidate lock"):
        proxy.validate(recovery, tmp_path / "does-not-exist.fa", outdir)
    assert not outdir.exists()

    candidate = _write(recovery / "recovery_candidates.tsv", "x\nlocked\n")
    locked = {"recovery_candidates.tsv": _sha256(candidate)}
    for name in ("recovery_validation.tsv", "recovery_loci.bed", "recruited_reads.tsv", "recovered_sequences.fasta"):
        path = _write(recovery / name, "original\n")
        locked[name] = _sha256(path)
    _write(lock, json.dumps({
        "status": "candidate_generation_complete",
        "outputs": locked,
    }))
    candidate.write_text("x\nchanged\n", encoding="utf-8")
    outdir = tmp_path / "out-mutated"
    with pytest.raises(ValueError, match="Locked recovery output changed: recovery_candidates.tsv"):
        proxy.validate(recovery, tmp_path / "still-does-not-exist.fa", outdir)
    assert not outdir.exists()


@pytest.mark.parametrize(("changed", "message"), (
    ("prepared.json", "Prepared recovery input lock changed"),
    ("recruitment_targets.fa", "Locked recruitment targets changed"),
))
def test_validate_rejects_post_lock_prepared_or_target_tampering_before_newer(
    tmp_path: Path, changed: str, message: str,
) -> None:
    recovery, _old, _newer, _logical_span = _locked_recovery(tmp_path)
    with (recovery / changed).open("a", encoding="utf-8") as handle:
        handle.write("\n")

    outdir = tmp_path / f"out-{changed}"
    with pytest.raises(ValueError, match=message):
        proxy.validate(recovery, tmp_path / "newer-must-not-be-opened.fa", outdir)
    assert not outdir.exists()


def test_validate_reverse_proxy_preserves_inputs_and_uses_frozen_catalog(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    recovery, old, newer, logical_span = _locked_recovery(tmp_path)
    candidate_before, old_before = _sha256(recovery / "recovery_candidates.tsv"), _sha256(old)

    def mock_mapping(_command: list[str], output: Path) -> None:
        # The left/right flanks are in reverse target order on a negative-strand
        # placement: right end (400) precedes left start (800).
        output.write_text(
            "left_anchor\t100\t0\t100\t-\tchrNew\t1300\t800\t900\t100\t100\t60\n"
            "right_anchor\t100\t0\t100\t-\tchrNew\t1300\t300\t400\t100\t100\t60\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(proxy, "run_mapping", mock_mapping)
    outdir = tmp_path / "post-lock-proxy"
    proxy.validate(recovery, newer, outdir)

    assessment = _rows(outdir / "proxy_validation.tsv")[0]
    assert assessment["proxy_status"] == "unique_flank_pair_reference_proxy"
    assert (assessment["newer_start"], assessment["newer_end"]) == ("400", "800")
    assert assessment["newer_span_bp"] == str(len(logical_span))
    # A real call to the frozen locate_record_arrays reports only F1 as the local
    # target despite F2 also occupying the complete frozen catalogue/span.
    assert int(assessment["newer_localized_repeat_bp"]) > 0
    assert int(assessment["newer_localized_repeat_bp"]) < len(logical_span)
    assert {row["family_id"] for row in _rows(outdir / "localized_repeat_intervals.tsv")} == {"F1"}
    assert _fasta(outdir / "comparison_spans.fa")["L0001_newer"] == logical_span
    assert _sha256(recovery / "recovery_candidates.tsv") == candidate_before
    assert _sha256(old) == old_before
