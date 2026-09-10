from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts.render_existing_family_report import render_existing_family_report


def _write(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8")


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inputs(root: Path) -> tuple[Path, Path, Path]:
    discover = root / "discover"
    _write(discover / "families.tsv", "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\nF1\tM1\t10\th\t0.5\t2\t20\t1\tfalse\thigh\t\n")
    _write(discover / "monomers.fa", ">family_id=F1;monomer_id=M1;length_bp=10;confidence=high\nACGTACGTAC\n")
    _write(discover / "candidate_reads.tsv", "read_id\tcandidate_id\nread1\tc1\n")
    quantify = root / "run" / "quantify"
    _write(quantify / "copy_number.tsv", "family_id\testimated_bp\tconfidence\twarning\nF1\t100\thigh\t\n")
    _write(quantify / "run_config.yaml", "command: quantify\n")
    _write(root / "run" / "automatic_defaults.json", json.dumps({"genome_size_bp": 1000, "genome_size_source": "explicit"}))
    arrays = root / "arrays.bed"
    _write(arrays, "chr1\t0\t30\tF1\t100\t+\tmedium\t\n")
    return discover, quantify / "copy_number.tsv", arrays


def test_render_existing_report_copies_sources_and_recomputes_only_comparison(tmp_path: Path) -> None:
    discover, copy_number, arrays = _inputs(tmp_path / "source")
    out = render_existing_family_report(discover=discover, copy_number=copy_number, arrays=arrays,
                                        recovery=None, outdir=tmp_path / "rendered")
    assert (out / "report.html").is_file()
    assert (out / "compare/assembly_vs_read_cn.tsv").is_file()
    manifest = json.loads((out / "source_reuse_manifest.json").read_text())
    assert manifest["raw_reads_accessed"] is False
    assert manifest["compare_recomputed_from_frozen_tables"] is True
    assert all(row["source_sha256"] == row["destination_sha256"] for row in manifest["files"])
    assert json.loads((out / "comparison_reuse_config.json").read_text())["collapse_threshold"] == 0.6
    assert (out / "automatic_defaults.json").is_file()


def test_render_existing_report_rejects_existing_outdir_and_changed_recovery_lock(tmp_path: Path) -> None:
    discover, copy_number, arrays = _inputs(tmp_path / "source")
    existing = tmp_path / "existing"; existing.mkdir()
    with pytest.raises(ValueError, match="must not already exist"):
        render_existing_family_report(discover=discover, copy_number=copy_number, arrays=arrays,
                                      recovery=None, outdir=existing)
    recovery = tmp_path / "recovery"
    outputs = {}
    for name in ("recovery_candidates.tsv", "recovery_validation.tsv", "recovery_loci.bed", "recruited_reads.tsv", "recovered_sequences.fasta"):
        _write(recovery / name, "original\n")
        outputs[name] = _digest(recovery / name)
    _write(recovery / "candidate_lock.json", json.dumps({"status": "candidate_generation_complete", "outputs": outputs}))
    _write(recovery / "recovery_candidates.tsv", "changed\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        render_existing_family_report(discover=discover, copy_number=copy_number, arrays=arrays,
                                      recovery=recovery, outdir=tmp_path / "changed")


def test_render_existing_report_copies_optional_final_recovery_assessment(tmp_path: Path) -> None:
    discover, copy_number, arrays = _inputs(tmp_path / "source")
    recovery = tmp_path / "recovery"
    outputs = {}
    for name in ("recovery_candidates.tsv", "recovery_validation.tsv", "recovery_loci.bed", "recruited_reads.tsv", "recovered_sequences.fasta"):
        _write(recovery / name, "original\n")
        outputs[name] = _digest(recovery / name)
    lock = {"status": "candidate_generation_complete", "outputs": outputs}
    _write(recovery / "candidate_lock.json", json.dumps(lock))
    _write(recovery / "recovery_assessment.json", json.dumps({"status": "assessment_complete", "validated_repeat_gain": False, "decision": "stop_recovery_expansion_negative_poc", "summary": "No validated repeat gain."}))
    _write(recovery / "recovery_validation_report.html", "<h1>Validation evidence</h1>")

    out = render_existing_family_report(discover=discover, copy_number=copy_number, arrays=arrays,
                                        recovery=recovery, outdir=tmp_path / "rendered")
    assert (out / "recovery/recovery_assessment.json").is_file()
    assert (out / "recovery/recovery_validation_report.html").is_file()
    assert json.loads((out / "recovery/candidate_lock.json").read_text()) == lock
    manifest = json.loads((out / "source_reuse_manifest.json").read_text())
    assert any(row["destination_path"].endswith("recovery/recovery_assessment.json") for row in manifest["files"])
