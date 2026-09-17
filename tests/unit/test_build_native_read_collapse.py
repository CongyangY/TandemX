"""Exact edit and fail-closed checks for native-read controlled collapse."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.scripts.build_native_read_collapse import generate, sha256_file, sha256_sequence


def fixture(tmp_path: Path) -> Path:
    source = tmp_path / "context.fa"
    sequence = "A" * 500 + "ACGT" * 4 + "T" * 500
    source.write_text(">C1_context\n" + sequence + "\n")
    reads = tmp_path / "native.fastq"
    reads.write_text("@native_1\n" + sequence + "\n+\n" + "I" * len(sequence) + "\n")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "schema_version": 1, "split": "development", "material_id": "toy",
        "raw_read_bundle": str(reads), "raw_read_bundle_sha256": sha256_file(reads),
        "read_pairing_status": "lineage_context_supported_donor_unverified",
        "independent_donor_count": 1,
        "contexts": [{"array_id": "C1", "reference_fasta": str(source),
                      "reference_fasta_sha256": sha256_file(source),
                      "reference_record_id": "C1_context", "source_genome_interval": [1000, 2016],
                      "array_start0": 500, "array_end0": 516, "operational_unit_bp": 4,
                      "array_sequence_sha256": sha256_sequence(sequence[500:516]),
                      "supporting_read_ids": ["native_1"]}]
    }))
    return config


def test_native_context_edits_preserve_reads_and_exact_delta(tmp_path: Path) -> None:
    config = fixture(tmp_path)
    out = tmp_path / "edited"
    receipt = generate(config, out)
    assert receipt["raw_reads_modified"] is False
    assert receipt["truth_scope"] == "injected_edit_delta_only"
    assert len(receipt["cases"]) == 9
    assert sha256_file(Path(receipt["raw_read_bundle"])) == receipt["raw_read_bundle_sha256"]
    for row in receipt["cases"]:
        edited = "".join((out / f"{row['case_id']}.fa").read_text().splitlines()[1:])
        assert len(edited) == 1016 - row["injected_deleted_bp"]
        assert row["original_array_bp"] - row["edited_array_bp"] == row["injected_deleted_bp"]
        assert row["source_deleted_interval"][1] - row["source_deleted_interval"][0] == row["injected_deleted_bp"]
        assert sha256_file(out / f"{row['case_id']}.fa") == row["edited_fasta_sha256"]
    assert next(row for row in receipt["cases"] if row["case_id"] == "C1_terminal_100")["edited_array_bp"] == 0


def test_rejects_changed_read_bundle_and_unsupported_pairing(tmp_path: Path) -> None:
    config = fixture(tmp_path)
    payload = json.loads(config.read_text())
    Path(payload["raw_read_bundle"]).write_text("@native_1\nACGT\n+\nIII\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        generate(config, tmp_path / "bad")
    fixture(tmp_path)
    payload = json.loads(config.read_text())
    payload["read_pairing_status"] = "verified_same_donor_haplotype"
    config.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="requires evidence"):
        generate(config, tmp_path / "also_bad")


def test_selects_declared_record_from_shared_context_fasta(tmp_path: Path) -> None:
    config = fixture(tmp_path)
    payload = json.loads(config.read_text())
    source = Path(payload["contexts"][0]["reference_fasta"])
    source.write_text(">other\n" + "C" * 1016 + "\n" + source.read_text())
    payload["contexts"][0]["reference_fasta_sha256"] = sha256_file(source)
    config.write_text(json.dumps(payload))
    receipt = generate(config, tmp_path / "selected")
    assert receipt["cases"][0]["source_record_id"] == "C1_context"
