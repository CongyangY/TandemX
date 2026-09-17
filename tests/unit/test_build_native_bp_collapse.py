"""Exact base-pair edit checks where monomer copy boundaries are unknown."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.scripts.build_native_bp_collapse import generate, sha256_file


def test_bp_editor_preserves_native_reads_and_exact_intervals(tmp_path: Path) -> None:
    sequence = "A" * 1000 + "ACGT" * 26 + "T" * 1000
    context = tmp_path / "context.fa"
    context.write_text(">E1|Chr1:10000-12104|array:11000-11104|family:F1\n" + sequence + "\n")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@native\nACGT\n+\nIIII\n")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "split": "development", "status": "development_source_enrolled_native_interval_supported",
        "sample": "toy", "source_pairing": "same_sample_extraction_unverified",
        "artifact_sha256": {context.name: sha256_file(context), reads.name: sha256_file(reads)},
        "array": {"chromosome": "Chr1", "context_start0": 10000, "context_end0": 12104,
                  "start0": 11000, "end0": 11104, "context_array_start0": 1000,
                  "context_array_end0": 1104, "length_bp": 104, "family_id": "F1"}
    }))
    out = tmp_path / "edited"
    result = generate(manifest, context, reads, out)
    assert len(result["cases"]) == 9
    assert result["original_native_fastq_sha256"] == sha256_file(reads)
    assert next(row for row in result["cases"] if row["case_id"] == "E1_terminal_100")["injected_deleted_bp"] == 104
    for row in result["cases"]:
        edited = (out / f"{row['case_id']}.fa").read_text().splitlines()[1]
        left, right = row["source_deleted_context_interval"]
        assert edited == sequence[:left] + sequence[right:]
        assert len(edited) == len(sequence) - row["injected_deleted_bp"]
        assert row["monomer_copy_truth"].startswith("unknown")
    context.write_text(context.read_text().replace("ACGT", "AAAA", 1))
    with pytest.raises(ValueError, match="hash mismatch"):
        generate(manifest, context, reads, tmp_path / "bad")
