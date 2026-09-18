import json

import pytest

from benchmarks.scripts.build_native_bp_collapse import generate
from benchmarks.scripts.verify_native_bp_collapse import verify


def test_independent_bp_edit_verifier_rejects_corrupt_sequence(tmp_path):
    context = tmp_path / "source.fa"
    reads = tmp_path / "reads.fastq.gz"
    manifest = tmp_path / "manifest.json"
    source = "C" * 1200 + "ACGT" * 75 + "G" * 1200
    context.write_text(">E1|Chr1:100-2800|array:1300-1600|family:F\n" + source + "\n")
    reads.write_bytes(b"unchanged-original-read-placeholder")
    import hashlib
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest.write_text(json.dumps({
        "split": "development", "status": "development_source_enrolled_native_interval_supported",
        "sample": "unit", "source_pairing": "test",
        "artifact_sha256": {context.name: digest(context), reads.name: digest(reads)},
        "array": {"chromosome": "Chr1", "context_start0": 100, "context_end0": 2800,
                  "start0": 1300, "end0": 1600, "context_array_start0": 1200,
                  "context_array_end0": 1500, "length_bp": 300, "family_id": "F"}}))
    output = tmp_path / "output"
    generate(manifest, context, reads, output)
    result = verify(manifest, context, reads, output)
    assert result["case_count"] == 9
    fasta = output / "E1_internal_050.fa"
    fasta.write_text(fasta.read_text().replace("G", "A", 1))
    with pytest.raises(ValueError, match="hash mismatch"):
        verify(manifest, context, reads, output)


def test_provisional_record_support_cannot_become_molecule_truth(tmp_path):
    import hashlib

    context = tmp_path / "source.fa"
    reads = tmp_path / "reads.fastq.gz"
    sequence = "C" * 1200 + "ACGT" * 75 + "G" * 1200
    context.write_text(">ctg:100-2800 selected\n" + sequence + "\n")
    reads.write_bytes(b"unaltered-source-placeholder")
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    data = {
        "split": "development", "status": "development_provisional_record_support",
        "support_tier": "full_reference_records_zmw_unverified", "sample": "toy",
        "source_pairing": "exact_extraction_unverified", "context_record_id": "ctg:100-2800",
        "case_prefix": "MJ1",
        "artifact_sha256": {context.name: digest(context), reads.name: digest(reads)},
        "array": {"chromosome": "ctg", "context_start0": 100, "context_end0": 2800,
                  "start0": 1300, "end0": 1600, "context_array_start0": 1200,
                  "context_array_end0": 1500, "length_bp": 300, "family_id": "F"},
    }
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(data))
    output = tmp_path / "generated"
    receipt = generate(manifest, context, reads, output)
    assert receipt["source_support_tier"] == "full_reference_records_zmw_unverified"
    assert receipt["independent_molecule_count"] is None
    assert verify(manifest, context, reads, output)["case_count"] == 9
    altered = json.loads((output / "receipt.json").read_text())
    altered["independent_molecule_count"] = 7
    (output / "receipt.json").write_text(json.dumps(altered))
    with pytest.raises(ValueError, match="promoted"):
        verify(manifest, context, reads, output)
