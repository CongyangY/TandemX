import json
from pathlib import Path

import pytest

from benchmarks.scripts.archive_fastq_qc_evidence import archive


def make_qc(root: Path) -> None:
    root.mkdir()
    receipt = {
        "complete": True,
        "fastq_records_valid": True,
        "gzip_trailer_checked": True,
        "exact_duplicate_read_ids": 0,
        "read_count": 3,
        "total_bases": 400,
        "input_sha256": "a" * 64,
    }
    (root / "qc.json").write_text(json.dumps(receipt) + "\n")
    (root / "length_histogram.tsv").write_text("length_bp\tread_count\n100\t2\n200\t1\n")
    (root / "joint_distribution.tsv").write_text(
        "length_bin_kb\tgc_bin_percent\tmean_quality_bin_phred\tread_count\n0\t40\t20\t2\n0\t42\t25\t1\n"
    )
    (root / "base_quality_histogram.tsv").write_text("phred\tbase_count\n20\t300\n25\t100\n")
    (root / "read_ids.sqlite").write_bytes(b"large local index")


def test_archive_fastq_qc_copies_validated_compact_evidence(tmp_path: Path) -> None:
    source = tmp_path / "source"
    make_qc(source)
    outdir = tmp_path / "archive"
    manifest = archive(source, outdir)
    assert len(manifest) == 4
    assert not (outdir / "read_ids.sqlite").exists()
    assert json.loads((outdir / "archive_manifest.json").read_text()) == manifest
    with pytest.raises(ValueError, match="already exists"):
        archive(source, outdir)


def test_archive_fastq_qc_rejects_incomplete_or_mismatched_evidence(tmp_path: Path) -> None:
    source = tmp_path / "source"
    make_qc(source)
    receipt_path = source / "qc.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["exact_duplicate_read_ids"] = 1
    receipt_path.write_text(json.dumps(receipt) + "\n")
    with pytest.raises(ValueError, match="duplicate-free"):
        archive(source, tmp_path / "duplicate")

    receipt["exact_duplicate_read_ids"] = 0
    receipt_path.write_text(json.dumps(receipt) + "\n")
    (source / "base_quality_histogram.tsv").write_text("phred\tbase_count\n20\t399\n")
    with pytest.raises(ValueError, match="Base-quality denominator"):
        archive(source, tmp_path / "mismatch")
