import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_sampling_evidence import archive


def write_source(root: Path) -> None:
    root.mkdir()
    (root / "sampling_plan.json").write_text('{"method":"test"}\n')
    receipt = {
        "complete": True,
        "plan_sha256": digest_file(root / "sampling_plan.json"),
        "samples": [{"sample_id": "sample_001", "read_count": 3, "status": "ok"}],
    }
    (root / "sampling_receipt.json").write_text(json.dumps(receipt) + "\n")
    (root / "sample_001.length_histogram.tsv").write_text(
        "length_bp\tread_count\n100\t2\n200\t1\n"
    )
    (root / "sample_001.joint_distribution.tsv").write_text(
        "length_bin_kb\tgc_bin_percent\tmean_quality_bin_phred\tread_count\n0\t50\t30\t3\n"
    )
    (root / "sample_001.fastq.gz").write_bytes(b"large data must not be archived")


def test_archive_sampling_evidence_copies_only_validated_compact_files(tmp_path: Path) -> None:
    source = tmp_path / "source"
    write_source(source)
    outdir = tmp_path / "archive"
    manifest = archive(source, outdir)
    assert len(manifest) == 4
    assert not (outdir / "sample_001.fastq.gz").exists()
    assert json.loads((outdir / "archive_manifest.json").read_text()) == manifest
    with pytest.raises(ValueError, match="already exists"):
        archive(source, outdir)


def test_archive_sampling_evidence_rejects_distribution_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "source"
    write_source(source)
    (source / "sample_001.length_histogram.tsv").write_text(
        "length_bp\tread_count\n100\t2\n"
    )
    with pytest.raises(ValueError, match="differs"):
        archive(source, tmp_path / "archive")


def test_archive_sampling_evidence_accepts_completed_empty_sample(tmp_path: Path) -> None:
    source = tmp_path / "source"
    write_source(source)
    receipt_path = source / "sampling_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["samples"][0].update(read_count=0, status="empty_random_sample")
    receipt_path.write_text(json.dumps(receipt) + "\n")
    (source / "sample_001.length_histogram.tsv").write_text("length_bp\tread_count\n")
    (source / "sample_001.joint_distribution.tsv").write_text(
        "length_bin_kb\tgc_bin_percent\tmean_quality_bin_phred\tread_count\n"
    )

    manifest = archive(source, tmp_path / "archive")
    assert len(manifest) == 4
