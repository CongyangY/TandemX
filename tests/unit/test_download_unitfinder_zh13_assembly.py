import gzip
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.download_unitfinder_zh13_assembly import download, fasta_qc


def test_fasta_qc_streams_records_and_composition(tmp_path: Path) -> None:
    path = tmp_path / "test.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">chr1\nACGTNN\n>chr2\nTTAA\n")
    assert fasta_qc(path) == {
        "record_count": 2,
        "base_count": 10,
        "acgt_bases": 8,
        "other_bases": 2,
        "other_fraction": 0.2,
    }


def test_download_verifies_frozen_size_md5_and_fasta(tmp_path: Path) -> None:
    source = tmp_path / "source.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">chr1\nACGTACGT\n")
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "experiment_id": "test",
                "assembly": {
                    "accession": "GWH-TEST",
                    "filename": "source.fa.gz",
                    "url": source.as_uri(),
                    "expected_bytes": source.stat().st_size,
                    "expected_md5": hashlib.md5(source.read_bytes()).hexdigest(),
                },
                "boundary": "input identity only",
            }
        )
    )
    receipt = download(config, tmp_path / "result")
    assert receipt["complete"] is True
    assert receipt["fate"] == "source_enrollment_passed"
    assert receipt["fasta_qc"]["base_count"] == 8
    assert not (tmp_path / "result/source.fa.gz.partial").exists()


def test_download_retains_checksum_failure(tmp_path: Path) -> None:
    source = tmp_path / "source.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">chr1\nACGT\n")
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "experiment_id": "test",
                "assembly": {
                    "accession": "GWH-TEST",
                    "filename": "source.fa.gz",
                    "url": source.as_uri(),
                    "expected_bytes": source.stat().st_size,
                    "expected_md5": "0" * 32,
                },
                "boundary": "input identity only",
            }
        )
    )
    receipt = download(config, tmp_path / "failed")
    assert receipt["complete"] is False
    assert receipt["fate"] == "source_enrollment_failure"
    assert receipt["partial_file"]["bytes"] == source.stat().st_size
