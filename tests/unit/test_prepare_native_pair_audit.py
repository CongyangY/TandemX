"""Input-integrity checks for bounded real-read copies."""

import gzip
import hashlib
import json

import pytest

from benchmarks.scripts.prepare_native_pair_audit import extract_selected, prepare


def test_prepare_checks_source_copy_and_complete_fastq(tmp_path):
    source = tmp_path / "source.fastq.gz"
    with gzip.open(source, "wb") as handle:
        handle.write(b"@r1\nACGT\n+\nIIII\n@r2\nTT\n+\nII\n")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"schema_version": 1, "subsets": [{
        "dataset": "demo", "source": str(source), "bytes": source.stat().st_size,
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "read_count": 2, "total_bases": 6}]}))
    receipt = prepare(config, tmp_path / "out", reserve_bytes=0)
    assert receipt["subsets"][0]["gzip_trailer_checked"]
    assert (tmp_path / "out" / "demo.fastq.gz").read_bytes() == source.read_bytes()


def test_prepare_rejects_changed_source(tmp_path):
    source = tmp_path / "source.fastq.gz"
    source.write_bytes(b"not gzip")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"schema_version": 1, "subsets": [{
        "dataset": "demo", "source": str(source), "bytes": source.stat().st_size,
        "sha256": "0" * 64, "read_count": 1, "total_bases": 1}]}))
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        prepare(config, tmp_path / "out", reserve_bytes=0)
    assert not (tmp_path / "out" / "demo.fastq.gz").exists()


def test_extract_selected_preserves_record_and_checks_complete_source(tmp_path):
    source = tmp_path / "source.fastq.gz"
    records = b"@r1 some metadata\nACGT\n+\nIIII\n@r2\nTT\n+\nII\n"
    with gzip.open(source, "wb") as handle:
        handle.write(records)
    target = tmp_path / "picked.fastq.gz"
    receipt = extract_selected(source, hashlib.sha256(source.read_bytes()).hexdigest(),
                               {"r1"}, target)
    assert receipt["source_read_count"] == 2
    assert gzip.open(target, "rb").read() == records.split(b"@r2")[0]
    with pytest.raises(ValueError, match="missing requested reads"):
        extract_selected(source, hashlib.sha256(source.read_bytes()).hexdigest(),
                         {"absent"}, target)
