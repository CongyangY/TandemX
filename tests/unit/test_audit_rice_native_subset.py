import gzip
import hashlib

import pytest

from benchmarks.scripts.audit_rice_native_subset import audit


def test_audit_valid_original_records(tmp_path) -> None:
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wb") as stream:
        stream.write(b"@r1\nACGT\n+\nIIII\n@r2\nTTA\n+\nHHH\n")
    result = audit(path, hashlib.sha256(path.read_bytes()).hexdigest())
    assert result["read_count"] == 2
    assert result["total_bases"] == 7
    assert result["unique_read_ids"] == 2


def test_audit_rejects_truncated_fastq(tmp_path) -> None:
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wb") as stream:
        stream.write(b"@r1\nACGT\n+\n")
    with pytest.raises(ValueError, match="sequence/quality"):
        audit(path, hashlib.sha256(path.read_bytes()).hexdigest())
