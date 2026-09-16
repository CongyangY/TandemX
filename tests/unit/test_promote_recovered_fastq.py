from __future__ import annotations

import errno
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts.promote_recovered_fastq import promote


def make_source(root: Path) -> tuple[Path, Path]:
    source = root / "spot_11_12.fastq.gz"
    with gzip.open(source, "wb") as handle:
        handle.write(b"@SRR123.11\nACGT\n+\nIIII\n@SRR123.12\nTTAA\n+\nIIII\n")
    receipt = root / "spot_11_12.receipt.json"
    receipt.write_text(json.dumps({
        "status": "complete", "vdb_dump_exit_code": 0,
        "run_accession": "SRR123", "spot_range": {"start": 11, "end": 12},
        "fastq_bytes": source.stat().st_size,
        "fastq_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "record_count": 2, "total_bases": 8,
    }))
    return source, receipt


def test_verified_promotion_keeps_source_and_copies_receipt(tmp_path: Path) -> None:
    source, receipt = make_source(tmp_path)
    target = tmp_path / "new_attempt" / "chunks"
    result = promote(source, receipt, target)
    assert result["status"] == "promoted_and_readback_verified"
    assert source.read_bytes() == (target / source.name).read_bytes()
    assert receipt.read_bytes() == (target / receipt.name).read_bytes()


def test_corrupt_source_refuses_target_creation(tmp_path: Path) -> None:
    source, receipt = make_source(tmp_path)
    source.write_bytes(source.read_bytes()[:-8] + b"corrupt!")
    target = tmp_path / "new_attempt" / "chunks"
    with pytest.raises(ValueError, match="source FASTQ hash"):
        promote(source, receipt, target)
    assert not target.exists()


def test_incomplete_receipt_refuses_target_creation(tmp_path: Path) -> None:
    source, receipt = make_source(tmp_path)
    data = json.loads(receipt.read_text())
    data["status"] = "failed_partial_retained"
    receipt.write_text(json.dumps(data))
    target = tmp_path / "new_attempt" / "chunks"
    with pytest.raises(ValueError, match="not a completed"):
        promote(source, receipt, target)
    assert not target.exists()


def test_target_disk_full_keeps_source_and_removes_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source, receipt = make_source(tmp_path)
    target = tmp_path / "new_attempt" / "chunks"

    def fail_copy(src: object, dst: object, length: int) -> None:
        del length
        dst.write(src.read(4))
        raise OSError(errno.ENOSPC, "injected disk full")

    monkeypatch.setattr("benchmarks.scripts.promote_recovered_fastq.shutil.copyfileobj", fail_copy)
    with pytest.raises(OSError, match="injected disk full"):
        promote(source, receipt, target)
    assert source.is_file() and receipt.is_file()
    assert list(target.iterdir()) == []
