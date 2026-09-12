import hashlib
import json

from benchmarks.scripts.audit_real_30x_recovery import _stream_sha256, audit_root


def receipt(*, complete=True, reads=12, bases=1200, sha=True):
    data = {"complete": complete, "read_count": reads, "total_bases": bases}
    if sha:
        data["input_sha256"] = "a" * 64
    return data


def test_zero_byte_block_receipt_with_chunks_is_recoverable_not_valid(tmp_path):
    block = tmp_path / "K30076" / "block_01_attempt_003"
    block.mkdir(parents=True)
    (block / "block_receipt.json").touch()
    (block / "chunk_00001.fastq.gz").write_bytes(b"nonzero recovered chunk")
    report = audit_root(tmp_path, [block])
    item = report["blocks"][0]
    assert item["status"] == "recoverable_invalid"
    assert "receipt_zero_bytes" in item["reasons"]
    assert item["nonzero_payload_file_count"] == 1


def test_complete_hashed_receipt_with_totals_is_valid(tmp_path):
    block = tmp_path / "YSD56" / "partition_00001"
    block.mkdir(parents=True)
    (block / "execution.json").write_text(json.dumps(receipt()))
    report = audit_root(tmp_path, [block])
    assert report["blocks"][0]["status"] == "valid"


def test_invalid_totals_with_payload_remain_recoverable_invalid(tmp_path):
    block = tmp_path / "V14167" / "partition_00001"
    block.mkdir(parents=True)
    (block / "run_receipt.json").write_text(json.dumps(receipt(reads=0)))
    (block / "trf.txt").write_text("partial tool output\n")
    report = audit_root(tmp_path, [block])
    item = report["blocks"][0]
    assert item["status"] == "recoverable_invalid"
    assert "missing_positive_record_total" in item["reasons"]


def test_declared_required_file_is_checked_before_validity(tmp_path):
    block = tmp_path / "YSD56" / "partition_00001"
    block.mkdir(parents=True)
    data = receipt()
    data["required_files"] = ["normalized.tsv"]
    (block / "execution.json").write_text(json.dumps(data))
    (block / "partial.tsv").write_text("not the declared artifact\n")
    report = audit_root(tmp_path, [block])
    item = report["blocks"][0]
    assert item["status"] == "recoverable_invalid"
    assert "declared_required_file_missing_or_zero:normalized.tsv" in item["reasons"]


def test_missing_root_is_explicitly_blocked(tmp_path):
    report = audit_root(tmp_path / "not-mounted")
    assert report["root_status"] == "blocked"
    assert report["blocks"] == []


def test_missing_receipt_without_payload_is_missing(tmp_path):
    block = tmp_path / "S245" / "block_01"
    block.mkdir(parents=True)
    report = audit_root(tmp_path, [block])
    item = report["blocks"][0]
    assert item["status"] == "missing"
    assert item["reasons"] == ["receipt_missing"]


def test_deep_files_are_not_inspected_past_payload_bound(tmp_path):
    block = tmp_path / "K30076" / "block_01"
    deep = block / "one" / "two" / "three"
    deep.mkdir(parents=True)
    (block / "block_receipt.json").touch()
    (deep / "chunk.fastq.gz").write_bytes(b"not inspected at depth one")
    report = audit_root(tmp_path, [block], payload_max_depth=1)
    item = report["blocks"][0]
    assert item["status"] == "corrupt"
    assert item["nonzero_payload_file_count"] == 0


def test_multiple_explicit_blocks_are_audited_without_root_discovery(tmp_path):
    valid = tmp_path / "YSD56" / "partition_00001"
    missing = tmp_path / "S245" / "block_01"
    valid.mkdir(parents=True); missing.mkdir(parents=True)
    (valid / "execution.json").write_text(json.dumps(receipt()))
    report = audit_root(tmp_path, [valid, missing])
    assert [item["status"] for item in report["blocks"]] == ["valid", "missing"]


def test_stream_sha256_does_not_use_read_bytes(monkeypatch, tmp_path):
    path = tmp_path / "largeish.bin"
    path.write_bytes(b"A" * (2 << 20))
    monkeypatch.setattr(type(path), "read_bytes", lambda *_: (_ for _ in ()).throw(AssertionError("read_bytes used")))
    assert _stream_sha256(path, 4096) == hashlib.sha256(b"A" * (2 << 20)).hexdigest()


def test_payload_entry_limit_is_blocked_not_silently_treated_as_empty(tmp_path):
    block = tmp_path / "K30076" / "block_02"
    block.mkdir(parents=True)
    (block / "block_receipt.json").touch()
    for index in range(3):
        (block / f"chunk_{index}.fastq.gz").write_bytes(b"payload")
    item = audit_root(tmp_path, [block], payload_max_entries=1)["blocks"][0]
    assert item["status"] == "blocked"
    assert "payload_entry_limit_exceeded" in item["reasons"]


def test_large_receipt_is_bounded_and_never_loaded_as_json(tmp_path):
    block = tmp_path / "V14167" / "block_02"
    block.mkdir(parents=True)
    (block / "block_receipt.json").write_bytes(b"{" + b"x" * ((10 << 20) + 1))
    item = audit_root(tmp_path, [block])["blocks"][0]
    assert item["status"] == "corrupt"
    assert item["reasons"] == ["receipt_exceeds_bounded_max_bytes"]
