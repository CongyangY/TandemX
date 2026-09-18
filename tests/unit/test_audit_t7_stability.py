"""The storage gate must reject a truncated gzip even with a matching file hash."""

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts import audit_t7_stability as audit


def test_input_gate_checks_gzip_payload_not_only_compressed_hash(tmp_path, monkeypatch):
    monkeypatch.setattr(audit, "T7_MOUNT", tmp_path)
    monkeypatch.setattr(audit, "mount_identity", lambda: (1, True))
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b"@r1\nACGT\n+\nIIII\n")

    def manifest() -> Path:
        spec = tmp_path / "inputs.json"
        spec.write_text(json.dumps({"inputs": [{
            "name": "reads", "path": str(path), "expected_bytes": path.stat().st_size,
            "expected_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }]}))
        return spec

    receipt = audit.check_inputs(manifest())
    assert receipt["status"] == "pass"
    assert receipt["results"][0]["gzip_complete"] is True

    path.write_bytes(path.read_bytes()[:-4])
    with pytest.raises((EOFError, OSError)):
        audit.check_inputs(manifest())
