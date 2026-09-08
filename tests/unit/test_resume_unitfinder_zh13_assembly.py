import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts.resume_unitfinder_zh13_assembly import validate_parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_validate_parent_binds_failure_diagnosis_and_partial(tmp_path: Path) -> None:
    parent = tmp_path / "parent"
    parent.mkdir()
    partial = parent / "genome.fa.gz.partial"
    partial.write_bytes(b"partial")
    receipt = parent / "run_receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "complete": False,
                "fate": "source_enrollment_failure",
                "partial_file": {"file": partial.name},
            }
        )
    )
    diagnosis = parent / "transport_diagnosis.json"
    diagnosis.write_text(
        json.dumps(
            {
                "classification": "truncated_transfer_before_expected_content_length",
                "safe_range_resume_offset": len(b"partial"),
            }
        )
    )
    config = {
        "parent_failure": {
            "directory": str(parent),
            "fate": "source_enrollment_failure",
            "classification": "truncated_transfer_before_expected_content_length",
            "resume_offset": len(b"partial"),
            "artifacts": {
                receipt.name: sha(receipt),
                diagnosis.name: sha(diagnosis),
                partial.name: sha(partial),
            },
        }
    }
    resolved, resolved_partial = validate_parent(config)
    assert resolved == parent
    assert resolved_partial == partial
    partial.write_bytes(b"changed")
    with pytest.raises(ValueError, match="artifact changed"):
        validate_parent(config)
