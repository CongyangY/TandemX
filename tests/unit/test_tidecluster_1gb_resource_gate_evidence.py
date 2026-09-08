import hashlib
import json
from pathlib import Path


def test_committed_tidecluster_1gb_resource_gate_is_hash_complete() -> None:
    root = (
        Path(__file__).resolve().parents[2]
        / "paper/evidence/tidecluster_1gb_resource_gate_v1"
    )
    manifest = json.loads((root / "archive_manifest.json").read_text())
    assert manifest["complete"] is True
    assert len(manifest["files"]) == 3
    for row in manifest["files"]:
        path = root / row["file"]
        assert path.stat().st_size == row["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
    receipt = json.loads((root / "resource_gate_receipt.json").read_text())
    assert receipt["complete"] is True
    assert receipt["execution_started"] is False
    assert receipt["fate"] == "resource_infeasible_preflight"
    assert receipt["runtime"]["peak_reference_fraction"] > 0.85
