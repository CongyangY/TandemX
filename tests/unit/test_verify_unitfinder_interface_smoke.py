import csv
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts.verify_unitfinder_interface_smoke import verify


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload) + "\n")


def test_verify_unitfinder_interface_smoke_checks_artifacts_independently(
    tmp_path: Path,
) -> None:
    result = tmp_path / "result"
    (result / "run").mkdir(parents=True)
    (result / "profile").mkdir()
    parent = result / "parent_build_snapshot/run_receipt.json"
    parent.parent.mkdir()
    parent.write_text("parent\n")
    parent_sha = hashlib.sha256(parent.read_bytes()).hexdigest()
    names = ["a.Units.fa", "a.Consensus.fa", "a.Units.uniq.fa", "a.Consensus.uniq.fa"]
    output_records = {}
    for name in names:
        path = result / "run" / name
        path.write_text(">x\nACGT\n")
        output_records[name] = {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    image_id = "sha256:" + "b" * 64
    config = tmp_path / "config.json"
    write_json(
        config,
        {
            "image_id": image_id,
            "parent_build_success": {
                "artifacts": {"run_receipt.json": parent_sha}
            },
            "acceptance": {"required_nonempty_outputs": names},
        },
    )
    write_json(
        result / "run_receipt.json",
        {
            "complete": True,
            "fate": "smoke_passed",
            "accuracy_available": False,
            "missing_required_outputs": [],
            "outputs": output_records,
        },
    )
    write_json(
        result / "environment.json",
        {
            "complete": True,
            "image": {"image_id": image_id},
            "parent_build_success": {"image_id": image_id},
        },
    )
    write_json(
        result / "profile/receipt.json",
        {"complete": True, "requested_stage_count": 2, "completed_stage_count": 2},
    )
    with (result / "profile/stages.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("stage", "exit_code"), delimiter="\t")
        writer.writeheader()
        writer.writerows(
            ({"stage": "help", "exit_code": 0}, {"stage": "smoke", "exit_code": 0})
        )
    verified = verify(config, result)
    assert verified["complete"] is True
    assert verified["accuracy_available"] is False
    assert verified["parent_build_artifact_checks"] == 1
    assert verified["required_fasta"][names[0]]["records"] == 1

    (result / "run" / names[0]).write_text(">x\nTGCA\n")
    with pytest.raises(ValueError, match="required output changed"):
        verify(config, result)
