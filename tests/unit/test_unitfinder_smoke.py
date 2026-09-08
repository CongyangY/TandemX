import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.run_unitfinder_interface_smoke import (
    validate_parent_build_failure,
    validate_parent_build_success,
)
from benchmarks.unitfinder.smoke import generate


ROOT = Path(__file__).resolve().parents[2]


def test_unitfinder_smoke_is_deterministic_and_has_frozen_outlier(tmp_path: Path) -> None:
    first = generate(tmp_path / "first")
    second = generate(tmp_path / "second")
    assert first == second
    assert first["target"] == {
        "array_id": "target_high_copy",
        "start": 180775,
        "end": 283375,
        "period": 171,
        "copies": 600,
    }
    assert first["sequence_length_bp"] == 308375
    assert first["files"]["unitfinder_smoke.fa"]["sha256"] == (
        "e7a372298273c1595eba39fe38930bd2be2000e71d33243fe71829c2efacf750"
    )
    assert first["files"]["truth_arrays.tsv"]["sha256"] == (
        "9ebbae8d6e17e8a027463f12003bec760f82f43ba4753b9a7c59f32d968dbb44"
    )
    assert first["scope"].endswith("not_accuracy_evidence")
    with pytest.raises(FileExistsError):
        generate(tmp_path / "first")


def test_unitfinder_smoke_config_is_frozen_before_execution() -> None:
    v1 = json.loads(
        (ROOT / "benchmarks/configs/unitfinder_interface_smoke_v1.json").read_text()
    )
    v2 = json.loads(
        (ROOT / "benchmarks/configs/unitfinder_interface_smoke_v2.json").read_text()
    )
    dockerfile = ROOT / "benchmarks/containers/unitfinder/Dockerfile"
    failed_dockerfile = (
        ROOT
        / "paper/evidence/unitfinder_container_build_v1_failure/results/source_snapshot/Dockerfile"
    )
    assert v1["status"] == "frozen_before_container_build_and_unitfinder_output"
    assert v1["dockerfile_sha256"] == hashlib.sha256(
        failed_dockerfile.read_bytes()
    ).hexdigest()
    assert v2["status"] == (
        "frozen_after_v1_container_failure_before_v2_build_or_smoke_output"
    )
    assert v2["source_commit"] == (
        "e80bff38a1a2718d0b6cf5fd769e6d484fc2ce28"
    )
    assert v2["dockerfile_sha256"] == hashlib.sha256(
        dockerfile.read_bytes()
    ).hexdigest()
    assert v2["input"]["sequence_length_bp"] == 308375
    assert v2["acceptance"]["no_accuracy_claim"] is True

    v3 = json.loads(
        (ROOT / "benchmarks/configs/unitfinder_interface_smoke_v3.json").read_text()
    )
    assert v3["status"] == (
        "frozen_after_v2_image_build_before_external_smoke_output"
    )
    assert v3["image_id"] == (
        "sha256:eadbae183cf29f24b987ea38db3aa5131b1b6b87b1074f45058c60b92929d221"
    )
    assert v3["input"] == v2["input"]
    assert v3["command_arguments"] == v2["command_arguments"]
    assert v3["acceptance"] == v2["acceptance"]


def test_unitfinder_v2_rehashes_parent_build_failure(tmp_path: Path) -> None:
    parent = tmp_path / "parent"
    parent.mkdir()
    receipt = parent / "run_receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "complete": False,
                "fate": "container_build_failure",
                "smoke_execution_started": False,
            }
        )
    )
    log = parent / "docker_build.log"
    log.write_text("failed\n")
    config = {
        "allowed_change_from_v1": (
            "run_tracked_file_sha256_from_opt_unitFinder_checkout_only"
        ),
        "parent_build_failure": {
            "directory": str(parent),
            "fate": "container_build_failure",
            "artifacts": {
                "run_receipt.json": digest_file(receipt),
                "docker_build.log": digest_file(log),
            },
        },
    }
    result = validate_parent_build_failure(config)
    assert result is not None
    assert result["directory"] == parent
    log.write_text("changed\n")
    with pytest.raises(ValueError, match="artifact changed"):
        validate_parent_build_failure(config)


def test_unitfinder_v3_rehashes_parent_build_success(tmp_path: Path) -> None:
    parent = tmp_path / "parent"
    parent.mkdir()
    image_id = "sha256:" + "a" * 64
    receipt = parent / "run_receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "complete": True,
                "fate": "container_build_passed",
                "smoke_execution_started": False,
                "image": {"image_id": image_id},
            }
        )
    )
    lock = parent / "conda_explicit.txt"
    lock.write_text("dependency-lock\n")
    config = {
        "image_id": image_id,
        "parent_build_success": {
            "directory": str(parent),
            "artifacts": {
                "run_receipt.json": digest_file(receipt),
                "conda_explicit.txt": digest_file(lock),
            },
        },
    }
    result = validate_parent_build_success(config)
    assert result is not None
    assert result["receipt"]["image"]["image_id"] == image_id
    lock.write_text("changed\n")
    with pytest.raises(ValueError, match="build-success artifact changed"):
        validate_parent_build_success(config)
