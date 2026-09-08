import hashlib
import json
from pathlib import Path

import pytest

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
    config = json.loads(
        (ROOT / "benchmarks/configs/unitfinder_interface_smoke_v1.json").read_text()
    )
    dockerfile = ROOT / "benchmarks/containers/unitfinder/Dockerfile"
    assert config["status"] == "frozen_before_container_build_and_unitfinder_output"
    assert config["source_commit"] == (
        "e80bff38a1a2718d0b6cf5fd769e6d484fc2ce28"
    )
    assert config["dockerfile_sha256"] == hashlib.sha256(
        dockerfile.read_bytes()
    ).hexdigest()
    assert config["input"]["sequence_length_bp"] == 308375
    assert config["acceptance"]["no_accuracy_claim"] is True
