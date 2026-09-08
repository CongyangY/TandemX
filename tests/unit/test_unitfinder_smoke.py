from pathlib import Path

import pytest

from benchmarks.unitfinder.smoke import generate


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
