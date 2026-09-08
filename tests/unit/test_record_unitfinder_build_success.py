import pytest

from benchmarks.scripts.record_unitfinder_build_success import validate_success_log


def test_validate_success_log_requires_exact_image_export() -> None:
    image = "sha256:" + "a" * 64
    validate_success_log(
        "python /opt/unitFinder/bin/unitFinder.py -h\n"
        f"writing image {image}\n"
        "naming to docker.io/tandemx/unitfinder:e80bff38-v2 done\n",
        image,
    )
    with pytest.raises(ValueError, match="inconsistent"):
        validate_success_log("ERROR: failed to solve\n", image)
