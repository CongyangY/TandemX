import pytest

from benchmarks.scripts.record_unitfinder_build_failure import validate_failure_log


def test_validate_failure_log_classifies_relative_path_bug() -> None:
    payload = validate_failure_log(
        "Transaction finished\n"
        "HEAD is now at e80bff3 Update README.md\n"
        "sha256sum: LICENSE: No such file or directory\n"
        "sha256sum: bin/unitFinder.py: No such file or directory\n"
        "did not complete successfully: exit code: 123\n",
        "e80bff38a1a2718d0b6cf5fd769e6d484fc2ce28",
    )
    assert payload["classification"] == (
        "container_provenance_hash_working_directory_error"
    )
    assert payload["missing_relative_path_count"] == 2
    assert payload["dependency_install_reached"] is True
    assert payload["upstream_interface_started"] is False


def test_validate_failure_log_rejects_unrelated_failure() -> None:
    with pytest.raises(ValueError, match="expected failure evidence"):
        validate_failure_log("network timeout\n", "e80bff38")
