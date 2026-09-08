from benchmarks.scripts.record_unitfinder_source_failure import classify


def test_classifies_short_invalid_gzip_with_full_remote_length_as_truncation() -> None:
    assert (
        classify(260_472_580, 286_375_748, 1, "286375748")
        == "truncated_transfer_before_expected_content_length"
    )


def test_does_not_guess_when_remote_length_differs() -> None:
    assert classify(260, 286, 1, "260") == "unresolved_source_enrollment_failure"
