"""Regression for the evidence-bound M1 backend promotion gate."""

from benchmarks.m1_shared_signature.check_backend_eligibility import audit


def test_existing_mapping_backend_fails_public_read_and_endpoint_gate() -> None:
    result = audit()
    checks = result["checks"]
    assert checks["frozen_acd_source_matches_receipt"]
    assert checks["frozen_occupancy_source_matches_receipt"]
    assert checks["same_12_development_inputs"]
    assert checks["mapping_refuses_identical_family_attribution"]
    assert not checks["mapping_accepts_full_reads"]
    assert not checks["independent_full_read_copy_bp_validation_in_m1_receipts"]
    assert not checks["paired_public_backend_runtime_rss_in_m1_receipts"]
    assert result["status"] == "not_eligible"
