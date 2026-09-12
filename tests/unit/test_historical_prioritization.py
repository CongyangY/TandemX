from benchmarks.scripts.evaluate_historical_prioritization import evaluate


def test_endpoint_keeps_na_separate_from_negative():
    ref=[{"family_id":"a","eligibility":"eligible","reference_state":"reference_collapse","old_assembly_bp":"20","read_estimated_bp":"100","observed_gain_bp":"80"},{"family_id":"b","eligibility":"eligible","reference_state":"reference_retained","old_assembly_bp":"90","read_estimated_bp":"100","observed_gain_bp":"10"}]
    result, summary=evaluate(ref, {"a":{"state":"ok","read_estimated_bp":"100"},"b":{"state":"no_unique_correspondence"}}, "srf")
    assert [x["outcome"] for x in result] == ["TP", "N/A"]
    assert summary["available_families"] == 1 and summary["unavailable_families"] == 1


def test_no_estimates_are_na_not_zero_counts_or_ranking():
    ref=[{"family_id":"a","eligibility":"eligible","reference_state":"reference_collapse","old_assembly_bp":"20","read_estimated_bp":"100","observed_gain_bp":"80"}]
    _, summary=evaluate(ref, None, "srf_k151")
    assert summary["reference_proxy_confusion_matrix"] == "N/A"
    assert summary["reference_proxy_recall"] == summary["reference_proxy_precision"] == "N/A"
    assert summary["ranking_concordance_state"] == "N/A_no_available_estimates"


def test_unresolved_reference_proxy_and_nonpositive_estimate_remain_na_not_negative():
    reference = [
        {"family_id":"a","eligibility":"eligible","reference_state":"unresolved","old_assembly_bp":"20","read_estimated_bp":"100","observed_gain_bp":"80"},
        {"family_id":"b","eligibility":"eligible","reference_state":"reference_retained","old_assembly_bp":"20","read_estimated_bp":"100","observed_gain_bp":"80"},
    ]
    rows, summary = evaluate(reference, {"b": {"state":"ok", "read_estimated_bp":"0"}}, "srf_k151")
    assert [row["outcome"] for row in rows] == ["N/A", "N/A"]
    assert summary["available_families"] == 0
