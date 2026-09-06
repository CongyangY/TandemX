from benchmarks.scripts.compare_discovery_snapshots import summarize_dataset


def test_changed_outputs_and_failed_runs_cannot_pass_performance_gate():
    first = dict(variant="baseline", exit_code=0, output_digest="old", runtime_seconds=2,
                 peak_rss_mib=3, cpu_user_seconds=1, cpu_system_seconds=0.1)
    second = {**first, "variant": "native_seed", "runtime_seconds": 1, "output_digest": "changed"}
    result = summarize_dataset("example", [first, second])
    assert result["speedup_baseline_over_native"] == 2
    assert not result["all_six_outputs_identical"]
    second["exit_code"] = 1
    result = summarize_dataset("example", [first, second])
    assert result["speedup_baseline_over_native"] == "NA"
    assert not result["all_successful"]
