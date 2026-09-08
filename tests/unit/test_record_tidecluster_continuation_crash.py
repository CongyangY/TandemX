import json

import pytest

from benchmarks.scripts.record_tidecluster_continuation_crash import record_failure


def test_record_failure_preserves_completed_external_stage(tmp_path) -> None:
    (tmp_path / "environment.json").write_text(
        json.dumps({"complete": False, "experiment_id": "continuation"})
    )
    run_dir = tmp_path / "seed6401/matched_period_sensitivity"
    run_dir.mkdir(parents=True)
    for name in ("tc_cmd_args.json", "tc_chunks.bed", "tc_tidehunter.gff3"):
        (run_dir / name).write_text(name)
    (run_dir / "tidehunter.gnu_time.txt").write_text(
        "User time (seconds): 1.0\n"
        "System time (seconds): 0.2\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00\n"
        "Maximum resident set size (kbytes): 3250036\n"
        "Exit status: 0\n"
    )
    profile = tmp_path / "profile"
    profile.mkdir()
    stage = "s6401_matched_period_sensitivity_tidehunter"
    (profile / f"{stage}.stdout.log").write_text("")
    (profile / f"{stage}.stderr.log").write_text("")
    payload = record_failure(tmp_path, stage)
    assert payload["complete"] is False
    assert payload["external_stage_observed_status"] == (
        "completed_after_profiler_exception"
    )
    assert payload["internal_gnu_time"]["maximum_rss_kb"] == 3250036
    with pytest.raises(FileExistsError):
        record_failure(tmp_path, stage)
