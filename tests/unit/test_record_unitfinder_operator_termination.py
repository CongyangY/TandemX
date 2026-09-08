import csv
import json
from pathlib import Path

import pytest

from benchmarks.scripts.record_unitfinder_operator_termination import record


def test_records_post_start_operator_termination_without_recasting_timeout(
    tmp_path: Path,
) -> None:
    result = tmp_path / "result"
    (result / "profile").mkdir(parents=True)
    (result / "run").mkdir()
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "acceptance": {
                    "required_nonempty_outputs": ["a.Units.fa", "a.Units.uniq.fa"]
                }
            }
        )
    )
    (result / "run_receipt.json").write_text(
        json.dumps(
            {
                "complete": False,
                "fate": "external_process_failure",
                "accuracy_available": False,
            }
        )
    )
    (result / "profile/receipt.json").write_text(json.dumps({"complete": False}))
    with (result / "profile/stages.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("stage", "exit_code", "wall_seconds"), delimiter="\t"
        )
        writer.writeheader()
        writer.writerow({"stage": "help", "exit_code": 0, "wall_seconds": 1})
        writer.writerow({"stage": "smoke", "exit_code": 137, "wall_seconds": 7267.9})
    (result / "run/a.Units.fa").write_text(">one\nACGT\n")
    recorded = record(result, config, 7200, "docker stop --time 30 test")
    assert recorded["formal_benchmark_timeout"] is False
    assert recorded["smoke_stage_exit_code"] == 137
    assert recorded["required_outputs"]["a.Units.fa"]["bytes"] > 0
    assert recorded["required_outputs"]["a.Units.uniq.fa"]["exists"] is False
    with pytest.raises(FileExistsError):
        record(result, config, 7200, "docker stop --time 30 test")
