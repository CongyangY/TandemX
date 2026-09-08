import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.tidecluster.resource_gate import evaluate_resource_gate, write_gate_receipt


def _fixture(tmp_path: Path, threshold: float = 0.85) -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "result.json").write_text(
        json.dumps(
            {
                "complete": True,
                "internal_gnu_time": {
                    "tidehunter": {"maximum_rss_kb": 400},
                    "clustering": {"maximum_rss_kb": 800},
                },
            }
        )
    )
    target = tmp_path / "one_gb.fa"
    target.write_text(">seq\nACGT\n")
    sampling = tmp_path / "sampling.json"
    sampling.write_text(
        json.dumps(
            {
                "samples": {
                    "large": {
                        "bytes": target.stat().st_size,
                        "total_bases": 4,
                        "sha256": digest_file(target),
                    }
                }
            }
        )
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "experiment_id": "gate",
                "scope": "same_runtime_only",
                "boundary": "not a cross-host tool limit",
                "completed_reference_run": str(run_dir),
                "sampling_receipt": str(sampling),
                "target_sample_id": "large",
                "target_sample_path": str(target),
                "decision_rule": {
                    "maximum_successful_stage_fraction_of_runtime_memory": threshold,
                    "if_at_or_above_threshold": "refuse",
                    "if_below_threshold": "eligible",
                },
            }
        )
    )
    return config


def test_resource_gate_refuses_without_starting_target(tmp_path: Path) -> None:
    config = _fixture(tmp_path)
    receipt = evaluate_resource_gate(config, runtime_memory_bytes=900 * 1024, project_commit="abc")
    assert receipt["fate"] == "resource_infeasible_preflight"
    assert receipt["execution_started"] is False
    assert receipt["reference_run"]["peak_stage"] == "clustering"


def test_resource_gate_passes_but_does_not_claim_execution(tmp_path: Path) -> None:
    config = _fixture(tmp_path)
    receipt = evaluate_resource_gate(config, runtime_memory_bytes=2000 * 1024, project_commit="abc")
    assert receipt["fate"] == "resource_gate_passed_not_executed"
    assert receipt["decision"] == "eligible"


def test_resource_gate_rechecks_target_hash_and_refuses_overwrite(tmp_path: Path) -> None:
    config = _fixture(tmp_path)
    payload = json.loads(config.read_text())
    Path(payload["target_sample_path"]).write_text("changed")
    with pytest.raises(ValueError, match="byte size"):
        evaluate_resource_gate(config, runtime_memory_bytes=1_000_000, project_commit="abc")
    outdir = tmp_path / "receipt"
    write_gate_receipt(outdir / "resource_gate_receipt.json", {"complete": True})
    with pytest.raises(FileExistsError):
        write_gate_receipt(outdir / "resource_gate_receipt.json", {"complete": True})
