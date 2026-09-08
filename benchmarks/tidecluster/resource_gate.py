"""Conservative preflight gate for a larger TideCluster execution."""
from __future__ import annotations

import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file


def evaluate_resource_gate(
    config_path: Path,
    runtime_memory_bytes: int,
    project_commit: str,
) -> dict[str, object]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if runtime_memory_bytes <= 0:
        raise ValueError("runtime_memory_bytes must be positive")
    rule = config["decision_rule"]
    threshold = float(rule["maximum_successful_stage_fraction_of_runtime_memory"])
    if not 0 < threshold < 1:
        raise ValueError("resource threshold must be in (0, 1)")

    run_dir = Path(config["completed_reference_run"])
    result_path = run_dir / "result.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("complete") is not True:
        raise ValueError("Reference TideCluster run is not complete")
    internal = result.get("internal_gnu_time", {})
    if not internal:
        raise ValueError("Reference run lacks internal GNU-time measurements")
    stage_rss_kb = {
        stage: int(values["maximum_rss_kb"])
        for stage, values in internal.items()
    }
    if any(value <= 0 for value in stage_rss_kb.values()):
        raise ValueError("Reference run has non-positive stage RSS")

    sampling_receipt_path = Path(config["sampling_receipt"])
    sampling = json.loads(sampling_receipt_path.read_text(encoding="utf-8"))
    target_id = config["target_sample_id"]
    samples = sampling.get("samples")
    if not isinstance(samples, dict):
        raise ValueError("Sampling receipt samples must be keyed by sample_id")
    if target_id not in samples:
        raise ValueError(f"Target sample is absent from sampling receipt: {target_id}")
    target = samples[target_id]
    target_path = Path(config["target_sample_path"])
    if target_path.stat().st_size != int(target["bytes"]):
        raise ValueError("Target sample byte size differs from sampling receipt")
    observed_target_sha256 = digest_file(target_path)
    if observed_target_sha256 != target["sha256"]:
        raise ValueError("Target sample SHA-256 differs from sampling receipt")

    peak_stage = max(stage_rss_kb, key=stage_rss_kb.__getitem__)
    peak_rss_kb = stage_rss_kb[peak_stage]
    peak_rss_bytes = peak_rss_kb * 1024
    fraction = peak_rss_bytes / runtime_memory_bytes
    execution_started = False
    if fraction >= threshold:
        fate = "resource_infeasible_preflight"
        decision = rule["if_at_or_above_threshold"]
    else:
        fate = "resource_gate_passed_not_executed"
        decision = rule["if_below_threshold"]

    return {
        "schema_version": 1,
        "experiment_id": config["experiment_id"],
        "complete": True,
        "fate": fate,
        "decision": decision,
        "execution_started": execution_started,
        "reference_run": {
            "path": str(run_dir),
            "result_sha256": digest_file(result_path),
            "stage_maximum_rss_kb": stage_rss_kb,
            "peak_stage": peak_stage,
            "peak_rss_kb": peak_rss_kb,
            "peak_rss_bytes": peak_rss_bytes,
        },
        "runtime": {
            "memory_limit_bytes": runtime_memory_bytes,
            "peak_reference_fraction": fraction,
            "minimum_free_fraction_required": 1 - threshold,
        },
        "target": {
            "sample_id": target_id,
            "path": str(target_path),
            "bytes": target_path.stat().st_size,
            "total_bases": int(target["total_bases"]),
            "sha256": observed_target_sha256,
        },
        "config_sha256": digest_file(config_path),
        "sampling_receipt_sha256": digest_file(sampling_receipt_path),
        "project_commit": project_commit,
        "scope": config["scope"],
        "boundary": config["boundary"],
    }


def write_gate_receipt(path: Path, receipt: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=False)
    with path.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
