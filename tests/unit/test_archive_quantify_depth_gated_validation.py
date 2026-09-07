from __future__ import annotations

import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_quantify_depth_gated_validation import recompute_decision


def test_recompute_depth_gated_decision_retains_branches_and_pairing() -> None:
    metrics = []
    paired = []
    executions = []
    for seed in (6401, 6402, 6403):
        for condition_index in range(1, 10):
            condition_id = f"condition_{condition_index:03d}"
            coverage = (1, 5, 20)[(condition_index - 1) // 3]
            selected_method = (
                "baseline_total_bases" if coverage == 1 else "empirical_controls"
            )
            executions.extend(
                {
                    "seed": str(seed),
                    "condition_id": condition_id,
                    "coverage": str(coverage),
                    "method": method,
                    "runtime_seconds": str(coverage / 5 + method_index / 10),
                    "peak_rss_mib": str(50 + coverage + method_index),
                }
                for method_index, method in enumerate(
                    ("baseline_total_bases", "empirical_controls")
                )
            )
            for family_index in range(55):
                family_id = f"f{family_index + 1:03d}"
                baseline_error = 0.4
                control_error = 0.5 if coverage == 1 else 0.2
                candidate_error = (
                    baseline_error
                    if selected_method == "baseline_total_bases"
                    else control_error
                )
                for method, error in (
                    ("baseline_total_bases", baseline_error),
                    ("empirical_controls", control_error),
                    ("depth_gated_controls", candidate_error),
                ):
                    metrics.append(
                        {
                            "seed": str(seed),
                            "condition_id": condition_id,
                            "coverage": str(coverage),
                            "family_id": family_id,
                            "method": method,
                            "absolute_relative_error": str(error),
                            "signed_relative_error": str(-error),
                        }
                    )
                paired.append(
                    {
                        "seed": str(seed),
                        "condition_id": condition_id,
                        "coverage": str(coverage),
                        "family_id": family_id,
                        "candidate_selected_method": selected_method,
                        "outcome": "equal" if coverage == 1 else "improved",
                    }
                )
    decision, resources = recompute_decision(metrics, paired, executions)
    assert decision["status"] == "passed_frozen_simulation_gates"
    assert decision["paired_outcomes"] == {"improved": 990, "equal": 495, "worse": 0}
    assert decision["selected_condition_counts"] == {
        "baseline_total_bases": 9,
        "empirical_controls": 18,
    }
    assert decision["aggregate_mare_reduction"] == 0.4 - (0.4 + 0.2 + 0.2) / 3
    assert len(resources) == 6


def test_committed_depth_gated_archive_manifest_is_hash_complete() -> None:
    root = Path(__file__).resolve().parents[2] / "paper/evidence/quantify_depth_gated_validation_v1"
    manifest = json.loads((root / "archive_manifest.json").read_text())
    assert len(manifest) == 67
    assert all(digest_file(root / row["file"]) == row["sha256"] for row in manifest)
