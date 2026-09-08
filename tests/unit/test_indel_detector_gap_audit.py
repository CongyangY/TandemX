import csv
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_indel_gap_audit_matches_frozen_validation_rows_and_hashes() -> None:
    root = Path(__file__).resolve().parents[2]
    audit = json.loads(
        (root / "paper/evidence/indel_detector_gap_audit_v1/audit.json").read_text()
    )
    summary_path = root / audit["source"]["summary_path"]
    paired_path = root / audit["source"]["paired_path"]
    assert digest(summary_path) == audit["source"]["summary_sha256"]
    assert digest(paired_path) == audit["source"]["paired_sha256"]
    with summary_path.open(newline="") as handle:
        summary = list(csv.DictReader(handle, delimiter="\t"))
    summary_index = {(row["scenario"], row["tool"]): row for row in summary}
    with paired_path.open(newline="") as handle:
        paired = list(csv.DictReader(handle, delimiter="\t"))
    paired_index = {row["scenario"]: row for row in paired}
    for scenario in audit["indel_scenarios"]:
        name = scenario["scenario"]
        for tool in ("tandemx", "tidehunter"):
            source = summary_index[(name, tool)]
            observed = scenario[tool]
            assert float(source["array_recall"]) == observed["array_recall"]
            assert float(source["array_precision"]) == observed["array_precision"]
            assert float(source["base_union_f1"]) == observed["base_union_f1"]
            assert float(source["matched_boundary_mae_bp"]) == observed["matched_boundary_mae_bp"]
            assert float(source["median_runtime_seconds"]) == observed["median_runtime_seconds"]
            assert float(source["median_peak_rss_mib"]) == observed["median_direct_child_peak_rss_mib"]
        paired_row = paired_index[name]
        assert float(paired_row["median_runtime_seconds_ratio"]) == scenario[
            "tandemx_to_tidehunter_runtime_ratio"
        ]
        assert float(paired_row["median_peak_rss_mib_ratio"]) == scenario[
            "tandemx_to_tidehunter_peak_rss_ratio"
        ]


def test_indel_gap_audit_keeps_scope_boundary() -> None:
    root = Path(__file__).resolve().parents[2]
    audit = json.loads(
        (root / "paper/evidence/indel_detector_gap_audit_v1/audit.json").read_text()
    )
    assert audit["accuracy_assessment"] == "not_clearly_behind_on_frozen_synthetic_indel_scenarios"
    assert audit["largest_observed_gap"] == "elapsed_time_not_read_local_accuracy"
    assert "real_read_local_interval_truth" in audit["unresolved"]
    assert "synthetic_distribution_level_evidence_only" in audit["warning"]
