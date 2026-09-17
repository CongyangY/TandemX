"""Keep the expert-facing native-source registry aligned with frozen receipts."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONTROLLED = ROOT / "benchmarks/controlled_collapse"


def _read(path: Path) -> dict:
    return json.loads(path.read_text())


def test_native_source_registry_agrees_with_source_receipts() -> None:
    registry = _read(ROOT / "benchmarks/inputs/native_source_registry_20260917.json")
    rows = {row["dataset"]: row for row in registry["datasets"]}
    assert len(rows) == len(registry["datasets"])
    assert registry["final_native_heldout_enrolled"] is False
    assert registry["independent_physical_copy_truth_enrolled"] is False
    for row in rows.values():
        assert row["split"].startswith("development")
        assert row["final_heldout"] is False
        assert row["biological_accuracy"] == "not_evaluated"
        assert (ROOT / row["source_package"]).is_dir()

    col = rows["Col_CEN_v1.2_ERR6210723"]
    col_config = _read(CONTROLLED / "native_edit_development_v1/source_config.json")
    col_source = _read(CONTROLLED / "native_pair_audit_20260917/eligibility_manifest.json")["datasets"][0]
    col_map = _read(CONTROLLED / "native_genomewide_audit_20260917/score/summary.json")
    assert col["assembly_source_sha256"] == col_config["source_assembly_sha256"] == col_source["assembly_full_sha256"]
    assert col["selected_original_reads_sha256"] == col_config["raw_read_bundle_sha256"]
    assert col["array_coordinates_0based_half_open"] == col_source["controlled_array_intervals"]
    assert col["frozen_case_count"] == 27 and col_map["full_reference_primary_spanner_count"] == 6

    ey = rows["Ey15_2_9994_ERR8666125"]
    ey_source = _read(CONTROLLED / "ey15_native_pair_audit_20260917/source_eligibility_manifest.json")
    ey_edits = _read(CONTROLLED / "ey15_native_bp_development_v1/receipt.json")
    ey_score = _read(CONTROLLED / "ey15_native_span_evidence_20260917/summary.json")
    assert ey["assembly_source_sha256"] == ey_source["assembly"]["sha256_current_full_readback"]
    assert ey["selected_original_reads_sha256"] == ey_source["reads"]["selected_fastq_sha256"]
    assert ey["selected_original_reads_sha256"] == ey_edits["original_native_fastq_sha256"]
    assert ey["array_coordinates_0based_half_open"]["E1"] == "Chr1:12829234-12832478"
    assert ey["frozen_case_count"] == len(ey_edits["cases"]) == ey_score["case_count"] == 9
    assert ey_score["frozen_m2_evaluated_cases"] == 0

    mo = rows["Mo17_SRR15447419_GCA_022117705.1"]
    mo_source = _read(CONTROLLED / "mo17_native_pair_audit_20260917/source_eligibility_manifest_v1.json")
    assert mo["assembly_source_sha256"] == mo_source["assembly"]["sha256"]
    assert mo["selected_original_reads_sha256"] == mo_source["read_screen"]["raw_diagnostic_fastq_sha256"]
    assert mo["screened_assembly_candidate_count"] == len(mo_source["pre_read_candidate_selection"]["region_ids_screened"])
    assert mo["geometric_spanner_count"] == mo_source["read_screen"]["geometric_spanner_count"] == 6
    assert mo["qualified_spanner_count"] == mo_source["read_screen"]["qualifying_spanner_count"] == 0
    assert mo["controlled_edit_status"] == "not_run"

    rice = rows["Nipponbare_AGIS1_SRR25241090"]
    rice_audit = _read(CONTROLLED / "rice_native_pair_audit_20260917/eligibility_audit.json")
    rice_selection = _read(CONTROLLED / "rice_native_pair_audit_20260917/selection_receipt.json")
    rice_subset = _read(CONTROLLED / "rice_native_pair_audit_20260917/subset_integrity_receipt.json")
    assert rice["assembly_source_sha256"] == rice_audit["source_assembly_current_sha256"]
    assert rice["read_subset_sha256"] == rice_subset["sha256"] == rice_audit["read_subset_current_sha256"]
    assert rice["screened_assembly_window_count"] == rice_selection["candidate_windows"] == rice_audit["window_denominator"]
    assert rice["ranked_window_review_count"] == rice_selection["ranked_candidates_checked"] == 1000
    assert rice_selection["selected"] is None and rice["passing_candidate_count"] == 0

    mac = rows["Macadamia_jansenii_SRR13557763_62"]
    mac_audit = _read(CONTROLLED / "macadamia_native_pair_audit_20260917/eligibility_audit.json")
    assert mac["assembly_source_sha256"] is None
    assert mac["assembly_historical_sha256"] == mac_audit["assembly_expected_sha256_historical"]
    assert mac["original_stream_compressed_bytes"] == [row["gzip_size_bytes"] for row in mac_audit["read_runs"]]
    assert mac_audit["current_full_input_hash_and_gzip_readback"] == "not_performed"
