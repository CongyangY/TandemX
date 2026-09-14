from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from benchmarks.challenge.unified_correspondence import match_native_catalogue
from benchmarks.challenge.unified_native_io import load_tandemx_catalogue
from benchmarks.scripts.run_historical_prioritization_srf_mapping import (
    ACK,
    PRIORITIZATION_FIELDS,
    assign_srf_abundance,
    kmc_command,
    mapping_command,
    normalize_amounts,
    parse_df_pk_free_bytes,
    plan_rows,
    score_estimates,
    srf_commands,
    validate_config,
    write_tsv,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "benchmarks/configs/historical_prioritization_srf_mapping_v1_20260913.json"


def config() -> dict:
    return json.loads(CONFIG.read_text())


def test_frozen_config_and_six_cell_order_validate_without_large_hashes() -> None:
    value = config()
    result = validate_config(value, ROOT, verify_files=False)
    assert result["verified"] is False
    assert [(row["species"], row["method"]) for row in plan_rows(value)] == [
        ("ey15_2", "srf_k151"),
        ("ey15_2", "srf_k101"),
        ("ey15_2", "competitive_mapping"),
        ("macadamia_jansenii", "srf_k151"),
        ("macadamia_jansenii", "srf_k101"),
        ("macadamia_jansenii", "competitive_mapping"),
    ]
    assert value["preregistration_guard"]["required_execution_ack"] == ACK
    catalogue = load_tandemx_catalogue(ROOT / value["species"][0]["catalogue"]["path"])
    assert len(catalogue) == 2133
    assert "TXF000001" in catalogue


def test_macadamia_fastq_order_and_mapping_rule_are_immutable() -> None:
    value = config()
    changed = copy.deepcopy(value)
    changed["species"][1]["reads"].reverse()
    with pytest.raises(ValueError, match="FASTQ order"):
        validate_config(changed, ROOT, verify_files=False)
    changed = copy.deepcopy(value)
    changed["competitive_mapping"]["minimum_identity"] = 0.91
    with pytest.raises(ValueError, match="mapping rule"):
        validate_config(changed, ROOT, verify_files=False)


def test_native_commands_use_fastq_list_both_files_and_fixed_presets(tmp_path: Path) -> None:
    value = config()
    species = value["species"][1]
    command = kmc_command(value, tmp_path / "reads.list", tmp_path, 151)
    assert command[1:8] == ["-fq", "-k151", "-t1", "-m2", "-sm", "-ci20", "-cs100000"]
    assert command[8] == f"@{tmp_path / 'reads.list'}"
    commands = dict((name, command) for name, command, _ in srf_commands(value, species, tmp_path, 101))
    assert commands["map"][-2:] == [row["path"] for row in species["reads"]]
    ordinary = mapping_command(value, species, tmp_path)
    assert ordinary[1:9] == ["-x", "map-hifi", "-c", "-N1000000", "-f1000", "-r100,100", "-t1", str(tmp_path / "templates.fa")]
    assert ordinary[-2:] == [row["path"] for row in species["reads"]]


def test_unique_correspondence_assigns_only_unique_native_abundance() -> None:
    tandemx = {"TXF1": "ACGT", "TXF2": "AAAA"}
    native = {"srf_unique": "ACGTACGT", "srf_unmatched": "CCCC"}
    correspondence = match_native_catalogue(native, tandemx, 0.9, 64)
    assigned, unassigned = assign_srf_abundance(
        {"srf_unique": 100, "srf_unmatched": 30}, correspondence, tandemx
    )
    assert assigned == {"TXF1": 100, "TXF2": 0}
    assert unassigned == [{"native_id": "srf_unmatched", "status": "unmatched",
                           "native_retained_read_bp": 30, "candidate_matches": ""}]


def test_depth_normalization_and_zero_are_distinct_from_technical_na() -> None:
    assert normalize_amounts({"TXF1": 20, "TXF2": 0}, 100, 1000) == {"TXF1": 200.0, "TXF2": 0.0}
    reference = [
        {"family_id": "TXF1", "old_assembly_bp": "50", "new_assembly_bp": "200",
         "reference_state": "reference_collapse", "eligibility": "eligible"},
        {"family_id": "TXF2", "old_assembly_bp": "100", "new_assembly_bp": "200",
         "reference_state": "reference_collapse", "eligibility": "eligible"},
    ]
    rows, summary = score_estimates(reference, {"TXF1": 200.0, "TXF2": 0.0}, "srf_k151")
    assert rows[0]["outcome"] == "TP"
    assert rows[1]["state"] == "zero_native_support"
    assert rows[1]["read_estimated_bp"] == 0.0
    assert rows[1]["predicted_proxy_positive"] == "N/A"
    assert summary["available_families"] == 1
    failed, failed_summary = score_estimates(reference, None, "srf_k151", "stage_failed:count")
    assert all(row["read_estimated_bp"] == "N/A" for row in failed)
    assert failed_summary["available_families"] == 0


def test_darwin_df_parser_uses_available_kib_column_for_exfat() -> None:
    output = (
        "Filesystem   1024-blocks       Used  Available Capacity Mounted on\n"
        "/dev/disk4s2  3906880640 1659918208 2246962432    43% /Volumes/T7\n"
    )
    assert parse_df_pk_free_bytes(output) == 2_300_889_530_368
    with pytest.raises(ValueError, match="Malformed"):
        parse_df_pk_free_bytes("Filesystem\n/dev/disk4s2 bad\n")


def test_technical_failure_rows_serialize_with_explicit_native_bp_field(tmp_path: Path) -> None:
    reference = [{"family_id": "TXF1", "old_assembly_bp": "50", "new_assembly_bp": "200",
                  "reference_state": "reference_collapse", "eligibility": "eligible"}]
    rows, _ = score_estimates(reference, None, "srf_k151", "CellFailure:stage_failed:assemble")
    assert set(rows[0]) <= set(PRIORITIZATION_FIELDS)
    output = tmp_path / "terminal_failure.tsv"
    write_tsv(output, rows, list(PRIORITIZATION_FIELDS))
    fields = output.read_text().splitlines()
    assert len(fields) == 2
    assert "native_retained_read_bp" in fields[0].split("\t")
    assert "CellFailure:stage_failed:assemble" in fields[1]
