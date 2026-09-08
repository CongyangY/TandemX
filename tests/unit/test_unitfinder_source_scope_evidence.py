import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "paper" / "evidence" / "unitfinder_source_scope_v1" / "audit.json"
SOURCE_SCAN = (
    ROOT / "paper" / "evidence" / "unitfinder_source_scope_v1" / "source_scan.json"
)


def test_unitfinder_source_scope_keeps_reproduction_distinct_from_truth() -> None:
    audit = json.loads(AUDIT.read_text())
    table = audit["supplementary_workbook"]["table_s1"]

    assert audit["status"] == "source_scope_only_no_execution_or_accuracy_result"
    assert table["chromosome_rows"] == 60
    assert table["length_consistency_failure_rows"] == []
    assert table["evidence_role"] == (
        "published_workflow_output_for_reproducibility_not_independent_truth"
    )
    assert "not_a_read_copy_number_estimator" in audit["comparison_boundary"]


def test_unitfinder_source_scope_records_exact_workbook_identity() -> None:
    audit = json.loads(AUDIT.read_text())
    workbook = audit["supplementary_workbook"]

    assert workbook["size_bytes"] == 102935
    assert workbook["sha256"] == (
        "747e240fe57c0d8e665b925cadea15ffa9eea8a1f39f7e0659aef6b9cdf7bf6b"
    )
    assert workbook["sheet_count"] == 14
    assert audit["failed_acquisition_retained"]["size_bytes"] == 1816
    assert audit["planned_real_inputs"]["zh13_t2t_assembly"]["biosample"] == "SAMC1127443"
    assert audit["planned_real_inputs"]["zh13_chip_seq"]["run_accessions"] == [
        f"CRR638{value}" for value in range(211, 217)
    ]
    assert "unresolved" in audit["planned_real_inputs"]["zh13_chip_seq"]["boundary"]


def test_unitfinder_source_scan_retains_reproducibility_failures() -> None:
    scan = json.loads(SOURCE_SCAN.read_text())
    assert scan["complete"] is True
    assert scan["has_dependency_manifest"] is False
    assert scan["has_automated_test_tree"] is False
    assert scan["python_syntax_failures"] == [
        {
            "error": "'(' was never closed",
            "file": "bin/clustering.1.py",
            "line": 99,
        }
    ]
    assert scan["os_system_call_count"] == 85
    assert scan["scope"].endswith("no_accuracy_execution")
