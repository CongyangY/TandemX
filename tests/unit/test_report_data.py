from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from tandemx.io.validators import ValidationError, validate_project
from tandemx.report_data import build_report_data, write_report_data


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _products(root: Path) -> None:
    _write(root / "discover/families.tsv", "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\nTXF1\tTXM1\t171\tmd5\t0.5\t3\t513\t0.99\tfalse\thigh\t\n")
    _write(root / "discover/monomers.fa", ">family_id=TXF1;monomer_id=TXM1;length_bp=171;confidence=high\nACGT\n")
    _write(root / "discover/monomer_membership.tsv", "read_id\tcandidate_id\tcluster_id\tfamily_id\trepresentative_sha256\tedit_distance_upper_bound\tsimilarity_lower_bound\tminimum_cluster_identity\tcompatible_cluster_count\talternative_cluster_ids\tstatus\twarning\nr1\tc1\tC1\tTXF1\thash\t0\t1\t.95\t1\t\tassigned\t\n")
    _write(root / "discover/candidate_reads.tsv", "read_id\tcandidate_id\n" "r1\tc1\n")
    _write(root / "discover/family_hierarchy.tsv", "hierarchy_edge_id\tshorter_family_id\tlonger_family_id\tshorter_length_bp\tlonger_length_bp\tnearest_integer_multiple\tlength_ratio\tmultiple_error\tlocal_identity\tlocal_overlap_fraction_shorter\tshared_kmer_fraction\torientation\tedge_type\tstatus\twarning\nTXH1\tTXF1\tTXF2\t171\t342\t2\t2\t0\t.9\t.9\t.9\tforward\tputative_period_multiple\tcandidate\theuristic_period_multiple_not_validated_hor\n")
    _write(root / "quantify/copy_number.tsv", "family_id\testimated_bp\twarning\nTXF1\t1000\tqwarn\n")
    _write(root / "locate/arrays.bed", "chr1\t0\t100\tTXF1\t10\t+\tmedium\tawarn\n")
    _write(root / "compare/assembly_vs_read_cn.tsv", "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\nTXF1\t1000\t700\t0.7\tpossible_collapse\tmedium\tcwarn\n")


def test_report_data_uses_actual_membership_and_preserves_candidate_architecture(tmp_path: Path) -> None:
    _products(tmp_path)
    data = write_report_data(tmp_path)
    family = data["families"][0]
    assert family["candidate_read_count"] == 1
    assert family["candidate_array_count"] == 1
    assert family["abundance_deficit_bp"] == 300
    assert family["discovery_confidence"] == "high"
    assert family["comparison_confidence"] == "medium"
    assert family["confidence"] == "medium"
    assert family["confidence_scope"] == "discovery+comparison"
    assert family["warning"] == "cwarn;qwarn"
    used_paths = {row["path"] for row in data["sources"] if row["used"]}
    assert {"discover/monomers.fa", "discover/candidate_reads.tsv", "discover/monomer_membership.tsv", "discover/family_hierarchy.tsv"} <= used_paths
    assert "putative_period_multiple" in (tmp_path / "families/repeat_architecture.tsv").read_text()
    network = (tmp_path / "families/family_network.graphml").read_text()
    assert 'source="TXF1"' in network and "validated" not in network
    with (tmp_path / "summary.tsv").open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t"))[0]["estimated_abundance_bp"] == "1000.0"
    assert json.loads((tmp_path / "summary.json").read_text())["schema_version"] == 1


def test_records_exclude_stale_outputs_from_failed_steps(tmp_path: Path) -> None:
    _products(tmp_path)
    records = [
        {"step": "discover", "exit_status": 0, "output_validated": True, "notes": ""},
        {"step": "quantify", "exit_status": 1, "output_validated": False, "notes": "failed"},
        {"step": "locate", "exit_status": 0, "output_validated": True, "notes": ""},
        {"step": "compare", "exit_status": 0, "output_validated": False, "notes": "failed"},
    ]
    family = build_report_data(tmp_path, records)["families"][0]
    assert family["estimated_abundance_bp"] is None
    assert family["assembly_representation_bp"] is None
    assert family["candidate_array_count"] == 1


def test_report_rejects_duplicate_families_and_malformed_arrays(tmp_path: Path) -> None:
    _products(tmp_path)
    families = tmp_path / "discover/families.tsv"
    families.write_text(families.read_text() + families.read_text().splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="Duplicate family_id"):
        build_report_data(tmp_path)
    _products(tmp_path / "arrays")
    _write(tmp_path / "arrays/locate/arrays.bed", "chr1\t0\t10\tTXF1\t100\t+\n")
    with pytest.raises(ValueError, match="expected 8"):
        build_report_data(tmp_path / "arrays")
    _products(tmp_path / "fasta")
    with (tmp_path / "fasta/discover/monomers.fa").open("a") as handle:
        handle.write(">family_id=TXF1;monomer_id=TXM2;length_bp=4;confidence=high\nACGT\n")
    with pytest.raises(ValueError, match="Duplicate family_id"):
        build_report_data(tmp_path / "fasta")


def test_report_preserves_zero_deficit_and_zero_observed_arrays(tmp_path: Path) -> None:
    _products(tmp_path)
    _write(tmp_path / "locate/arrays.bed", "")
    _write(tmp_path / "compare/assembly_vs_read_cn.tsv", "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\nTXF1\t1000\t1200\t1.2\tpossible_overexpansion\tlow\t\n")
    family = build_report_data(tmp_path)["families"][0]
    assert family["candidate_array_count"] == 0
    assert family["abundance_deficit_bp"] == 0
    assert family["confidence"] == "low"


def test_report_rejects_invalid_numeric_and_duplicate_comparison_rows(tmp_path: Path) -> None:
    _products(tmp_path)
    _write(tmp_path / "quantify/copy_number.tsv", "family_id\testimated_bp\tconfidence\twarning\nTXF1\tNaN\thigh\t\n")
    with pytest.raises(ValueError, match="non-finite estimated_bp"):
        build_report_data(tmp_path)
    _products(tmp_path / "duplicate")
    compare = tmp_path / "duplicate/compare/assembly_vs_read_cn.tsv"
    compare.write_text(compare.read_text() + compare.read_text().splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="Duplicate family_id"):
        build_report_data(tmp_path / "duplicate")


def test_report_rejects_canonical_comparison_numeric_conflicts(tmp_path: Path) -> None:
    _products(tmp_path)
    _write(tmp_path / "compare/assembly_vs_read_cn.tsv", "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\nTXF1\t999\t700\t0.7\tpossible_collapse\tmedium\t\n")
    with pytest.raises(ValueError, match="read_estimated_bp conflicts"):
        build_report_data(tmp_path)
    _write(tmp_path / "compare/assembly_vs_read_cn.tsv", "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\nTXF1\t1000\t700\t0.6\tpossible_collapse\tmedium\t\n")
    with pytest.raises(ValueError, match="assembly_read_ratio conflicts"):
        build_report_data(tmp_path)


def test_locate_fallback_preserves_assembly_but_marks_read_mismatch_not_comparable(tmp_path: Path) -> None:
    _products(tmp_path)
    (tmp_path / "compare/assembly_vs_read_cn.tsv").unlink()
    _write(tmp_path / "locate/assembly_vs_read_cn.tsv", "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\nTXF1\t0\t700\t0\tassembly_only\tlow\t\n")
    family = build_report_data(tmp_path)["families"][0]
    assert family["assembly_representation_bp"] == 700
    assert family["assembly_read_ratio"] is None
    assert family["comparison_status"] == "not_comparable"
    assert "missing_validated_comparison" in family["warning"]


def test_report_counts_canonical_statuses_and_metadata_warning(tmp_path: Path) -> None:
    _products(tmp_path)
    _write(tmp_path / "automatic_defaults.json", json.dumps({"genome_size_bp": 1000, "genome_size_source": "assembly_total_length_provisional"}))
    data = build_report_data(tmp_path)
    assert data["summary"]["possible_underrepresented_count"] == 1
    assert data["summary"]["recommended_review_family_ids"] == ["TXF1"]
    assert data["data_availability_warnings"][0]["code"] == "assembly_length_proxy"
    assert data["sources"][0]["sha256"]


def test_graphml_declares_edge_keys_and_retains_isolated_family_node(tmp_path: Path) -> None:
    _products(tmp_path)
    hierarchy = tmp_path / "discover/family_hierarchy.tsv"
    hierarchy.write_text(hierarchy.read_text().splitlines()[0] + "\n")
    write_report_data(tmp_path)
    graph = (tmp_path / "families/family_network.graphml").read_text()
    assert '<key id="edge_type" for="edge"' in graph
    assert '<node id="TXF1"/>' in graph


def test_report_derivatives_have_explicit_schemas_and_recovery_is_status_only(tmp_path: Path) -> None:
    _products(tmp_path)
    _write(tmp_path / "recovery/recovery_candidates.tsv", "locus_id\tfamily_id\tchromosome\tstart\tend\toriginal_assembly_repeat_bp\tread_derived_abundance_bp\tfamily_abundance_deficit_bp\tleft_flank_uniqueness\tright_flank_uniqueness\trecruited_read_count\tflank_anchored_read_count\tdual_flank_read_count\tmaximum_read_span_bp\trecovered_bp\trecovery_status\tconfidence\tfailure_reason\twarning\nL1\tTXF1\tchr1\t0\t10\t10\t100\t90\t0.9\t0.9\t3\t2\t1\t1000\t20\tpartially_resolved\tmedium\t\tbounded_candidate\n")
    write_report_data(tmp_path)
    report_root = tmp_path / "validated_report"
    shutil.copytree(tmp_path / "families", report_root / "families")
    shutil.copytree(tmp_path / "recovery", report_root / "recovery")
    shutil.copy(tmp_path / "summary.tsv", report_root / "summary.tsv")
    _write(report_root / "figures/family_hierarchy.tsv", "family_id\tvalue\nTXF1\t1\n")
    assert validate_project(report_root)
    data = build_report_data(tmp_path)
    assert data["summary"]["recovery_status_counts"] == {"partially_resolved": 1}
    assert data["summary"]["recovery_candidate_count"] == 1

    invalid = tmp_path / "invalid" / "discover" / "families.tsv"
    _write(invalid, (tmp_path / "families" / "families.tsv").read_text())
    with pytest.raises(ValidationError, match="missing required field"):
        validate_project(invalid.parent.parent)
