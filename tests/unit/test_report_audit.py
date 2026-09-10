"""Independent semantic checks for report inputs and exported GraphML."""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from tandemx.report_data import build_report_data, write_report_data


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _fixture(root: Path) -> None:
    _write(root / "discover/families.tsv", "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\nF1\tM1\t10\tx\t0.5\t2\t20\t1\tfalse\thigh\t\n")
    _write(root / "discover/monomers.fa", ">family_id=F1;monomer_id=M1;length_bp=10;confidence=high\nACGTACGTAC\n")
    _write(root / "discover/candidate_reads.tsv", "read_id\tcandidate_id\nread\tcandidate\n")
    _write(root / "discover/family_hierarchy.tsv", "hierarchy_edge_id\tshorter_family_id\tlonger_family_id\tshorter_length_bp\tlonger_length_bp\tnearest_integer_multiple\tlength_ratio\tmultiple_error\tlocal_identity\tlocal_overlap_fraction_shorter\tshared_kmer_fraction\torientation\tedge_type\tstatus\twarning\n")
    _write(root / "quantify/copy_number.tsv", "family_id\testimated_bp\tconfidence\twarning\nF1\t100\tlow\t\n")
    _write(root / "locate/arrays.bed", "")
    _write(root / "compare/assembly_vs_read_cn.tsv", "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\nF1\t100\t20\t0.2\tpossible_collapse\tmedium\t\n")


def test_composite_confidence_and_graphml_are_conservative_and_parseable(tmp_path: Path) -> None:
    _fixture(tmp_path)
    data = write_report_data(tmp_path)
    family = data["families"][0]
    assert family["confidence"] == "low"
    assert family["abundance_deficit_bp"] == 80
    assert family["candidate_array_count"] == 0
    graph = ET.parse(tmp_path / "families/family_network.graphml")
    assert graph.getroot().tag.endswith("graphml")


def test_failed_stage_records_exclude_existing_stage_files(tmp_path: Path) -> None:
    _fixture(tmp_path)
    records = [
        {"step": "discover", "exit_status": 0, "output_validated": True, "notes": ""},
        {"step": "quantify", "exit_status": 1, "output_validated": False, "notes": "failed"},
        {"step": "locate", "exit_status": 1, "output_validated": False, "notes": "failed"},
        {"step": "compare", "exit_status": 1, "output_validated": False, "notes": "failed"},
    ]
    family = build_report_data(tmp_path, records)["families"][0]
    assert family["estimated_abundance_bp"] is None
    assert family["assembly_representation_bp"] is None
    assert family["candidate_array_count"] is None
