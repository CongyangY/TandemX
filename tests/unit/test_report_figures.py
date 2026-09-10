import json
from pathlib import Path

from tandemx.report_figures import write_report_figures


def _data():
    return {
        "schema_version": "report.v1",
        "run": {"sample": "toy"},
        "summary": {"family_count": 4},
        "families": [
            {"family_id": "TXF000001", "monomer_length_bp": 171, "gc_fraction": .5, "support_read_count": 8, "candidate_array_count": 1, "estimated_abundance_bp": 1_000_000, "assembly_representation_bp": 800_000, "assembly_read_ratio": .8, "abundance_deficit_bp": 200_000, "confidence": "high", "warning": ""},
            {"family_id": "TXF000002", "monomer_length_bp": 342, "gc_fraction": .6, "support_read_count": 4, "candidate_array_count": 1, "estimated_abundance_bp": 0, "assembly_representation_bp": 0, "assembly_read_ratio": None, "abundance_deficit_bp": 0, "confidence": "unresolved", "warning": "zero_estimate"},
            {"family_id": "TXF000003", "monomer_length_bp": 684, "gc_fraction": None, "support_read_count": 2, "candidate_array_count": 0, "estimated_abundance_bp": 2_000_000_000, "assembly_representation_bp": 10_000, "assembly_read_ratio": .000005, "abundance_deficit_bp": 1_999_990_000, "confidence": "medium", "warning": "assembly_missing"},
            {"family_id": "TXF000004", "monomer_length_bp": None, "gc_fraction": None, "support_read_count": None, "candidate_array_count": None, "estimated_abundance_bp": None, "assembly_representation_bp": None, "assembly_read_ratio": None, "abundance_deficit_bp": None, "confidence": "unresolved", "warning": "metadata_blocked"},
        ],
        "architecture_edges": [{"hierarchy_edge_id": "TXH1", "shorter_family_id": "TXF000001", "longer_family_id": "TXF000002", "length_ratio": 2.0, "edge_type": "putative_period_multiple", "status": "candidate", "warning": "heuristic_period_multiple_not_validated_hor"}],
    }


def test_write_report_figures_exports_all_formats_and_audit_files(tmp_path: Path):
    result = write_report_figures(tmp_path, _data())
    for name in ("summary", "family_abundance_vs_assembly", "top_underrepresented", "family_landscape", "family_hierarchy", "family_evidence_cards"):
        for suffix in ("svg", "pdf", "png", "tsv", "receipt.json"):
            assert (tmp_path / "figures" / f"{name}.{suffix}").is_file()
        receipt = json.loads((tmp_path / "figures" / f"{name}.receipt.json").read_text())
        assert receipt["svg_node_counts"]["text"] > 0
        assert receipt["svg_node_counts"]["raster"] == 0
        assert f"{name}_svg" in result
    top_rows = (tmp_path / "figures" / "top_underrepresented.tsv").read_text().splitlines()
    assert all("TXF000002" not in row for row in top_rows[1:])
    assert "read abundance" in (tmp_path / "figures" / "top_underrepresented.svg").read_text()
    assert "assembly representation" in (tmp_path / "figures" / "top_underrepresented.svg").read_text()
    assert "r=2.00" in (tmp_path / "figures" / "family_hierarchy.svg").read_text()


def test_empty_data_is_a_valid_auditable_report(tmp_path: Path):
    result = write_report_figures(tmp_path, {"schema_version": "report.v1", "families": [], "architecture_edges": []})
    assert len([key for key in result if key.endswith("_svg")]) == 6
    assert "No quantitative data available" in (tmp_path / "figures" / "family_abundance_vs_assembly.svg").read_text()


def test_long_warnings_are_summarized_visibly_but_preserved_in_source(tmp_path: Path):
    data = _data()
    warning = ";".join(["first_warning", "second_warning", "third_warning"] + [f"detail_{i:04d}" for i in range(100)])
    data["families"][0]["warning"] = warning
    write_report_figures(tmp_path, data)
    svg = (tmp_path / "figures" / "family_evidence_cards.svg").read_text()
    source = (tmp_path / "figures" / "family_evidence_cards.tsv").read_text()
    assert "Notes: First warning;Second warning;Third warning (100 further warnings in table)" in svg
    assert warning in source
