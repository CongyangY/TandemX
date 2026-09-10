from __future__ import annotations

import json
from pathlib import Path

from tandemx.report_html import _recommended_cards, _summary_cards, write_html_report


def test_html_report_is_standalone_and_escapes_family_values(tmp_path: Path) -> None:
    discover = tmp_path / "discover"; discover.mkdir()
    (discover / "families.tsv").write_text("family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\nX<script>\tM1\t10\tx\t0.5\t1\t10\t1\tfalse\tlow\t<&\n")
    (discover / "monomers.fa").write_text(">family_id=X<script>;monomer_id=M1\nACGT\n")
    write_html_report(tmp_path)
    report = (tmp_path / "report.html").read_text()
    assert "X&lt;script&gt;" in report
    assert "Filter families" in report
    assert (tmp_path / "summary.json").is_file()
    assert report.index("Evidence overview") < report.index("Full family table")
    assert "families recommended for evidence review" in report.lower()
    assert "possible under-representation" in report
    assert "canonical possible-collapse" not in report


def test_html_report_shows_final_recovery_assessment_without_upgrading_partial_candidate(tmp_path: Path) -> None:
    discover = tmp_path / "discover"; discover.mkdir()
    (discover / "families.tsv").write_text("family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\nF1\tM1\t10\tx\t0.5\t1\t10\t1\tfalse\tlow\t\n")
    (discover / "monomers.fa").write_text(">family_id=F1;monomer_id=M1\nACGT\n")
    recovery = tmp_path / "recovery"; recovery.mkdir()
    (recovery / "recovery_candidates.tsv").write_text("recovery_status\npartially_resolved\n")
    (recovery / "recovery_assessment.json").write_text(json.dumps({"status": "assessment_complete", "validated_repeat_gain": False, "decision": "stop_recovery_expansion_negative_poc", "summary": "No validated repeat gain."}))

    write_html_report(tmp_path)
    report = (tmp_path / "report.html").read_text()
    assert "No validated repeat gain." in report
    assert "partially resolved is not a resolved assembly result" in report


def test_recommendation_cards_use_biological_language_not_machine_status() -> None:
    cards = _recommended_cards({"recommended_families": [{
        "family_id": "F1", "comparison_status": "possible_collapse", "confidence": "medium",
        "abundance_deficit_bp": 300,
    }]})
    assert "Possible under-representation" in cards
    assert "confidence: medium" in cards
    assert "possible_collapse" not in cards


def test_summary_cards_show_candidate_relationships_and_only_existing_recovery_assessment_counts() -> None:
    cards = _summary_cards({
        "family_count": 1, "high_confidence_family_count": 1, "possible_underrepresented_count": 0,
        "unresolved_family_count": 0, "low_confidence_family_count": 0, "architecture_edge_count": 2,
        "recovery_assessment": {"candidate_sequences": 1, "validated_recovery_successes": 0},
    })
    assert "candidate family relationships" in cards
    assert "recovery candidate sequences" in cards and ">1<" in cards
    assert "validated recovery successes" in cards and ">0<" in cards
    without_assessment = _summary_cards({"family_count": 1, "architecture_edge_count": 0})
    assert "validated recovery successes" not in without_assessment
