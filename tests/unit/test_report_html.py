from __future__ import annotations

from pathlib import Path

from tandemx.report_html import write_html_report


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
