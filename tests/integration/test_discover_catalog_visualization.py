from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_render_discover_catalog_report_outputs_figures_and_summary(tmp_path: Path) -> None:
    families = tmp_path / "families.tsv"
    family_similarity = tmp_path / "family_similarity.tsv"
    monomers = tmp_path / "monomers.fa"
    outdir = tmp_path / "discover_viz"

    families.write_text(
        "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\tsupport_read_count\t"
        "support_span_bp\tmean_identity\tlow_complexity_flag\tconfidence\twarning\n"
        "TXF000001\tTXM000001\t120\tmd5a\t0.4500\t10\t2400\t0.9500\tfalse\thigh\tpossible_higher_order_or_partial:TXF000001-TXF000002\n"
        "TXF000002\tTXM000002\t360\tmd5b\t0.4700\t5\t1200\t0.9200\tfalse\tmedium\tpossible_higher_order_or_partial:TXF000001-TXF000002\n"
        "TXF000003\tTXM000003\t6\tmd5c\t0.3300\t8\t600\t0.8800\ttrue\tmedium\tlow_complexity_family\n",
        encoding="utf-8",
    )
    family_similarity.write_text(
        "family_a\tfamily_b\tlength_a_bp\tlength_b_bp\tkmer_jaccard\tshared_kmer_fraction\tlocal_identity\t"
        "local_overlap_bp\tlocal_overlap_fraction_shorter\tlength_ratio\torientation\trelationship\t"
        "redundant_candidate\tnotes\n"
        "TXF000001\tTXF000002\t120\t360\t0.3200\t0.7000\t0.9300\t120\t1.0000\t3.0000\tforward\t"
        "possible_higher_order_or_partial\tfalse\tPossible higher-order relationship.\n"
        "TXF000001\tTXF000003\t120\t6\t0.0000\t0.0000\t0.5000\t6\t1.0000\t20.0000\tforward\tdistinct\tfalse\t\n",
        encoding="utf-8",
    )
    monomers.write_text(
        ">family_id=TXF000001;monomer_id=TXM000001;length_bp=120;confidence=high\n"
        + "ACGT" * 30
        + "\n>family_id=TXF000002;monomer_id=TXM000002;length_bp=360;confidence=medium\n"
        + "ACGT" * 90
        + "\n>family_id=TXF000003;monomer_id=TXM000003;length_bp=6;confidence=medium\nACGTAC\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "benchmarks/scripts/render_discover_catalog_report.py",
            "--families",
            str(families),
            "--family-similarity",
            str(family_similarity),
            "--monomers",
            str(monomers),
            "--outdir",
            str(outdir),
        ],
        cwd=Path(__file__).resolve().parents[2],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    for stem in (
        "discover_family_abundance",
        "discover_length_distribution",
        "discover_abundance_vs_length",
        "discover_quality_overview",
        "discover_similarity_space",
        "discover_flagged_pairs_heatmap",
    ):
        assert (outdir / f"{stem}.pdf").is_file()
    assert (outdir / "top_families.tsv").is_file()
    assert (outdir / "flagged_pairs.tsv").is_file()
    assert (outdir / "family_label_map.tsv").is_file()
    report = (outdir / "discover_summary.md").read_text(encoding="utf-8")
    assert "TXF000001-TXF000002" in report
