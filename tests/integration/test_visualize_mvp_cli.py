from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "tandemx.cli", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def test_visualize_mvp_outputs_svg_and_pdf(tmp_path: Path) -> None:
    copy_number = tmp_path / "copy_number.tsv"
    comparison = tmp_path / "assembly_vs_read_cn.tsv"
    probes = tmp_path / "probes.rank.tsv"
    fish = tmp_path / "in_silico_fish.tsv"
    catalog = tmp_path / "monomers.fa"
    outdir = tmp_path / "visualize"

    catalog.write_text(">family_id=TXF000001;length_bp=100\nACGT\n", encoding="utf-8")
    copy_number.write_text(
        "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
        "haploid_depth\testimated_copy_number\testimated_bp\tconfidence\twarning\n"
        "TXF000001\t100\t10\t10\t1\t10\t1000\thigh\t\n",
        encoding="utf-8",
    )
    comparison.write_text(
        "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning\n"
        "TXF000001\t1000\t300\t0.3\tpossible_collapse\tmedium\t\n",
        encoding="utf-8",
    )
    probes.write_text(
        "probe_id\tfamily_id\tsequence_length\tgc_content\ttm\testimated_copy_number\t"
        "arrayiness_score\tspecificity_score\toff_target_hits\tpredicted_regions\tprobe_score\tconfidence\twarning\n"
        "TXP000001\tTXF000001\t100\t0.5\t300\t10\t1\t1\t0\tchr1:10-100\t0.9\thigh\t\n",
        encoding="utf-8",
    )
    fish.write_text(
        "probe_id\tchrom\tstart\tend\tpredicted_signal\tconfidence\twarning\n"
        "TXP000001\tchr1\t10\t100\t0.9\thigh\t\n",
        encoding="utf-8",
    )

    result = run_cli(
        "visualize",
        "--catalog",
        str(catalog),
        "--copy-number",
        str(copy_number),
        "--comparison",
        str(comparison),
        "--probes",
        str(probes),
        "--fish",
        str(fish),
        "--outdir",
        str(outdir),
    )

    assert result.returncode == 0, result.stderr
    for stem in ("catalogue_summary", "assembly_vs_read", "in_silico_fish"):
        assert (outdir / f"{stem}.svg").is_file()
        assert (outdir / f"{stem}.pdf").is_file()
