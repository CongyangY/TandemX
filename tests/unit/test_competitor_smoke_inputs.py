"""Guard the tiny comparator fixture against silent truth drift."""

from pathlib import Path
import subprocess
import sys


def test_competitor_smoke_fixture_is_deterministic(tmp_path: Path) -> None:
    script = (
        Path(__file__).resolve().parents[2]
        / "benchmarks/scripts/generate_competitor_smoke_inputs.py"
    )
    out = tmp_path / "inputs"
    subprocess.run([sys.executable, str(script), str(out)], check=True)
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    subprocess.run([sys.executable, str(script), str(out)], check=True)
    assert first == {p.name: p.read_bytes() for p in out.iterdir()}

    a = "ACGTGCAATGTCAGTACCGTACGATCGTTA"
    b = "TTGACCGATGCTAGACCTGAGTCATCGTAC"
    assert (out / "cendetecthor_ab.fa").read_text() == ">toy.chr1\n" + (a + b) * 12 + "\n"
    assert (out / "cendetecthor_pipeline.fa").read_text() == (
        ">toy:0-5700\n" + ((a + b) * 9 + a) * 10 + "\n"
    )
    bed = (out / "cendetecthor_ab.bed").read_text().splitlines()
    assert len(bed) == 24
    assert bed[0].split("\t")[:3] == ["toy.chr1", "0", "30"]
    assert bed[-1].split("\t")[:3] == ["toy.chr1", "690", "720"]
