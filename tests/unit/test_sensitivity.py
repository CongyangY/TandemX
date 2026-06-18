from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from tandemx.sensitivity import check_known_repeats
from tandemx.simulate.toy import reverse_complement


def test_known_repeat_checker_matches_forward_and_reverse(tmp_path: Path) -> None:
    catalog = tmp_path / "monomers.fa"
    known = tmp_path / "known.fa"
    output = tmp_path / "known_repeat_matches.tsv"
    sequence = "AACCGTTAGGCTACGATTCG"
    catalog.write_text(
        f">family_id=TXF000001;monomer_id=TXM000001;length_bp=20;confidence=high\n{sequence}\n",
        encoding="utf-8",
    )
    known.write_text(
        f">known_forward\n{sequence}\n>known_reverse\n{reverse_complement(sequence)}\n",
        encoding="utf-8",
    )

    matches = check_known_repeats(catalog, known, output, kmer_size=7)

    assert [match.best_family_id for match in matches] == ["TXF000001", "TXF000001"]
    assert matches[0].orientation == "forward"
    assert matches[1].orientation == "reverse"
    assert all(match.similarity_score == 1.0 for match in matches)
    with output.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["interpretation"] == "strong_post_hoc_match"
    assert rows[1]["shared_kmer_fraction"] == "1.000000"


def test_known_repeat_checker_script(tmp_path: Path) -> None:
    catalog = tmp_path / "monomers.fa"
    known = tmp_path / "known.fa"
    output = tmp_path / "matches.tsv"
    catalog.write_text(
        ">family_id=TXF000009;monomer_id=TXM000009;length_bp=12;confidence=high\nAACCGGTTACGA\n",
        encoding="utf-8",
    )
    known.write_text(">known9\nAACCGGTTACGA\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            "benchmarks/scripts/check_known_repeats_against_catalog.py",
            "--catalog",
            str(catalog),
            "--known",
            str(known),
            "--out",
            str(output),
            "--kmer-size",
            "5",
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert output.is_file()
    assert "wrote 1 post hoc matches" in result.stdout
