from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from tandemx.annotation import annotate_repeat_catalog
from tandemx.simulate.toy import reverse_complement


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_catalog(path: Path, sequence: str, family_id: str = "TXF000001") -> None:
    path.write_text(
        f">family_id={family_id};monomer_id=TXM000001;length_bp={len(sequence)};confidence=high\n"
        f"{sequence}\n",
        encoding="utf-8",
    )


def test_annotate_repeats_reports_strong_known_match(tmp_path: Path) -> None:
    sequence = "ACGTTCAGGACTAACCGTGA" * 8
    catalog = tmp_path / "monomers.fa"
    known = tmp_path / "known.fa"
    output = tmp_path / "repeat_annotation.tsv"
    write_catalog(catalog, sequence)
    known.write_text(f">known_repeat\n{sequence}\n", encoding="utf-8")

    annotations = annotate_repeat_catalog(catalog, known, output, kmer_size=11)

    rows = read_tsv(output)
    assert len(annotations) == 1
    assert rows[0]["family_id"] == "TXF000001"
    assert rows[0]["best_known_id"] == "known_repeat"
    assert rows[0]["best_orientation"] == "forward"
    assert rows[0]["annotation_status"] == "strong_known_match"
    assert float(rows[0]["dice"]) == 1.0


def test_annotate_repeats_matches_reverse_complement_known_repeat(tmp_path: Path) -> None:
    sequence = "ACGTTCAGGACTAACCGTGA" * 8
    catalog = tmp_path / "monomers.fa"
    known = tmp_path / "known.fa"
    output = tmp_path / "repeat_annotation.tsv"
    write_catalog(catalog, sequence)
    known.write_text(f">known_reverse\n{reverse_complement(sequence)}\n", encoding="utf-8")

    annotate_repeat_catalog(catalog, known, output, kmer_size=11)

    row = read_tsv(output)[0]
    assert row["best_known_id"] == "known_reverse"
    assert row["best_orientation"] == "reverse"
    assert row["annotation_status"] == "strong_known_match"


def test_annotate_repeats_reports_no_known_match(tmp_path: Path) -> None:
    catalog = tmp_path / "monomers.fa"
    known = tmp_path / "known.fa"
    output = tmp_path / "repeat_annotation.tsv"
    write_catalog(catalog, "ACGTTCAGGACTAACCGTGA" * 8)
    known.write_text(f">unrelated\n{'T' * 80 + 'G' * 80}\n", encoding="utf-8")

    annotate_repeat_catalog(catalog, known, output, kmer_size=11)

    row = read_tsv(output)[0]
    assert row["annotation_status"] == "no_known_match"
    assert row["best_known_id"] == "unrelated"


def test_annotate_repeats_cli_writes_run_files(tmp_path: Path) -> None:
    sequence = "ACGTTCAGGACTAACCGTGA" * 8
    catalog = tmp_path / "monomers.fa"
    known = tmp_path / "known.fa"
    output = tmp_path / "out" / "repeat_annotation.tsv"
    write_catalog(catalog, sequence)
    known.write_text(f">known_repeat\n{sequence}\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "tandemx.cli",
            "annotate-repeats",
            "--catalog",
            str(catalog),
            "--known",
            str(known),
            "--out",
            str(output),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
    assert output.is_file()
    assert (output.parent / "run_config.yaml").is_file()
    assert (output.parent / "run.log").is_file()
