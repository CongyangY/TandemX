"""Independent checks of the exact engineered competitor control."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.build_assembly_hor_comparator_v1 import build


def read_fasta(path: Path) -> tuple[str, str]:
    lines = path.read_text().splitlines()
    assert len(lines) == 2 and lines[0].startswith(">")
    return lines[0][1:], lines[1]


def test_engineered_copy_truth_matches_bases_and_coordinates(tmp_path: Path) -> None:
    output = tmp_path / "control"
    receipt = build(output)
    identifier, assembly = read_fasta(output / "control.chr1.fasta")
    assert identifier == f"control:0-{len(assembly)}"
    assert hashlib.sha256((output / "control.chr1.fasta").read_bytes()).hexdigest() == receipt["input_sha256"]
    monomer_text = (output / "candidate_monomers.fa").read_text().splitlines()
    monomers = dict(zip((line[1:] for line in monomer_text[::2]), monomer_text[1::2], strict=True))
    with (output / "truth_monomers.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 201
    assert rows[0]["start0"] == "1000"
    assert int(rows[-1]["end0"]) == receipt["array_interval_0based_halfopen"][1]
    for previous, row in enumerate(rows):
        start, end = int(row["start0"]), int(row["end0"])
        assert end - start == 171
        if previous:
            assert start == int(rows[previous - 1]["end0"])
        emitted = assembly[start:end]
        source = monomers[row["label"]]
        assert sum(a != b for a, b in zip(emitted, source, strict=True)) == 1
    with (output / "truth_hors.tsv").open() as handle:
        hors = list(csv.DictReader(handle, delimiter="\t"))
    assert len(hors) == 40
    variants = [row for row in hors if row["canonical"] == "0"]
    assert [int(row["hor_index"]) for row in variants] == [7, 18, 29, 35]
    assert [row["label_path"] for row in variants] == ["ABCE", "ABCCDE", "ABDCE", "ABCDEE"]
    assert json.loads((output / "receipt.json").read_text())["truth_type"] == "exact_engineered_only"
