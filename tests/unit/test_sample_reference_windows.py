from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.challenge.adapters import read_fasta
from benchmarks.scripts.sample_reference_windows import sample_reference


def write_reference(path: Path) -> dict[str, str]:
    sequences = {
        "chr1": "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
        "chr2": "ACGT" * 8,
        "chr3": "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA",
    }
    path.write_text(
        "".join(f">{name} descriptive text\n{sequence[:13]}\n{sequence[13:]}\n" for name, sequence in sequences.items())
    )
    return sequences


def test_nested_reference_windows_preserve_source_sequence(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    source = write_reference(reference)
    outdir = tmp_path / "samples"
    receipt = sample_reference(
        reference,
        outdir,
        ["small=16", "large=40"],
        window_size=8,
        seed=17,
    )

    assert receipt["complete"] is True
    assert receipt["reference"]["total_bases"] == 96
    assert receipt["samples"]["small"]["window_count"] == 2
    assert receipt["samples"]["large"]["window_count"] == 5
    small = read_fasta(outdir / "small.fa")
    large = read_fasta(outdir / "large.fa")
    assert set(small) < set(large)
    assert sum(map(len, small.values())) == 16
    assert sum(map(len, large.values())) == 40

    with (outdir / "windows.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 5
    for row in rows:
        rank = int(row["rank"])
        expected_samples = "small,large" if rank <= 2 else "large"
        assert row["sample_ids"] == expected_samples
        expected = source[row["source_sequence"]][
            int(row["source_start"]) : int(row["source_end"])
        ]
        assert large[row["record_id"]] == expected

    persisted = json.loads((outdir / "receipt.json").read_text())
    assert persisted == receipt


def test_reference_sampling_rejects_nonmultiple_scale(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    write_reference(reference)
    with pytest.raises(ValueError, match="not divisible"):
        sample_reference(
            reference,
            tmp_path / "samples",
            ["bad=17"],
            window_size=8,
            seed=1,
        )
