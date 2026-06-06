from __future__ import annotations

from pathlib import Path

import pytest

from tandemx.simulate.toy import (
    ToySimulationConfig,
    generate_toy_dataset,
    parse_int_list,
    reverse_complement,
)


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    sequence_parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current is not None:
                lengths[current] = len("".join(sequence_parts))
            current = line[1:]
            sequence_parts = []
        else:
            sequence_parts.append(line.strip())
    if current is not None:
        lengths[current] = len("".join(sequence_parts))
    return lengths


def test_reverse_complement() -> None:
    assert reverse_complement("ACGTN".replace("N", "A")) == "TACGT"


def test_parse_int_list() -> None:
    assert parse_int_list("566,350", "--monomer-lengths") == (566, 350)
    with pytest.raises(ValueError, match="positive"):
        parse_int_list("566,0", "--monomer-lengths")


def test_generate_toy_dataset_outputs_truth_files(tmp_path: Path) -> None:
    config = ToySimulationConfig(
        outdir=tmp_path,
        seed=11,
        num_reads=8,
        read_length=900,
        background_length=1200,
        monomer_lengths=(566, 350),
        copies=(6, 4),
        error_rate=0.01,
    )

    generate_toy_dataset(config)

    expected = {
        "reads.fa",
        "assembly.fa",
        "truth_monomers.fa",
        "truth_arrays.bed",
        "truth_copy_number.tsv",
        "simulation_config.yaml",
    }
    assert expected == {path.name for path in tmp_path.iterdir()}

    monomer_lengths = read_fasta_lengths(tmp_path / "truth_monomers.fa")
    assert sorted(monomer_lengths.values()) == [350, 566]

    copy_number = (tmp_path / "truth_copy_number.tsv").read_text(encoding="utf-8")
    assert "family_id\tmonomer_id\tmonomer_length_bp" in copy_number
    assert "TXF000001\tTXM000001\t566\t6\t2\t3396\t1132\tsimulated_under_assembly" in copy_number
    assert "TXF000002\tTXM000002\t350\t4\t4\t1400\t1400\ttruth_like" in copy_number

    arrays = (tmp_path / "truth_arrays.bed").read_text(encoding="utf-8").splitlines()
    assert len(arrays) == 2
    first = arrays[0].split("\t")
    assert first[0] == "toy_chr1"
    assert first[3] == "TXF000001"
    assert first[8] == "compressed"

    total_size = sum(path.stat().st_size for path in tmp_path.iterdir())
    assert total_size < 1_000_000
