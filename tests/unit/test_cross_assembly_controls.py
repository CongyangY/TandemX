from __future__ import annotations

from collections import Counter
from pathlib import Path

from benchmarks.scripts.build_cross_assembly_single_copy_controls import (
    build_controls,
    candidate_words,
    selection_score,
)
from tandemx.utils.kmers import canonical_kmer, iter_circular_canonical_kmers


def write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    repeated = "ACGTTGCACTGATCGAACCTG"
    background = (
        "GACCTAGTTCGATCGGATCCATGCTAGCATCGTACGATGGCATTCGACCTAGGCTAACGT"
        "TGACCGTATGCCATGATCGTACCTGATGCTAGCATTCGGTACGATCGATGCCATGGTCA"
    )
    old = tmp_path / "old.fa"
    new = tmp_path / "new.fa"
    catalogue = tmp_path / "monomers.fa"
    old.write_text(f">old_a\n{background}\n>old_repeat\n{repeated * 3}\n")
    new.write_text(f">new_a\n{background}\n>new_repeat\n{repeated * 5}\n")
    catalogue.write_text(
        f">family_id=TXF000001;monomer_id=M1;length_bp={len(repeated)};confidence=high\n{repeated}\n"
    )
    return old, new, catalogue


def test_candidates_are_deterministic_and_exclude_catalogue(tmp_path: Path) -> None:
    _, new, catalogue = write_inputs(tmp_path)
    first, audit = candidate_words(
        new, catalogue, k=15, stride=1, candidate_limit=25, seed=6101
    )
    second, _ = candidate_words(
        new, catalogue, k=15, stride=1, candidate_limit=25, seed=6101
    )
    excluded = set(iter_circular_canonical_kmers("ACGTTGCACTGATCGAACCTG", 15))
    assert first == second
    assert first == sorted(first, key=lambda word: (selection_score(word, 6101), word))
    assert not set(first).intersection(excluded)
    assert audit["catalogue_overlap_sampled_words"] > 0


def test_build_controls_requires_one_copy_in_both_assemblies(tmp_path: Path) -> None:
    old, new, catalogue = write_inputs(tmp_path)
    output = tmp_path / "controls.tsv"
    receipt = tmp_path / "receipt.json"
    result = build_controls(
        old,
        new,
        catalogue,
        output,
        receipt,
        k=15,
        desired=5,
        stride=1,
        candidate_limit=80,
        seed=7,
    )
    rows = [line.split("\t") for line in output.read_text().splitlines()[1:]]
    old_sequence = "".join(
        line.strip() for line in old.read_text().splitlines() if not line.startswith(">")
    )
    new_sequence = "".join(
        line.strip() for line in new.read_text().splitlines() if not line.startswith(">")
    )
    old_counts = Counter(
        canonical_kmer(old_sequence[index : index + 15])
        for index in range(len(old_sequence) - 14)
    )
    new_counts = Counter(
        canonical_kmer(new_sequence[index : index + 15])
        for index in range(len(new_sequence) - 14)
    )
    assert len(rows) == 5
    assert all(row[1:4] == ["1", "1", "1"] for row in rows)
    assert all(old_counts[row[0]] == new_counts[row[0]] == 1 for row in rows)
    assert result["selected_controls"] == 5
    assert receipt.is_file()
