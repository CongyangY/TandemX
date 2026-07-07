from __future__ import annotations

from collections import Counter
from pathlib import Path

from tandemx.quantify.mvp import count_selected_read_kmers_and_bases


def test_count_selected_read_kmers_and_bases_merges_multiple_files(tmp_path: Path) -> None:
    first = tmp_path / "first.fa"
    second = tmp_path / "second.fa"
    first.write_text(">r1\nACGTAC\n", encoding="utf-8")
    second.write_text(">r2\nGTACGT\n", encoding="utf-8")

    counts, total_bases, read_count, max_read_len = count_selected_read_kmers_and_bases(
        (first, second),
        4,
        {"ACGT", "CGTA", "GTAC"},
        "python",
    )

    assert counts == Counter({"ACGT": 2, "CGTA": 2, "GTAC": 2})
    assert total_bases == 12
    assert read_count == 2
    assert max_read_len == 6
