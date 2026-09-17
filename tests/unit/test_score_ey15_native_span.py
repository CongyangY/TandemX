"""Guard the Ey15 diagnostic against scoring from its edit ledger alone."""

from __future__ import annotations

import pytest

from benchmarks.scripts import score_ey15_native_span as scorer


def test_edited_assembly_span_is_measured_independently_of_receipt(
    monkeypatch: pytest.MonkeyPatch, tmp_path,
) -> None:
    """A changed in-memory FASTA must fail even while the ledger remains frozen."""
    original_reader = scorer.fasta_records

    def spoof_one_edited_sequence(path):
        records = original_reader(path)
        if path.name == "E1_terminal_025.fa":
            sequence = records["E1_terminal_025"]
            records["E1_terminal_025"] = sequence[:4000] + "A" * 20 + sequence[4000:]
        return records

    monkeypatch.setattr(scorer, "fasta_records", spoof_one_edited_sequence)
    with pytest.raises(ValueError, match="Edited assembly span disagrees with exact edit ledger"):
        scorer.run(tmp_path / "score")
