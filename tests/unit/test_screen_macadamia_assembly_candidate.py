from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "benchmarks/scripts/screen_macadamia_assembly_candidate.py"
spec = importlib.util.spec_from_file_location("screen_macadamia_assembly_candidate", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_streamed_contexts_and_frozen_rank(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">ctg1 metadata\nAAAA\nCGCG\nCGCG\nTTTT\n>ctg2\nAAAACCCC\n")
    candidates = [
        {"contig": "ctg1", "start": 4, "end": 12, "family_id": "F1", "length_bp": 8, "period_bp": 2, "flank_bp": 4},
        {"contig": "ctg2", "start": 2, "end": 6, "family_id": "F2", "length_bp": 4, "period_bp": 2, "flank_bp": 2},
    ]
    module.extract_contexts(fasta, candidates)
    assert candidates[0]["context"] == "AAAACGCGCGCGTTTT"
    assert candidates[1]["context"] == "AAAACCCC"
    scored = module.score_candidates(candidates)
    assert scored[0]["contig"] == "ctg1"
    assert scored[0]["period_shift_identity"] == 1.0
    assert scored[0]["period_shift_matches"] == 6


def test_truncated_candidate_context_fails_closed(tmp_path: Path) -> None:
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">ctg1\nACGT\n")
    candidates = [{"contig": "ctg1", "start": 1, "end": 4, "family_id": "F1", "length_bp": 3, "period_bp": 1, "flank_bp": 1}]
    with pytest.raises(ValueError, match="truncated FASTA context"):
        module.extract_contexts(fasta, candidates)
