import pytest
import hashlib
import json

from benchmarks.scripts.score_rice_native_context import parse_paf, qualifies, score


def paf(start: int, end: int, match: int = 9000, columns: int = 9000, tag: str = "tp:A:P") -> str:
    return f"read1\t15000\t100\t9100\t+\tchr1:0-11000\t11000\t{start}\t{end}\t{match}\t{columns}\t60\t{tag}"


def test_requires_primary_identity_and_both_natural_flanks() -> None:
    good = parse_paf(paf(1999, 9001))
    assert qualifies(good, "chr1:0-11000", 3000, 8000, 1000, 0.99)
    assert not qualifies(parse_paf(paf(2001, 9001)), "chr1:0-11000", 3000, 8000, 1000, 0.99)
    assert not qualifies(parse_paf(paf(1999, 8999)), "chr1:0-11000", 3000, 8000, 1000, 0.99)
    assert not qualifies(parse_paf(paf(1999, 9001, 8900, 9000)), "chr1:0-11000", 3000, 8000, 1000, 0.99)
    assert not qualifies(parse_paf(paf(1999, 9001, tag="tp:A:S")), "chr1:0-11000", 3000, 8000, 1000, 0.99)


def test_rejects_malformed_coordinates() -> None:
    with pytest.raises(ValueError, match="target coordinates"):
        parse_paf(paf(-1, 9000))


def test_distinct_spanner_gate_deduplicates_split_alignments(tmp_path) -> None:
    protocol = tmp_path / "protocol.json"
    selected = tmp_path / "selected.json"
    context = tmp_path / "context.fa"
    paf_path = tmp_path / "reads.paf"
    context.write_text(">chr1:0-11000 array=3000-8000 period=155\n" + "A" * 11000 + "\n")
    protocol.write_text(json.dumps({
        "candidate_generation_before_any_read_alignment": {"natural_flank_bp_each_side": 3000},
        "read_alignment_only_after_selected_context_committed": {
            "spanner_min_natural_flank_bp_each_side": 1000,
            "spanner_min_identity_nmatch_over_alignment_columns": 0.99,
            "spanner_min_distinct_original_records": 3,
        },
    }))
    selected.write_text(json.dumps({
        "selected_context.fa_sha256": hashlib.sha256(context.read_bytes()).hexdigest(),
        "selected": {"array_start_0": 3000, "array_end_0": 8000},
    }))
    paf_path.write_text(paf(0, 11000) + "\n" + paf(0, 11000) + "\n" + paf(0, 11000).replace("read1", "read2") + "\n")
    result = score(protocol, selected, context, paf_path, tmp_path / "out")
    assert result["qualifying_distinct_spanners"] == 2
    assert result["pass_for_full_reference_step"] is False
    with paf_path.open("a") as stream:
        stream.write(paf(0, 11000).replace("read1", "read3") + "\n")
    result = score(protocol, selected, context, paf_path, tmp_path / "out")
    assert result["qualifying_distinct_spanners"] == 3
    assert result["pass_for_full_reference_step"] is True
