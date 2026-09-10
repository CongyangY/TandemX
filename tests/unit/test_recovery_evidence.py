from pathlib import Path

import pytest

from tandemx.recovery.evidence import (
    Alignment, anchor_uniqueness, extract_intervals, historical_loci, make_flanks,
    read_paf, select_spanning_candidate, spanning_interval,
)


def alignment(query="r", strand="+", start=0, end=2000, target="a"):
    return Alignment(query, 10000, start, end, strand, target, 2000, 0, 2000, 2000, 2000, 60)


def test_stream_extract_and_edge_flanks(tmp_path: Path):
    path = tmp_path / "old.fa"
    path.write_text(">Chr1\nACGT\nACGT\n>Chr2\nTTTT\n")
    intervals = [dict(anchor_id="a", chromosome="Chr1", start=2, end=7),
                 dict(anchor_id="b", chromosome="Chr1", start=-3, end=3),
                 dict(anchor_id="c", chromosome="Chr2", start=1, end=10)]
    assert extract_intervals(path, intervals) == {"a": "GTACG", "c": "TTT"}
    loci = [dict(locus_id="L1", family_id="f", chromosome="Chr1", start=2, end=5)]
    assert len(make_flanks(loci)) == 8
    assert make_flanks(loci)[0]["start"] == -1998


def test_no_locus_and_overlap_union(tmp_path: Path):
    bed = tmp_path / "arrays.bed"
    bed.write_text("Chr1\t1\t5\tf\nChr1\t3\t7\tf\nChr1\t10\t12\tg\n")
    assert historical_loci(bed, {"missing"}) == []
    rows = historical_loci(bed, {"f"})
    assert len(rows) == 1 and (rows[0]["start"], rows[0]["end"]) == (1, 7)
    bed.write_text("Chr1\t-1\t2\tf\n")
    with pytest.raises(ValueError):
        historical_loci(bed, {"f"})


def test_uniqueness_requires_expected_and_no_competitor():
    anchor = dict(chromosome="Chr1", start=100, end=2100)
    hit = Alignment("a", 2000, 0, 2000, "+", "Chr1", 10000, 100, 2100, 2000, 2000, 60)
    assert anchor_uniqueness(anchor, "ACGT"*500, [hit])["eligible"]
    competitor = Alignment("a", 2000, 0, 2000, "+", "Chr2", 10000, 200, 2200, 1950, 2000, 0)
    assert not anchor_uniqueness(anchor, "ACGT"*500, [hit, competitor])["eligible"]
    assert not anchor_uniqueness(anchor, "N"*2000, [hit])["eligible"]
    assert not anchor_uniqueness(anchor, "ACGT", [hit])["eligible"]
    assert not anchor_uniqueness(anchor, "ACGT"*500, [])["eligible"]


def test_spans_require_same_read_strand_order():
    assert spanning_interval(alignment(), alignment(start=5000, end=7000)) == (2000, 5000, "+")
    assert spanning_interval(alignment(strand="-", start=5000, end=7000), alignment(strand="-")) == (2000, 5000, "-")
    assert spanning_interval(alignment(), alignment(query="other", start=5000, end=7000)) is None
    assert spanning_interval(alignment(), alignment(strand="-", start=5000, end=7000)) is None
    assert spanning_interval(alignment(start=5000, end=7000), alignment()) is None


def test_three_distinct_reads_and_discordant_lengths():
    rows = [dict(read_id=f"r{i}", start=1000, end=6000+i, strand="+") for i in range(3)]
    chosen, status = select_spanning_candidate(rows)
    assert status == "partially_resolved" and chosen == rows[1]
    assert select_spanning_candidate(rows[:1]*3)[1] == "insufficient_read_support"
    assert select_spanning_candidate(rows + [dict(rows[0], end=8000)])[1] == "unresolved_conflicting_paths"
    assert select_spanning_candidate(rows + [dict(rows[0], read_id="fourth", end=8000)])[1] == "unresolved_conflicting_paths"


def test_paf_empty_valid_and_invalid(tmp_path: Path):
    path = tmp_path / "alignment.paf"
    path.write_text("")
    assert list(read_paf(path)) == []
    path.write_text("a\t2000\t0\t2000\t+\tChr1\t10000\t100\t2100\t2000\t2000\t60\n")
    assert list(read_paf(path))[0].identity == 1
    path.write_text("a\t2000\t-1\t2000\t+\tChr1\t10000\t100\t2100\t2000\t2000\t60\n")
    with pytest.raises(ValueError):
        list(read_paf(path))
