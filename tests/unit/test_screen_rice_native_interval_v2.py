import csv
import hashlib

import pytest

from benchmarks.scripts.screen_rice_native_interval_v2 import (
    rank_key,
    score_periods,
    seeds_from_v1,
)


def test_all_periods_scored_and_exact_repeat_wins() -> None:
    scores = list(score_periods(b"AACCGGTT" * 5, 3, 12))
    assert len(scores) == 10
    assert next(x for x in scores if x[0] == 8)[2] == 1.0


def test_ranking_prioritizes_flank_uniqueness() -> None:
    low_unique_high_identity = ("a", 0, 1000, 100, 899, 0.999, 1.9, 0.3, True)
    high_unique_lower_identity = ("b", 0, 1000, 150, 850, 0.85, 1.9, 0.8, True)
    assert sorted([low_unique_high_identity, high_unique_lower_identity], key=rank_key)[0] == high_unique_lower_identity


def test_seed_table_hash_and_negative_history_enforced(tmp_path) -> None:
    path = tmp_path / "v1.tsv"
    with path.open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "start_0", "end_0", "shift_identity", "flank_unique_31mer_fraction", "selected"])
        for i in range(1000):
            writer.writerow(["chr1", i * 500, i * 500 + 4000, "0.9", f"{i/1000:.3f}", "false"])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    seeds = seeds_from_v1(path, digest, 2)
    assert [int(x["start_0"]) for x in seeds] == [499500, 499000]
    with pytest.raises(ValueError, match="SHA-256"):
        seeds_from_v1(path, "0" * 64, 2)
