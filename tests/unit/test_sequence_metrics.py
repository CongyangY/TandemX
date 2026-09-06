import random

import pytest

from benchmarks.challenge.sequence_metrics import cyclic_edit_similarity, score_cyclic_recovery
from benchmarks.challenge.evaluate import score_base_coverage
from benchmarks.challenge.schema import ArrayRecord


def distance_reference(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[-1] + 1,
                               previous[j - 1] + int(x != y or x == "N")))
        previous = current
    return previous[-1]


def test_edlib_circular_metric_against_exhaustive_independent_dp():
    rng = random.Random(181)
    for _ in range(40):
        a = "".join(rng.choices("ACGTN", k=rng.randint(1, 15)))
        b = "".join(rng.choices("ACGTN", k=rng.randint(1, 15)))
        reverse = b.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        distance = min(distance_reference(a, s[i:] + s[:i]) for s in (b, reverse) for i in range(len(s)))
        assert cyclic_edit_similarity(a, b) == pytest.approx(1 - distance / max(len(a), len(b)))


def test_indel_rotation_and_unknown_base_scoring():
    sequence = "ACGCTTGACCCATAG"
    mutated = sequence[:6] + "A" + sequence[6:]
    assert cyclic_edit_similarity(sequence, mutated[5:] + mutated[:5]) == pytest.approx(15 / 16)
    assert cyclic_edit_similarity("NNN", "NNN") == 0
    assert cyclic_edit_similarity(sequence, "") == 0
    with pytest.raises(ValueError):
        cyclic_edit_similarity(sequence, "invalid")


def test_one_prediction_cannot_recover_two_planted_units():
    seq = "ACGCTTGACCCATAG"
    scores, _ = score_cyclic_recovery([seq, seq[2:] + seq[:2]], {"a": seq, "b": seq})
    assert scores["cyclic_monomer_recall"] == 0.5
    assert scores["homologous_consensus_fraction"] == 1
    assert scores["distinct_consensus_count"] == 1


def test_union_metrics_distinguish_duplicate_calls_from_wrong_bases():
    truth = [ArrayRecord("r", 10, 110, 20)]
    predictions = [truth[0], truth[0], ArrayRecord("r", 100, 130, 20), ArrayRecord("n", 0, 10, 3)]
    metrics = score_base_coverage(predictions, truth)
    assert metrics["base_union_recall"] == 1
    assert metrics["base_union_precision"] == pytest.approx(100 / 130)
    assert metrics["duplicated_prediction_bp"] == 110
