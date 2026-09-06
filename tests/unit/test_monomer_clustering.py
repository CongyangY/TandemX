from dataclasses import replace
import random

import edlib
import pytest

from benchmarks.challenge.sequence_metrics import cyclic_edit_similarity
from tandemx.discover.clustering import cluster_monomers, resolve_clustering_method
from tandemx.discover.distance import bounded_edit_distance, cyclic_merge_evidence
from tandemx.discover.mvp import CandidateRepeat, orient_monomer
from tandemx.discover.rust_backend import rust_backend_available


BACKENDS = ["python", pytest.param("rust", marks=pytest.mark.skipif(
    not rust_backend_available(), reason="compiled extension unavailable"))]


def candidate(sequence, index, read_id=None):
    return CandidateRepeat(read_id or f"r{index}", f"TXC{index:06d}", sequence, 0, len(sequence) * 10,
                           "+", len(sequence), len(sequence) * 10, 10, 1.0, False, "medium", "")


def change(sequence, positions):
    bases = list(sequence)
    for i in positions:
        bases[i] = "ACGT"[("ACGT".index(bases[i]) + 1) % 4]
    return "".join(bases)


@pytest.mark.parametrize("backend", BACKENDS)
def test_thresholded_distance_matches_independent_edlib(backend):
    rng = random.Random(73)
    pairs = [("", ""), ("NN", "NN"), ("ACGT", ""), ("aCgT", "acgt")]
    pairs += [("".join(rng.choices("ACGTN", k=rng.randrange(35))),
               "".join(rng.choices("ACGTN", k=rng.randrange(35)))) for _ in range(60)]
    for a, b in pairs:
        distance = edlib.align(a.upper().replace("N", "X"), b.upper().replace("N", "Y"),
                              mode="NW")["editDistance"]
        for limit in (0, 1, 5, 35):
            assert bounded_edit_distance(a, b, limit, backend) == min(limit + 1, distance)
    with pytest.raises(ValueError, match="million"):
        bounded_edit_distance("A" * 100, "A" * 100, 1_000_000, backend)


@pytest.mark.parametrize("backend", BACKENDS)
def test_circular_prefilter_has_no_false_rejections_against_oracle(backend):
    rng = random.Random(978)
    for _ in range(55):
        a = "".join(rng.choices("ACGT", k=rng.randint(5, 180)))
        b = list(a)
        for _ in range(rng.randint(0, max(1, len(a) // 8))):
            position = rng.randrange(len(b))
            operation = rng.randrange(3)
            if operation == 0:
                b[position] = rng.choice("ACGTN")
            elif operation == 1:
                b.insert(position, rng.choice("ACGT"))
            elif len(b) > 1:
                b.pop(position)
        b = "".join(b)
        offset = rng.randrange(len(b))
        b = b[offset:] + b[:offset]
        if rng.randrange(2):
            b = b.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        oracle = cyclic_edit_similarity(a, b)
        for threshold in (0.9, 0.95, 1.0):
            result = cyclic_merge_evidence(a, b, threshold, backend)
            assert (result is not None) == (oracle + 1e-12 >= threshold), (a, b, threshold, oracle)
            if result:
                assert threshold - 1e-12 <= result[1] <= oracle + 1e-12


@pytest.mark.parametrize("backend", BACKENDS)
def test_clusters_keep_related_monomers_and_avoid_transitive_bridges(backend):
    rng = random.Random(956)
    a = "".join(rng.choices("ACGT", k=100))
    b = change(a, [6, 26, 56, 76])
    c = change(b, [16, 36, 66, 86])
    candidates = [candidate(a, i) for i in range(4)]
    candidates += [candidate(b, i) for i in (4, 5)] + [candidate(c, i) for i in (6, 7)]
    families, rows = cluster_monomers(candidates, 2, backend=backend)
    assert len(families) == 2
    assert {f.monomer_sequence for f in families} == {orient_monomer(a), orient_monomer(c)}
    assert sorted(f.support_read_count for f in families) == [2, 6]
    shuffled = candidates[:]
    rng.shuffle(shuffled)
    assert cluster_monomers(shuffled, 2, backend=backend) == (families, rows)
    assert len(cluster_monomers(candidates, 2, 0.9, backend)[0]) == 1


def test_candidate_evidence_survives_support_filters_and_unknown_bases():
    a = "ACGATTCGCCATGTAACTGCGTACATTCGGACGTGCATTCGATCGATCGTGA"
    candidates = [candidate(a, 1), candidate(a[11:] + a[:11], 2, "r1"),
                  candidate("NNNNN", 3), candidate("A" * 30, 4)]
    families, rows = cluster_monomers(candidates, 2)
    assert not families  # two array observations in one read are one supporting read
    assert len(rows) == len(candidates)
    assert rows[0]["status"] == "below_minimum_support"
    assert rows[2]["status"] == "unresolved_sequence"
    assert rows[2]["family_id"] == "NA"
    one_n = replace(candidate(a, 1), sequence=a[:-1] + "N")
    _, rows = cluster_monomers([one_n], 1)
    assert rows[0]["edit_distance_upper_bound"] == 1
    assert rows[0]["similarity_lower_bound"] < 1
    assert cluster_monomers([], 1) == ([], [])
    assert resolve_clustering_method("auto", "elastic") == "sequence"
    assert resolve_clustering_method("auto", "legacy") == "legacy"
    for value in (0, 1.1, float("nan")):
        with pytest.raises(ValueError):
            cluster_monomers(candidates, 1, value)


def test_multiple_compatible_cluster_assignments_are_exposed():
    rng = random.Random(982)
    a = "".join(rng.choices("ACGT", k=100))
    b = change(a, [6, 26, 56, 76])
    c = change(b, [16, 36, 66, 86])
    candidates = [candidate(a, i) for i in range(5)]
    candidates += [candidate(c, i) for i in range(5, 9)] + [candidate(b, 9)]
    _, rows = cluster_monomers(candidates, 1)
    assert rows[-1]["compatible_cluster_count"] == 2
    assert rows[-1]["alternative_cluster_ids"]
    assert rows[-1]["warning"] == "multiple_compatible_clusters"
