import random

import pytest

from tandemx.discover.rust_backend import rust_backend_available, seed_spacing_histogram
from tandemx.discover.spacing import extract_repeated_kmer_positions, build_spacing_histogram


@pytest.mark.skipif(not rust_backend_available(), reason="compiled extension unavailable")
def test_native_histogram_matches_bounded_reference_over_parameters():
    rng = random.Random(287831)
    sequences = ["", "N" * 100, "A" * 1000, "AT" * 500]
    for _ in range(45):
        unit = "".join(rng.choices("ACGT", k=rng.randint(13, 183)))
        sequences.append((unit * rng.randint(2, 12) + "N" + unit.lower() * 3))
        sequences.append("".join(rng.choices("ACGTN", k=rng.randint(1, 700))))
    for sequence in sequences:
        for k in (3, 7, 11, 31):
            minimum, maximum = rng.choice([(20, 1000), (20, 19), (30, 31), (1, 10**9)])
            occurrences, cap = rng.choice([(2, 100), (4, 5), (2, 1)])
            positions, overflow = extract_repeated_kmer_positions(sequence, k, occurrences, cap)
            histogram = build_spacing_histogram(positions, minimum, maximum, cap)
            assert seed_spacing_histogram(sequence, k, minimum, maximum, occurrences, cap) == (histogram, overflow)


@pytest.mark.skipif(not rust_backend_available(), reason="compiled extension unavailable")
def test_native_histogram_invalid_parameters_are_explicit():
    for k in (0, 32):
        with pytest.raises(ValueError):
            seed_spacing_histogram("ACGT", k, 20, 100, 2, 100)
    with pytest.raises(ValueError):
        seed_spacing_histogram("DNAé", 3, 20, 100, 2, 100)
