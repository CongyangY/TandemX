import copy
import math
import random

import pytest

from tandemx.quantify.joint_moments import JointReadMoments


def test_joint_variance_matches_independent_per_read_influences():
    ks = (15, 21, 27, 31)
    rng = random.Random(9012)
    reads = []
    moments = JointReadMoments(ks, {'f': (10, 10, 10, 10)})
    for i in range(100):
        length = rng.randrange(60, 200)
        amount = rng.uniform(.1, 10) if i % 3 else 0
        # Shared read support and variable loss generate correlated k values.
        loss = rng.uniform(.005, .03)
        y = [amount*math.exp(-loss*k) for k in ks]
        reads.append((length, y))
        moments.add(length, {'f': y} if amount else {})
    estimate = moments.estimates(10000, 2)[0]
    mean_k = sum(ks)/4
    weights = [1/4-mean_k*(k-mean_k)/sum((x-mean_k)**2 for x in ks) for k in ks]
    sx = [sum(length-k+1 for length, y in reads) for k in ks]
    sy = [sum(y[j] for length, y in reads) for j in range(4)]
    influences = [sum(weights[j]*(y[j]/sy[j]-(length-k+1)/sx[j]) for j, k in enumerate(ks))
                  for length, y in reads]
    expected = len(reads)/(len(reads)-1)*sum(g*g for g in influences)
    assert estimate.log_sampling_variance == pytest.approx(expected, rel=1e-11)
    assert estimate.sampling_interval_low < estimate.estimated_copy_number < estimate.sampling_interval_high
    assert estimate.positive_reads_by_k == (66, 66, 66, 66)
    assert 'correlated' in estimate.warning


def test_short_zero_and_missing_diagnostics_remain_explicit():
    moments = JointReadMoments((3, 5, 7), {'absent': (2, 2, 2), 'ambiguous': (0, 0, 0), 'rare': (2, 2, 2)})
    moments.add(2, {})
    moments.add(30, {'rare': (10, 8, 6)})
    by_family = {r.family_id: r for r in moments.estimates(100)}
    assert by_family['absent'].estimated_copy_number is None
    assert 'zero_support' in by_family['absent'].fit_status
    assert by_family['ambiguous'].fit_status == 'missing_diagnostic_or_exposure'
    assert by_family['rare'].interval_status == 'insufficient_effective_reads'
    assert by_family['rare'].sampling_interval_low is None


def test_invalid_read_does_not_partially_mutate_moments():
    moments = JointReadMoments((3, 5, 7), {'f': (1, 1, 1)})
    moments.add(50, {'f': (5, 4, 3)})
    before = copy.deepcopy(moments.__dict__)
    for length, contributions in [(50, {'unknown': (1, 1, 1)}), (50, {'f': (1, 2)}),
                                  (50, {'f': (1, float('nan'), 1)}), (4, {'f': (1, 1, 1)})]:
        with pytest.raises(ValueError):
            moments.add(length, contributions)
        assert moments.__dict__ == before
    with pytest.raises(ValueError, match='increasing'):
        JointReadMoments((7, 3, 5), {'f': (1, 1, 1)})


def test_identical_read_contributions_do_not_create_spurious_tiny_interval():
    moments = JointReadMoments((15, 21, 27, 31), {'f': (10, 10, 10, 10)})
    for _ in range(100):
        moments.add(100, {'f': (10., 9., 8., 7.)})
    estimate = moments.estimates(10000)[0]
    assert estimate.interval_status == 'variance_below_moment_resolution_uncalibrated'
    assert estimate.log_sampling_variance is None
    assert estimate.sampling_interval_low is None
