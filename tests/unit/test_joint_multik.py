from collections import defaultdict
import random

import pytest

from tandemx.quantify.joint_multik import estimate_joint_multik
from tandemx.quantify.multik import estimate_multik
from tandemx.quantify.mvp import MonomerRecord
from tandemx.utils.kmers import canonical_kmer, iter_linear_canonical_kmers, reverse_complement


@pytest.mark.parametrize('k', [1, 3, 15, 31])
def test_native_sparse_weighted_counts_match_naive_word_oracle(k):
    core = pytest.importorskip('tandemx._rust_core')
    rng = random.Random(941+k)
    sequences = [''.join(rng.choices('ACGT', k=200)) for _ in range(8)]
    sequences += [s[:60]+'NNN'+reverse_complement(s[60:]) for s in sequences]
    sequences += ['N'*100, 'A', 'acgt'*20]
    words = sorted({canonical_kmer(s[i:i+k]) for s in sequences[:8] for i in range(0, 200-k+1, 7)})
    targets = {w: (i % 3, 1/(i+1)) for i, w in enumerate(words)}
    counter = core.WeightedKmerCounter(k, 3, [(w, i, v) for w, (i, v) in targets.items()])
    expected = []
    for sequence in sequences:
        values = defaultdict(float)
        for word in iter_linear_canonical_kmers(sequence, k):
            if word in targets:
                family, weight = targets[word]
                values[family] += weight
        expected.append(dict(values))
    for native, oracle in zip(counter.count_sequences(sequences), expected):
        assert dict(native) == pytest.approx(oracle)
    assert counter.count_sequences([]) == []
    with pytest.raises(ValueError, match='ACGTN'):
        counter.count_sequences(['ACGT', 'ACRY'])
    with pytest.raises(ValueError, match='ACGTN'):
        counter.count_sequences([''])


@pytest.mark.parametrize('backend', ['python', 'rust'])
@pytest.mark.parametrize('batch_bases', [1, 8_000_000])
def test_joint_estimator_points_match_word_aggregate_with_rc_n_and_short_reads(backend, batch_bases):
    if backend == 'rust':
        pytest.importorskip('tandemx._rust_core')
    rng = random.Random(1952)
    sequence = ''.join(rng.choices('ACGT', k=73))
    catalogue = [MonomerRecord('target', sequence), MonomerRecord('no_diagnostics', 'A'*41)]
    reads = [(sequence*8)[i % 70:i % 70+200+i] for i in range(40)]
    reads += [reverse_complement(s).lower() for s in reads] + ['AC', 'N'*100, sequence+'NN'+sequence]
    reference = estimate_multik(reads, catalogue, 10000, backend='python')
    joint = estimate_joint_multik(iter(reads), catalogue, 10000, backend=backend,
                                  batch_bases=batch_bases, minimum_effective_reads=2)
    assert (joint.read_count, joint.total_bases, joint.ambiguous_bases) == (
        reference.read_count, reference.total_bases, reference.ambiguous_bases)
    for a, b in zip(joint.per_k, reference.per_k):
        assert a.keys() == b.keys()
        for key in a:
            assert a[key] == pytest.approx(b[key]) if isinstance(b[key], float) else a[key] == b[key]
    for a, b in zip(joint.estimates, reference.estimates):
        assert a.estimated_copy_number == pytest.approx(b['extrapolated_copy_number'])
        assert a.fit_status == b['status']
    assert joint.estimates[0].log_sampling_variance > 0
    assert joint.estimates[1].sampling_interval_low is None
    assert 'ambiguous_base_windows' in joint.warning


def test_joint_collector_rejects_empty_invalid_and_duplicate_inputs():
    cat = [MonomerRecord('a', 'ACGTGCAT'*10)]
    for reads in [[], ['N'], ['ACGT', 'ACRY'], ['']]:
        with pytest.raises(ValueError):
            estimate_joint_multik(reads, cat, 100, backend='python')
    for catalogue in [[], cat*2, [MonomerRecord('a', 'ACRY')]]:
        with pytest.raises(ValueError):
            estimate_joint_multik(['ACGT'*10], catalogue, 100, backend='python')
    result = estimate_joint_multik(['ACGT'*20]*30,
        [MonomerRecord('a', 'ACGTGCAT'*10), MonomerRecord('b', 'ACGTGCAT'*10)], 100, backend='python')
    assert all(e.fit_status == 'missing_diagnostic_or_exposure' for e in result.estimates)
