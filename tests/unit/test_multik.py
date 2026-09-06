import math
import random
from collections import Counter

import pytest

from tandemx.quantify.multik import estimate_multik, fit_attenuation
from tandemx.quantify.mvp import MonomerRecord


def test_log_survival_fit_known_parameters_and_leave_one_out():
    ks = (31,15,27,21)
    values = [83 * math.exp(-.024*k) for k in ks]
    fit = fit_attenuation(ks,values)
    assert fit.extrapolated_copy_number == pytest.approx(83)
    assert fit.log_slope_per_base == pytest.approx(-.024)
    assert fit.effective_word_loss_probability == pytest.approx(1-math.exp(-.024))
    assert fit.max_absolute_log_residual < 1e-14
    assert fit.leave_one_k_out_log_range < 1e-14
    assert fit_attenuation((15,21,31),[2,3,4]).status == 'positive_slope_inconsistent_with_loss'
    assert fit_attenuation((15,21,31),[0,0,0]).extrapolated_copy_number is None
    assert fit_attenuation((15,21,31),[2,None,3]).status == 'missing_diagnostic_or_exposure'
    with pytest.raises(ValueError):
        fit_attenuation((15,21,31),[2,float('nan'),3])
    with pytest.raises(ValueError):
        fit_attenuation((15,21,21),[2,3,4])
    with pytest.raises(ValueError):
        fit_attenuation((15,21,31),[2,3])


def _canonical(word):
    return min(word,word.translate(str.maketrans('ACGT','TGCA'))[::-1])


def test_multik_counting_matches_independent_naive_oracle_and_backends():
    rng=random.Random(20032)
    founder=''.join(''.join(rng.sample('ACGT',4)) for _ in range(16))
    # Double founder exercises canonical-word multiplicities of two.
    motif=founder*2
    assert all(len(set((motif+motif)[i:i+k]))==4 for k in (15,21,27,31) for i in range(len(motif)))
    reads=[founder*8,'N'*6+founder*4,founder[:8],(founder*3).translate(str.maketrans('ACGT','TGCA'))[::-1]]
    catalogue=[MonomerRecord('F',motif)]
    python=estimate_multik(iter(reads),catalogue,2000,backend='python',batch_bases=500)
    rust=estimate_multik(iter(reads),catalogue,2000,backend='rust',batch_bases=300)
    assert python == rust
    assert python.read_count==4 and python.total_bases==sum(map(len,reads)) and python.ambiguous_bases==6
    for row in python.per_k:
        k=row['k']
        target=Counter(_canonical((motif+motif)[i:i+k]) for i in range(len(motif)))
        counts=Counter(_canonical(read[i:i+k]) for read in reads for i in range(len(read)-k+1)
                       if 'N' not in read[i:i+k])
        exposure=sum(max(0,len(read)-k+1) for read in reads)
        expected=sum(counts[w]/n for w,n in target.items())/len(target)
        assert row['mean_corrected_count']==pytest.approx(expected)
        assert row['uncorrected_copy_number']==pytest.approx(2000*expected/exposure)
    assert 'ambiguous_base_windows' in python.estimates[0]['warning']
    assert estimate_multik((s.lower() for s in reads),catalogue,2000,backend='python')==python


def test_multik_ambiguity_empty_low_complexity_and_invalid_inputs():
    a='ACGGTCTATGCTCGAGGATCGAGTCATGCAAGCTGACTACGTAGCATGCTAGGCTATG'
    duplicate=[MonomerRecord('F',a),MonomerRecord('G',a)]
    result=estimate_multik([a*2],duplicate,1000,backend='python')
    assert all(r['extrapolated_copy_number'] is None for r in result.estimates)
    result=estimate_multik(['N'*100],[MonomerRecord('F',a)],1000,backend='python')
    assert result.estimates[0]['status']=='zero_support_at_one_or_more_k_not_absence'
    low=estimate_multik(['A'*100],[MonomerRecord('F','A'*50)],1000,backend='rust')
    assert low.estimates[0]['status']=='missing_diagnostic_or_exposure'
    for reads in ([],[''],['ACGTZ'],['A'*10]):
        with pytest.raises(ValueError):
            estimate_multik(reads,[MonomerRecord('F',a)],1000,backend='python')
    with pytest.raises(ValueError):
        estimate_multik([a],duplicate[:1],float('nan'),backend='python')
