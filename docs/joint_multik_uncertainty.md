# Experimental joint-read multi-k uncertainty

`tandemx.quantify.joint_multik.estimate_joint_multik` retains the fixed
15/21/27/31 point estimator from [the multi-k prototype](multik_quantification.md).
It adds an approximate sampling interval with each read as one independent
sampling cluster. No truth, error rate or planted coordinates enter inference.
This is a research API; the default public `quantify` output is unchanged.

For read i and k, define X_ik=max(0,L_i-k+1) and
Y_ifk=mean over exclusive diagnostic words j of count_ijk/multiplicity_fjk.
The point estimate at k is theta_fk=G*sum_i(Y_ifk)/sum_i(X_ik).
OLS fits log(theta_fk)=a_f+b_f*k. The intercept has fixed weights
w_k=1/K-mean(k)*(k-mean(k))/sum((k-mean(k))²).
Consequently the read's first-order influence on log(copy number) is

```text
g_if = sum_k w_k * (Y_ifk / sum_i Y_ifk - X_ik / sum_i X_ik)
v_f  = N/(N-1) * sum_i g_if²
sampling interval = exp(a_f ± 1.959963984540054 * sqrt(v_f))
```

All k axes share reads and their cross-products are retained. The k values are
not independent replicates. Streaming first and cross moments of X and Y give
the same variance as an explicit independent per-read influence calculation.
State scales with target words plus families*K². Sparse native read counters
use canonical two-bit words, merge reverse complements and reset at N. Sequence
batches are bounded at512 reads or8 Mb; one unusually long read can exceed8 Mb.
The Python collector remains a differential reference. No new dependency is used.

Zero support and missing diagnostic words give no extrapolation, not absence.
The fixed support guard requires at least20 effective reads on every k axis,
where effective reads=(sum Y)²/sum Y². Sparse support, variance below numerical
moment resolution and numerical failures give explicitly missing intervals.
The cancellation guard is1e-12 times the sum of absolute quadratic terms;
it prevents identical read contributions producing spurious tiny intervals.
Positive fitted slopes
retain their model inconsistency flag even if sampling variance is estimable.

Genome size and the supplied catalogue are treated as fixed. This interval does
not include extrapolation bias, unknown background sharing, catalogue discovery
error, genome-size error, biological pooling or correlated technical reads.
Reported base qualities are not inferred empirical error rates. N-overlapping
windows remain in exposure and are flagged. A sampling interval is not total
biological uncertainty, and a nominal95% label is not evidence of95% coverage.

## Reproducible calibration

First run the documented independent factorial generator, conditional baseline
and fixed multi-k replay. Then, inside `tandemx-dev`:

```bash
python -m benchmarks.scripts.evaluate_joint_multik \
  --previous /Volumes/T7/Codex/TandemX/results/factorial_multik_replay_v1_20260906 \
  --outdir /Volumes/T7/Codex/TandemX/results/factorial_joint_multik_v1_20260906
pytest -q tests/unit/test_joint_moments.py tests/unit/test_joint_multik.py tests/unit/test_joint_calibration_summary.py
```

The controller requires completed development inputs, verifies source/read/
catalogue/previous-result hashes, snapshots executable source, retains native
resources and checks point-estimate parity within floating-point tolerance.
It evaluates the27 pre-existing genome/coverage/error conditions; the three
independent genome seeds6301–6303 remain separate in summaries. All55 families
remain in denominators; missing intervals and available-but-missing-truth
intervals are separate outcomes. Conditions share reads/genomes, and families
are dependent, so pooled rows are not used as independent binomial trials.
Held-out seeds7301–7303 are untouched. Calibration results must be reviewed
before considering promotion into the public estimator.

## Minimal in-memory example

```python
from tandemx.quantify.joint_multik import estimate_joint_multik
from tandemx.quantify.mvp import MonomerRecord

monomer = 'ACGTCAGTGCATCGTAGCTAGGCTA'
result = estimate_joint_multik(
    [monomer * n for n in range(4, 44)],
    [MonomerRecord('example', monomer)], genome_size=10000,
    backend='python',
)
print(result.estimates[0])
```

This is a synthetic API demonstration, not a biological benchmark or a source
of paper results. Output fields and unavailable-value conventions are described
in [file_formats.md](file_formats.md#experimental-joint-read-multi-k-calibration).
