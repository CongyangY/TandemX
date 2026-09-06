# Experimental read-cluster quantification

This Python reference in `tandemx/quantify/read_moments.py` is a development
experiment, not the default `tandemx quantify`. Existing TSVs and median/k-mer
spread estimates remain reproducible. No new statistical principle or calibrated
biological confidence is claimed merely from implementing a ratio standard error.

For D unique diagnostic canonical k-mers, m_j is multiplicity in the circular
monomer and N_ij is its count in read i:

```
Y_i = sum_j(N_ij / m_j) / D
X_i = max(0, read_length_i - k + 1)
R = sum(Y_i) / sum(X_i)
estimated_copy_number = supplied_genome_size * R
SE = supplied_genome_size / sum(X_i)
     * sqrt(n / (n - 1) * sum_i((Y_i - R * X_i)^2))
```

All overlapping words within a read contribute to the same read-level residual.
The exposure corrects finite read-end opportunity loss under uniform sampling;
it does not correct disruption of circular monomer k-mers at array boundaries.
The point estimate is a weighted mean, not the old median estimator.

Stored moments are n, sum(X), sum(X²), and each family's sum(Y), sum(Y²), sum(XY),
and positive-read count. Zero contributions are implicit; there is no full
read-by-family matrix. Memory scales with targets/families and current-read
contributions. Python rolling scans support k=1..31; the replay uses the existing
catalogue's shared-k-mer and low-complexity filters. Background uniqueness remains
unverified. Genome size is treated as fixed. N-containing opportunities remain
in exposure with a warning rather than silently inflating observed depth.

Normal intervals use estimate ±1.959963984540054 SE, lower bound clipped at zero.
At least 20 effective positive reads, `(sum Y)^2/sum(Y²)`, are required. This is
a disclosed development support guard, not a validated coverage guarantee.
No targets produce a missing estimate; no observed hits produce estimate zero
with a missing interval and `no_observed_support_not_absence`. Sparse support
and zero observed variance also yield missing intervals. Report availability
alongside coverage among available intervals; never replace missingness by [0,0].

Independent uniform read sampling is an assumption. Sequencing errors, biological
unit divergence, library/GC/length selection, donor mixtures, organelles, genome
size error and duplicate molecules are not corrected. Neither reported FASTQ
quality nor planted simulation error rates enter the estimator. It therefore
cannot solve the baseline mutation bias by uncertainty modeling alone. Further
real/semisynthetic calibration and resource validation precede CLI integration.

## Development replay and toy check

```bash
conda run --no-capture-output -n tandemx-dev python -m benchmarks.scripts.evaluate_read_moments \
  --baseline /path/to/completed/abundance_experiment \
  --outdir /path/to/new/read_cluster_replay
conda run --no-capture-output -n tandemx-dev pytest -q tests/unit/test_read_moments.py
```

The replay consumes only development inputs from the independent abundance
experiment. The estimator sees reads, supplied catalogue and genome size; only
the evaluator sees planted counts, read-coordinate oracle and error labels.
The worker executes from a source snapshot with input hashes and receipts.
Raw estimates, errors, missing-interval statuses and conditional summaries are
retained. Its resource observations are diagnostics, not external-tool rankings.

Tests use an independent naive sequence scan, direct residual calculation,
reverse strands/N/short/empty inputs, invalid targets, and interval coverage
under an independent Bernoulli model. An independently generated small genome
checks replay and changed-input/held-out rejection. That mathematical model
test does not establish biological calibration on real satellites.
