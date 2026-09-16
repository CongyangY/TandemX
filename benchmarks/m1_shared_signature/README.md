# M1 fixed-catalogue development tournament

This is a bounded research prototype. It does not change `tandemx quantify`,
the frozen estimator, its historical outputs, or the manuscript. Run it in
`tandemx-dev` on an internal-disk output directory:

```bash
conda run -n tandemx-dev python -m benchmarks.m1_shared_signature.run --outdir /tmp/tandemx-m1-dev
conda run -n tandemx-dev pytest -q tests/unit/test_m1_shared_signature.py
```

The deterministic generator creates a known catalogue of three 80-bp units:
two present families (320 and 180 unit-sized reads) and one zero-abundance
decoy, plus 100 random-sequence negative reads. Related families differ at
three positions; the unidentifiable case uses identical units. Unit-sized
reads are randomly rotated and independently substituted. Three seeds and
four conditions are *development* trials, not independent organisms or an
unseen holdout. The family-error-shift condition raises the substitution rate
for family `f2` by 0.20; it is not a coverage-bias simulation.

Both methods receive the exact same catalogue, read list, 80-bp truth, and
cyclic-Hamming alignment gate (at most 18 mismatches). The ordinary mapping
baseline assigns a read to its unique best aligned family and leaves ties
unassigned. The shared-signature candidate fits circular 5-mer counts from
accepted reads to catalogue signature columns, with nonnegative weights
constrained to sum to the accepted-read count. It cannot secretly use truth
labels. Exactly equal columns and a smallest/largest singular-value ratio
below 0.05 trigger an input-only individual-family refusal; only the combined
group count is returned. The 0.05 safeguard is a development choice, not a
validated biological identifiability threshold. Forty read resamples give
exploratory percentile intervals, conditional on this catalogue and gate.

`results.json` records source and input SHA-256 values, every case, mapping and
regression estimates, all rejected/tied/unknown-source counts, zero-decoy
assignment, individual refusal, per-case wall time, and process peak RSS.
`all_positive_mare_refusal_penalized` counts an individual refusal as 100%
relative error, keeping all planted positive families in the denominator.
`identifiable_positive_mare` excludes refused families and must not be used
alone to select a method. The interval-coverage indicators are descriptive;
read bootstrap samples do not create independent biological replicates.

Current boundaries: exact unit-sized reads only, substitutions without indels,
forward-strand rotations, three equal-length catalogue units, no full-genome
background or GC bias, no unknown-catalogue discovery, no physical copy-number
calibration, and no assembly deficit. The prototype uses NumPy already
available in the development environment; no production dependency is added.
