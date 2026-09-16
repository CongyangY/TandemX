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

## Second candidate: per-read probability and unknown rejection

`run_probabilistic.py` replays those exact 12 cases, checks both catalogue and
read SHA-256 values against the first tournament, and compares all three
methods without changing the first evidence directory:

```bash
conda run -n tandemx-dev python -m benchmarks.m1_shared_signature.run_probabilistic \
  --outdir /tmp/tandemx-m1-prob
conda run -n tandemx-dev pytest -q tests/unit/test_m1_probabilistic.py
```

The candidate computes the best circular-Hamming distance per read and
catalogue family, estimates a bounded per-read substitution rate from the
closest alignment, then fits family and uniform-DNA unknown mixing weights
using EM. A read is assigned only if its best posterior is at least 0.95 and
the common 18-mismatch gate passes. Others are marked `gate`, `unknown`, or
`ambiguous`. Exact cyclic-equivalent units get one group estimate; no
individual count or interval is invented. The settings were fixed before the
first development run; they are not validated calibration constants.

The runner writes source, prior-evidence and input hashes; full source-to-call
confusion; positive-family MARE, MAE, bias and refusal-penalized error; decoy
assignment; 40-read-bootstrap percentile intervals and coverage; and fit and
bootstrap time. It uses source truth only in scoring. The read-derived error
estimate and uniform unknown model can be misspecified, especially for indels,
platform errors, near-family decoys, or real full-genome background. The
descriptive read bootstrap does not validate 95% biological coverage.

## Bounded A/C/D completion

`run_ablations.py` adds three fixed development candidates and checks that
the catalogue/read hashes and B/E point estimates still match both earlier
evidence files:

```bash
conda run -n tandemx-dev python -m benchmarks.m1_shared_signature.run_ablations \
  --outdir /tmp/tandemx-m1-acd
conda run -n tandemx-dev pytest -q tests/unit/test_m1_ablations.py
```

**A** weights only catalogue-exclusive circular 5-mers, correcting their
expected count by a fixed 8% substitution survival factor. It has no
per-read attribution beyond the common gate. **C** fits nonnegative family
counts under an error-aware Poisson 5-mer count model, constrained to sum to
accepted reads. **D** applies EM to each accepted read's composite 5-mer
likelihood and assigns only posterior ≥0.95; overlapping 5-mer windows are
not independent, so this is a pseudo-likelihood. A/C/D use the same
predeclared 8% substitution assumption and the same known catalogue; no
truth labels are input to inference. D can add ambiguous refusals. Exact
equal-signature families are group-resolved only. A/C/D also use a fixed
input-design singular-value ratio floor of 0.05; a more degenerate design
returns only an aggregate group count, with no individual estimates.

The runner reexecutes ordinary mapping, B, and E in the same process, checks
their point estimates against the prior records, and reports all six methods
plus the prior F read-bootstrap outcomes. A/C/D are point-estimate screens;
they have no 95% interval calibration result. The prior F read bootstraps for
B/E do not establish biological coverage. This closes one bounded M1 search
round rather than opening a parameter or held-out search.
