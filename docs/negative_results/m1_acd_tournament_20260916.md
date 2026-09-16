# M1 A-F bounded development tournament: no estimator beats ordinary mapping

Status: **one bounded development round complete; stop new M1 estimator search
under this protocol**. The 12 cases are the same three toy seeds across four
conditions, with the exact same catalogue/read SHA-256 and planted truth as
the earlier B/E trials. They are not 12 independent biological materials.
`benchmarks/m1_shared_signature/evidence_acd_20260916/results.json` contains
all cases, source and input hashes, point estimates, confusion, rejection,
timing, memory, and the prior read-bootstrap reference. The script refuses to
run if B/E source files or prior input hashes differ.

The ordinary baseline is unique-best cyclic-Hamming mapping with a shared
18-mismatch gate. Candidate **A** uses catalogue-exclusive circular 5-mers
with a fixed 8% error-survival correction. **B** is the prior nonnegative
shared-signature regression. **C** fits an error-aware 5-mer Poisson count
model with nonnegative abundances summing to accepted reads. **D** fits a
per-read composite 5-mer mixture by EM and assigns posterior ≥0.95. **E** is
the prior alignment-likelihood per-read EM plus unknown rejection. **F** is
the prior 40-read-bootstrap check for B/E; A/C/D have no interval calibration.
The 8% error assumption, posterior cutoffs and model forms were fixed before
this A/C/D execution. No true source label is used by an estimator.
An input-only singular-value ratio below 0.05 forces group-level refusal in
A/C/D; this safeguard is a development setting, not calibrated certainty.

### All-case primary endpoint

The table reports mean positive-family MARE for every seed; individual
refusal incurs 100% relative error. The planted-zero decoy is scored
separately, without division by zero.

| Method | 3% error | 12% error | `f2` extra 20% error | Identical `f1`/`f2` |
| --- | ---: | ---: | ---: | ---: |
| Ordinary mapping | **0.0034** | **0.0144** | **0.3796** | 1.0000; 500 ties |
| A exclusive signature | 0.3645 | 0.0652 | 0.4634 | 1.0000; group only |
| B shared-signature regression | 0.0680 | 0.2359 | 0.4712 | 1.0000; group only |
| C Poisson count | 0.0209 | 0.0382 | 0.4807 | 1.0000; group only |
| D 5-mer EM | 0.0049 | 0.0207 | 0.3966 | 1.0000; group only |
| E per-read alignment probability | 0.0043 | 0.0686 | 0.4046 | 1.0000; group only |

No A-E candidate beats ordinary mapping in any of the three identifiable
conditions. Identical units permit a combined group estimate, not individual
family abundance; the family-level refusal penalty is retained for all six
methods. A group-level success cannot be substituted for that endpoint.

### Error, false assignment, rejection, and cost

Across the 3% / 12% / family-error-shift conditions respectively:

| Method | Mean positive MAE, reads | Mean positive bias, reads | Mean zero-decoy assignment, reads | Mean rejected reads / 600 | Mean fit seconds / case |
| --- | --- | --- | --- | --- | --- |
| Ordinary mapping | 0.83 / 3.83 / 68.83 | -0.50 / -3.83 / -68.50 | 0 / 0 / 0 | 101.00 / 107.67 / 237.00 | 0.022 / 0.023 / 0.024 |
| A exclusive signature | 89.77 / 17.93 / 94.95 | +89.77 / -14.93 / -42.18 | 12.15 / 39.36 / 18.74 | 100.00 / 101.00 / 231.67 | 0.148 / 0.146 / 0.142 |
| B shared-signature | 15.89 / 54.60 / 93.11 | -15.89 / -54.60 / -93.11 | 31.78 / 108.21 / 54.56 | 100.00 / 101.00 / 231.67 | 0.011 / 0.011 / 0.008 |
| C Poisson | 4.82 / 8.90 / 92.28 | 0.00 / -8.90 / -65.95 | 0.00 / 16.80 / 0.23 | 100.00 / 101.00 / 231.67 | 0.162 / 0.153 / 0.144 |
| D 5-mer EM | 1.17 / 5.83 / 71.83 | -0.17 / -5.83 / -69.83 | 0 / 0.33 / 0.33 | 100.33 / 111.33 / 239.33 | 0.152 / 0.148 / 0.162 |
| E per-read probability | 1.00 / 13.67 / 73.33 | -0.67 / -13.67 / -73.00 | 0 / 0 / 0 | 101.33 / 127.33 / 246.00 | 0.024 / 0.025 / 0.025 |

The 100 random negative reads per case failed the common gate; none entered
A-D or received an E family assignment. D made 50 wrong-family calls among
the nine identifiable cases, despite near-zero aggregate decoy calls. A/C
have only aggregate feature counts, so a per-read wrong-family attribution
rate is **N/A**, not zero. A's estimated totals may exceed accepted reads;
its strong overcount at low error shows the exclusive-word survival correction
is misspecified. C assigns ~16.8 reads/case to a truly absent decoy at 12%
error. High error and family-specific error produce large common-gate losses;
more aggressive refusal does not recover true abundance.

F read-bootstrap coverage was B 2/18 and E 14/18 across the identifiable
positive-family cells. These 18 family/seed cells are correlated development
observations, not 18 independent donors. A/C/D have no intervals and cannot
claim calibration. Wall times are one-process development measurements on this
host, include A/C/D design construction, and exclude toy generation; RSS is
recorded as macOS raw `ru_maxrss` bytes. They are not production or large
genome scaling results.

Decision: preserve all A-F results and keep ordinary competitive mapping as
the strong fixed-catalogue development baseline. No candidate meets the
preregistered direction of improvement in this small common-input screen.
Stop this round of M1 estimator search; do not change the frozen production
quantifier or use a final held-out set to select new thresholds. The toy setup
still lacks indels, real full-genome background, donor variation, GC bias,
long mixed-family reads, and physical copy-number truth. These gaps prevent
both a universal rejection of future abundance models and a Methods accuracy
claim for the current ones.
