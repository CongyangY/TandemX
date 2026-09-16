# M1 per-read probability: development loss against ordinary mapping

Status: **second M1 candidate rejected on the fixed development tournament**.
This study reused exactly the 12 catalogue/read/truth hashes from the prior
shared-signature trial. Both experimental candidates and the ordinary
competitive mapping baseline saw the same input catalogue and reads. Truth
labels entered only the score. Full cases, timing, confusion, bootstrap values
and source hashes are in
`benchmarks/m1_shared_signature/evidence_probabilistic_20260916/results.json`.

The new model estimates a per-read substitution probability from its nearest
cyclic alignment, fits family plus uniform unknown weights with EM, and
requires posterior ≥0.95 plus the same 18-mismatch gate for assignment. Exact
cyclic-equivalent families produce a group count and individual refusal. These
are development settings, not calibrated confidence guarantees.

| Condition, 3 seeds | Ordinary mapping positive MARE | Shared 5-mer fit positive MARE | Per-read probability positive MARE | Per-read probability mean positive MAE (reads) | Per-read probability mean bias (reads) | Mean rejected / 600 | Zero-decoy assignment |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Related, 3% substitutions | 0.0034 | 0.0680 | 0.0043 | 1.00 | -0.67 | 101.33 | 0 |
| Related, 12% substitutions | 0.0144 | 0.2359 | 0.0686 | 13.67 | -13.67 | 127.33 | 0 |
| `f2` extra 20% substitutions | 0.3796 | 0.4712 | 0.4046 | 73.33 | -73.00 | 246.00 | 0 |
| Identical `f1`/`f2` | 1.0000 | 1.0000 | 1.0000 | 250.00 (individual refusal penalty) | -250.00 | 100.00 | 0 |

The new model prevented false assignment to the zero-abundance decoy and
rejected all 100 random negative reads in each case. Ordinary mapping also
assigned zero reads to that decoy, so this is not an advantage over the strong
baseline. Wrong-family attribution across the nine identifiable cases totaled
7 reads; unknown-source attribution was 0. The identical-unit family-group
total was 500/500 in all three seeds, but individual family counts are not
identifiable and remain refused.

Read-bootstrap intervals covered 6/6 positive truth values at low error,
5/6 at high error, and 3/6 under family-specific error: 14/18 total. This
small conditional development count is insufficient for nominal 95% interval
calibration, and the poor high-error performance makes calibration unlikely
under this model. Mean single-fit wall time was about 0.02 s/case; 40 bootstrap
fits took about 0.7--0.8 s/case on this host. One-process peak `ru_maxrss`
and per-case timings are in the JSON, with no cross-tool whole-pipeline
resource claim.

The most direct failure is extra refusal of real family reads: at 12% error,
the model rejected about 27 additional reads/case beyond random negatives;
under family-specific error it rejected about 146 additional true reads/case.
The shared 18-mismatch gate alone already removes many high-error `f2` reads,
and the posterior cutoff removes more. The model does not beat ordinary mapping
on the refusal-penalized positive-family endpoint in any identifiable
condition. There are only three random seeds of one synthetic design, no
independent donors, indels, realistic background, or held-out genomes.

Decision: retain this second negative result and stop M1 promotion. Do not
modify the frozen production estimator, tune to the final holdout, or claim a
calibrated abundance/assembly-deficit method. A future M1 revival would need a
new preregistered error/unknown model and independent validation, not an
after-the-fact threshold adjustment on these cases.
