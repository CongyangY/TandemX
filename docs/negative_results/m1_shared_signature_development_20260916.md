# M1 shared-signature constrained regression: development loss

Status: **development negative result; stop this candidate**. This is a
fixed-catalogue, three-seed toy tournament, not an external or biological
validation. The executable protocol and full per-case record are in
`benchmarks/m1_shared_signature/README.md` and
`benchmarks/m1_shared_signature/evidence_20260916/results.json`.

The same generated catalogue, reads, cyclic-Hamming gate and source truth were
used for ordinary competitive mapping and the shared 5-mer constrained fit.
All positive families remain in the `all_positive_mare_refusal_penalized`
denominator. A zero-abundance decoy and 100 random negative reads are recorded
separately. The exact-equal family case requires individual refusal and a
combined group count.

| Condition (3 seeds each) | Mapping positive MARE | Shared fit positive MARE | Shared fit mean false assignment to zero decoy (reads) | Mean rejected reads / 600 | Shared fit positive bootstrap coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| Related, 3% substitutions | 0.0034 | 0.0680 | 31.78 | 100.00 | 2 / 6 |
| Related, 12% substitutions | 0.0144 | 0.2359 | 108.21 | 101.00 | 0 / 6 |
| `f2` extra 20% substitutions | 0.3796 | 0.4712 | 54.56 | 231.67 | 0 / 6 |
| Identical `f1`/`f2` | 1.0000 (ties left unassigned) | individual estimate refused; group reported | 36.36 | 100.00 | N/A |

The shared fit loses to the ordinary-mapping baseline in every identifiable
condition. It gives a positive estimate to the planted-zero decoy in all
12 cases. Its read-bootstrap intervals cover only 2 of 18 identifiable
positive truth values in this small development set. High family-specific
error causes many true `f2` reads to fail the common gate, so both methods
underestimate it; shared regression does not repair that failure. In the
identical-unit case the shared fit appropriately refuses individual family
values, but its group total is also biased downward by decoy allocation.

Mechanistic interpretation is limited: substitution noise moves observed
5-mers away from error-free catalogue signatures, so a constrained mixture
can explain noisy counts by misallocating reads to a related decoy. The
prototype has no error-aware exposure model. This explains a plausible failure
mode; it is not proof that every shared-feature model must fail.

Decision: retain the code and failed evidence for audit, do not promote the
model into production, do not tune it on a held-out set, and do not claim a
calibrated abundance or missing-bp Methods advance from this trial. A future
distinct model would require a preregistered error/exposure model, stronger
background and indel tests, locked baseline and power plan, and an independent
organism/family holdout. The current study cannot set or pass the G1 gate.
