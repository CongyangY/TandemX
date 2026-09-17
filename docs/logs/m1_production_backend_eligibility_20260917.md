# M1 competitive mapping as a production quantify backend: eligibility gate

Decision: **not eligible for substitution or a public opt-in backend on current
evidence**. This is a backend decision, separate from the prior M1 Gate A
NO-GO for an absolute-abundance Methods novelty claim. The frozen diagnostic
k-mer implementation and `copy_number.tsv` contract remain unchanged.

The executable audit is
`benchmarks/m1_shared_signature/check_backend_eligibility.py`; its checked
receipt is `docs/evidence/m1_backend_eligibility_20260917.json`. Run with:

```bash
conda run -n tandemx-dev python -m benchmarks.m1_shared_signature.check_backend_eligibility --output /tmp/m1_backend_eligibility.json
```

It verified both archived source-hash maps and equal input hashes across all
12 M1 development cases. Ordinary mapping positive-family MARE is 0.0034,
0.0144, 0.3796 and 1.0 across low error, high error, family-specific error
and identical-family cases. The revised occupancy prototype gives 0.0034,
0.0164, 0.3658 and 1.0. Its advantage in the family-specific-error condition
does not dominate the other identifiable conditions. These counts are from
three toy seeds of unit-sized, substitution-only reads, not independent
biological materials. The identical-family result is an individual-family
refusal, not a successful individual estimate.

The actual ordinary mapping function requires every read to equal every
catalogue unit in length. The audit presents a two-unit read to this unchanged
function and confirms that it raises `ValueError`. It also confirms that
identical family units yield a tie with no individual attribution. Thus the
existing mapping baseline cannot ingest the full-length, mixed-content reads
accepted by public `quantify`. The occupancy prototype is a separate fixed
80-bp-window, eight-family-capped research implementation; its optional edlib
speed and ~35.7 MB toy-process RSS are not paired production measurements.

The endpoint also differs. Public `quantify` counts diagnostic k-mers over
whole reads, applies an optional survival correction, chooses a named
haploid-depth normalization, and emits estimated copies and bp, uncertainty
fields and warnings. The M1 mapping result is assigned unit-sized read counts;
dividing assigned read bp by a depth would be a new, uncalibrated model. The
M1 receipts contain no independent full-read physical-bp calibration, no
assembly-deficit validation for this backend, and no paired public-backend
runtime/RSS comparison. Neither the mapping baseline nor the occupancy
prototype can be placed behind the existing `--kmer-backend` option, which
selects an exact counting implementation rather than a scientific estimator.

## Eligibility for a future, separately authorized backend

Before exposing any mapping-based backend, freeze its complete-read assignment
and abstention rules using a development split only. Require an input-only
identifiability decision for shared/identical units, explicit unassigned and
unknown fractions, and a depth/length calibration that produces the same
physical units as `copy_number.tsv` without imputing ambiguous bases to a
family. Then compare the frozen candidate with the current production backend
and ordinary mapping on independent donor/lineage data, including mixed-family
reads, background, indels, divergent members, zero-abundance decoys and
shared-sequence families. Report all-family error with refusals included,
wrong-family assignment, false decoy abundance, continuous copy-bp bias and
absolute error, uncertainty coverage, and assembly-deficit magnitude and
binary classification separately. Measure equal-input wall time and peak RSS
at realistic catalogue/read scale. Do not select parameters on a final holdout
or use the already consumed M1 toy cases as independent validation.

This gate is falsifiable: a frozen complete-read backend with independent
matched-endpoint validation and resource receipts can clear it. The current
research implementations fail before that final validation gate, so no
unconsumed holdout or T7 large-read run was opened here.
