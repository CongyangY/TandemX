# FPR scope reconciliation: A3 development versus historical classifier validation

## Decision

The A3 synthetic FPR and the historical 0.0164609053 false-positive rate have
different units of analysis, truth construction, denominators, decision inputs,
and simulation structures.  They must remain separately named and must not be
compared as absolute false-positive rates.  In particular, the A3 result must
not be written as a result against the historical 0.016 rate, nor as a generic
biological false-positive rate.

This is a read-only reconciliation.  It did not rerun a simulation, inspect a
new hold-out, or alter a threshold.

## A3 development synthetic FPR

The final A3 assessment is
`benchmarks/evidence/phase_and_scaling_gate_20260910/accuracy_r2_assessment.json`.
For `A3_round_two` it records 120 family-condition rows, sensitivity 0.95,
precision 0.890625 and FPR 0.1166666667.  The computation is explicit in
`benchmarks/scripts/evaluate_phase_round_two.py`:

```python
selected = [row for row in rows if row["model"] == model]
tp = sum(row["binary_call"] and row["binary_truth"] for row in selected)
fp = sum(row["binary_call"] and not row["binary_truth"] for row in selected)
pos = sum(row["binary_truth"] for row in selected)
FPR = fp / (len(selected) - pos)
```

Thus the denominator is 60 synthetic **F2 family-condition** rows and the
numerator is 7 calls on those rows (`7 / 60 = 0.1166666667`).  It is not a
base-level false-positive fraction, a negative-read rate, a no-array control
rate, or a rate over independent organisms.  The 120 rows are 2 development
seeds × 10 conditions × 3 coverages × 2 families.  As recorded in
`accuracy_r1_protocol.json`, the only seeds are `dev-phase-gate-01` and
`dev-phase-gate-02`; the evidence itself says holdout was not generated or
inspected.

`generate_case` in `benchmarks/challenge/phase_gate_fixture.py` creates a
22,000-bp source with an F1 array (six copies, or two in `abundance_low`) and
an F2 array (12 copies).  It samples 1,000-bp reads at 2/5/10×.  Conditions
include substitution and indel cases, partial arrays, a 300-bp interruption,
F1 80-bp fragment background, close families, heterogeneous units, and low F1
abundance.  The fragment background is explicitly not complete-family truth;
it is not a collection of absent-family genomes.

The classification convention is also fixture-specific.  The protocol states
that source assembly amount is `0.3 * truth` for F1 and `1 * truth` for F2;
the runner estimates source bp by dividing sampled-repeat bp by observed
coverage, then calls when `A / R < 0.6`.  `binary_truth` is simply `family ==
"F1"`.  The assessed A3 configuration has maximum seed gap 45, maximum drift
20, phase fraction 0.7, at least two units, and the second round only bridges
short gaps in already accepted chains.  These settings define the A3 interval
eligibility experiment; they are not the historical classifier's eligibility
or confidence rules.

The round-two assessment carries a development gate named
`binary_FPR_at_most_002`, whose Boolean result is false.  That historical Boolean is retained only as an executed-record artifact.
Its use as an acceptance criterion is withdrawn after this metric audit; it
must not be used to say that A3 failed an FPR<=0.02 gate. The stop decision
uses within-development relative evidence instead.

## Historical 0.0164609053 classifier FPR

The historical value is in
`paper/evidence/abundance_classifier_depth_gated_validation_v1/heldout/validation.json`.
The selected depth-gated classifier has TP/FN/FP/TN = 953/505/16/956 across
2,430 paired family conditions, so its FPR is `16 / (16 + 956) =
0.01646090534979424`.  The single-k21 baseline has the same 16/956 rate.  The
three held-out seeds are 5801--5803; the validation receipt labels the result
complete and held-out.

Here a row is a known-catalogue family in the conditional assembly simulation,
not an F1/F2 role.  The archived held-out `run_config.json` specifies three
monomer periods (61/171/421), copy counts (20/80/200), three unit-substitution
levels (1/3/5%), one or three array fragments, 1/5/20× read coverage, three
read substitution levels, and five assembly-retention levels (1/.75/.5/.25/0).
The truth label is `truth_assembly_read_ratio < 0.6`; a call is generated from
the assembly/read estimate under the frozen depth-gated rule.  The evaluator
in `benchmarks/abundance/evaluate_multik_collapse.py` defines FPR as
`FP / (FP + TN)`, and the figure legend records the same formula.

Unlike A3, this endpoint includes assembly localization, read abundance and a
depth-gated classifier: below estimated depth 2 it uses single k=21 at 0.6;
otherwise it uses the preselected alpha-0.5 blend at 0.5.  All 2,430 rows were
available.  There are 1,458 positive rows and 972 negative rows.  The result is
still an IID-substitution, known-catalogue simulation; the archive explicitly
does not call it biological collapse validation or a universal performance
claim.

## Why a numerical comparison is invalid

| Property | A3 synthetic FPR | Historical classifier FPR |
| --- | --- | --- |
| Split | Two development source labels; no hold-out | Three separately reserved held-out seeds |
| Unit counted | F2 role in a sampled-read family-condition | Simulated assembly/read family-condition |
| Negative definition | F2 has source assembly amount equal to source repeat truth | Known assembly/read ratio is at least 0.6 |
| Denominator | 60 F2 rows | 972 classifier-negative rows |
| Call evidence | A3 eligible sampled-read interval, scaled by observed coverage | Frozen k-mer abundance plus assembly localization under depth gate |
| Simulation | 22-kb two-family source with designed periodic-context decoys | Factorial three-family conditional assembly/read matrix |
| Threshold use | Fixed `A/R < 0.6` with F1 at 0.3 and F2 at 1.0 | Truth ratio threshold 0.6; decision threshold is 0.6 or 0.5 by frozen depth rule |

Both rates are internally well-defined for their own records.  They do not
estimate the same probability, and similar notation does not make them a
fair method comparison.  The useful A3 findings are relative within its
protocol: round two improved indel interval recall over strict A2
(0.916558 versus 0.107885) while its development MARE remained above A1
(0.074693 versus 0.056828).  Those are development observations only.

## Reporting language

Use “A3 synthetic F2-row false-call rate, 7/60 under the phase-gate fixture”
when the value must be named.  Use “held-out conditional assembly-classifier
FPR, 16/972” for the historical result.  Keep the phase-gate result separate
from the frozen classifier evidence and do not turn either into a biological
false-positive-rate claim.
