# M2 simple baseline, frozen before B1 v3 predictions

Source commit `35be831` established a paired-flank array-length baseline. A
subsequent sign correction aligned its optional delta field with the preexisting
B1 ledger convention (`edited assembly - original/read span`). No threshold or
event decision changed. The v3 input-only bundle and v3 truth had not been
generated or opened at that point.

On B1 v2's 13 cases from one engineered material, the baseline produced 13
decisions: TP 9, FN 2, FP 0, TN 2; sensitivity 9/11, intact FPR 0/2. The two
missed positives are equal-length structure changes, which span comparison
cannot detect. Its signed bp difference is exact on these error-free synthetic
full-span reads (MAE 0 bp). This says nothing about real read errors, partial
spans, breakpoint placement, family or HOR annotation, or independent biological
array truth. The three read IDs per case repeat one synthetic sequence.

The two M2 routes each issued one abstention and therefore have null all-case
primary metrics under the frozen v2 rule. Comparing baseline's 9/11 with either
route's eligible-only 10/10 as if they were the same denominator would be
invalid. The baseline remains an explicit stop-loss comparator for v3 and future
held-out work. It may outperform M2 on length-only arrays.

Reproduce with `benchmarks/m2_routes/baselines/README.md` and the v2 scorer.
Predictions and scorer receipt are under `benchmarks/m2_routes/baselines/evidence/`
and `benchmarks/controlled_collapse/v2/flank_length_score/`.

Before opening B1 v3 truth, the frozen baseline was run on its 39 public
input-only rows. It returned 18 `DISCORDANT` and 21 `SUPPORTED`, no abstentions;
the raw predictions are archived as `v3_predictions.jsonl`. These are decision
counts, not performance metrics. v3 still contains exact synthetic full-length
reads and engineered flanks, so even a high scored result cannot validate
native molecule anchoring.
