# B1 v2: injected structure development cases

This is a frozen **component development** bundle. The only biological sequence
source is an archived 785 bp YSD56 candidate; four 24 bp tiles from it are
repeated into an engineered ABCD array. The HOR, copy count, edit coordinates,
flanks, and reads are engineered, not observed native array structure. The
three `raw_read_sequences` per case are identical exact synthetic paths under
different IDs. They provide zero independent molecules or sampled error
profiles; read-support thresholds and score calibration cannot be inferred.
There is one source material and no validation or held-out set.

`protocol.json` fixes the 13 case operations, source SHA-256, source tile
offsets, random seeds, and a 0.5 uncalibrated event-decision threshold.
`development_bundle/inputs.jsonl` is the input-only interface shared with M2
routes. `truth.jsonl` contains the causal edit ledger and must stay separate
from route training/inference. `receipt.json` hashes the protocol, generator,
inputs, truth, and each assembly FASTA. Case IDs are opaque labels, not secure
blinding: the protocol and truth are committed for audit. Coordinates are
0-based, half-open bp intervals in the full engineered molecule.

There are **13 causal cases but only 10 distinct edited assembly sequences**.
`copy_loss`, `HOR_copy_deletion`, and `boundary_truncation` produce exactly the
same sequence because they remove different copies of an identical HOR; their
exact deletion locus and event subtype are not identifiable from these inputs.
The two intact controls also have identical sequences and are not independent
replicates. The scorer keeps exact causal truth while marking the three mixed
event-type rows `not_identifiable_from_available_sequences` and their subtype
matches `null`. An exact event subtype error is not attributed to a route in
those rows. The 13-to-10 equivalence classes are exposed in per-case outputs.

`metrics_frozen.json` was fixed before common-route predictions. Its primary
technical event-detection denominator includes all 13 case IDs. Every
prediction must have `status=ok` and a finite `event_score` in [0,1] for
TP/FN/FP/TN, sensitivity, FPR, and PPV to be reported. Any abstention or
missing result keeps the denominator and blocks those aggregate metrics with
`null`; it is never imputed as 0. Route scores of 1/0 are uncalibrated decisions,
not probabilities. Secondary fields are likewise blocked when unavailable.
No AUPRC, confidence interval, original-read comparison, or biological
collapse accuracy is reported.

Reproduce using the `tandemx-dev` environment:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.build_b1_structure_cases \
  --protocol benchmarks/controlled_collapse/v2/protocol.json \
  --repo-root . --outdir /tmp/b1v2_replay
conda run -n tandemx-dev python -m benchmarks.scripts.score_b1_structure_cases \
  --bundle benchmarks/controlled_collapse/v2/development_bundle \
  --metrics benchmarks/controlled_collapse/v2/metrics_frozen.json \
  --predictions benchmarks/m2_routes/alignment/evidence/development_b1v2_predictions.jsonl \
  --outdir /tmp/b1v2_alignment_replay
conda run -n tandemx-dev pytest -q tests/unit/test_b1_structure_cases.py
```

`alignment_score/` and `graph_score/` are independent scoring receipts for
the frozen input-only route predictions. Both have 12 `ok` and one `abstain`;
both formal all-case primary aggregates are blocked. The TSV
`route_comparison_per_case.tsv` is a per-case audit, not a route ranking.
The perfect 11 TP/2 TN fixture in unit tests is constructed directly from the
ledger solely to test scorer arithmetic; it is a tautological smoke test and
is never model performance.

Source eligibility: no original reads were matched to this archived YSD56
candidate. Even a source-qualified assembly sequence cannot by itself prove
physical missing copies. The injected edit is exact technical truth only.
