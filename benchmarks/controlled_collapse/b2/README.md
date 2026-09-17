# B2: synthetic structural hold-out

B2 is a separately seeded **synthetic structural hold-out**, not a biological
validation set or an original-read/assembly comparison. Its design borrows the
100/200/300/400 bp monomer length range and nested HOR motivation from
[HiCAT (Gao et al., Genome Biology 2023)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02900-5).
It is an independent generator and parameter set, not a reproduction of the
published HiCAT simulation. B2 was defined after the Col-CEN development
results but before any B2 route prediction; it cannot support a claim of
independent biological replication.

`protocol.json` fixes four distinct synthetic founders, 5 monomers per HOR,
four canonical HOR copies, an inserted `BCBC` nested segment in the source,
0.5% array-only read substitutions, 0.1% single-base array-only indels,
protected engineered flanks, and 13 cases per founder. A separate source
generator writes `source_manifest.json` with all source copy sequences and
their hashes. The case generator then writes 52 input and causal truth rows:
44 edited positives and eight intact controls. The 52 causal rows contain 48
distinct assembly sequences because each founder has two identical intact
assemblies under different read-support/error settings. They are not 52
biological replicates. Source, case and read-error seeds are distinct from
Col-CEN v3. No real genome or raw-read sequence enters B2.

The input-only `held_out_bundle/inputs.jsonl` uses schema version 2 and
`split=synthetic_held_out`. It supplies assembly sequence, simulated reads,
monomer templates, unique engineered flanks, error profile and read support.
Twenty-four cases have only one simulated read; 28 have three. Error events
are independently generated per read in noisy cases, and exact controls have
identical copies. The separate `truth.jsonl` contains the exact source and
edited label paths, copy mapping, orientation, bp delta, coordinate chain and
per-read error event ledger. All coordinates are 0-based half-open bp in the
full engineered molecule. `receipt.json` hashes the source manifest,
generator, inputs, truth and each assembly FASTA. The truth file is separated
from route inputs for masked execution; committing it is auditability rather
than secure blinding.

`metrics_frozen.json` fixes a 0.5 score threshold before prediction. In the
all-case intent-to-diagnose denominator, a positive refusal or missing row is
FN; an intact refusal is `unresolved_negative`. FPR counts explicit false
positive calls over all eight intact controls, while the negative failure rate
also counts unresolved controls. Eligible-only confusion and coverage are
reported separately. The scorer preserves N/A for missing subtype, signed-bp,
breakpoint, HOR-period or order predictions. For B2 only, the source canonical
HOR period of five and synthetic monomer label path are known exactly, so
those secondary metrics can be computed if a route reports the required
fields. Runtime and peak RSS require measured route receipts. Scores of 1/0
are uncalibrated decisions, not probabilities; no AUPRC or biological copy
accuracy is inferred.

The current graph route adapter explicitly requires `split=development`, so
it may refuse B2 at the input gate. B2 must keep its truthful held-out split;
an interface refusal is reported as such, without modifying the frozen route.
Some arrays exceed the alignment route's declared 4,096 bp ceiling. Those
cases remain in the denominator as prespecified operating-limit probes.

Reproduce in `tandemx-dev`:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.generate_b2_source \
  --protocol benchmarks/controlled_collapse/b2/protocol.json \
  --out-manifest /tmp/b2_source_replay.json
conda run -n tandemx-dev python -m benchmarks.scripts.build_b2_cases \
  --protocol benchmarks/controlled_collapse/b2/protocol.json \
  --source-manifest benchmarks/controlled_collapse/b2/source_manifest.json \
  --outdir /tmp/b2_bundle_replay
conda run -n tandemx-dev pytest -q tests/unit/test_b2_cases.py
conda run -n tandemx-dev python -m benchmarks.scripts.score_b2_cases \
  --bundle benchmarks/controlled_collapse/b2/held_out_bundle \
  --metrics benchmarks/controlled_collapse/b2/metrics_frozen.json \
  --predictions /path/to/frozen_route_predictions.jsonl \
  --outdir /tmp/b2_score_replay
```

The test fixture that predicts 44 TP and eight TN reads `truth.jsonl`
directly. It is a tautological scorer arithmetic check, never route
performance. Route predictions must be generated only after this protocol,
bundle, metrics, scorer and tests are committed.
