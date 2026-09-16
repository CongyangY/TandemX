# B1 toy software verification receipt

This is a synthetic two-locus software fixture. It has no real reads or
biological truth. The first eight-copy array is pure monomer A; the second
eight-copy array mixes A and one-base-different B units. The 24-bp nonrepeat
deletion equals their combined 25% contraction (12 + 12 bp). The donor-swap
control leaves the FASTA unchanged and marks pairing invalid. Generated
FASTAs, exact family ledgers, event ledgers, chains, and SHA-256 receipts are in
`evidence_v1/`; `score_smoke_v1/` uses a header-only prediction TSV to check
that no call is represented as `not_reported`, not as zero.

Commands run from repository root, inside the dedicated environment:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.controlled_collapse generate \
  --source benchmarks/controlled_collapse/toy_assembly.fa \
  --plan benchmarks/controlled_collapse/toy_plan.json \
  --outdir benchmarks/controlled_collapse/evidence_v1
conda run -n tandemx-dev python -m benchmarks.scripts.controlled_collapse score \
  --generated benchmarks/controlled_collapse/evidence_v1 \
  --predictions benchmarks/controlled_collapse/empty_predictions.tsv \
  --decision-threshold 0.5 \
  --outdir benchmarks/controlled_collapse/score_smoke_v1
conda run -n tandemx-dev pytest -q tests/unit/test_controlled_collapse.py
conda run -n tandemx-dev pytest -q
```

Results on 2026-09-16: focused tests 13 passed after the last scoring change;
the previous full repository regression passed **914 tests** before this
bounded scoring change. The
receipt has 10 cases × 2 loci = 20 fixed-denominator rows. Twelve rows carry
injected positive deletion, including two completely absent arrays; all 20
smoke predictions are `not_reported`, and classification/MAE and primary AUPRC
are null. Total ledger denominator is **20**; the valid-pair edit-metric
denominator is **18**, excluding the two invalid donor-swap rows while keeping
them in the total and status counts. A separate full-`ok` test copies the known
ledger truth into predictions and produces TP=12, FN=0, FP=0, TN=6,
sensitivity=1, FPR=0, and injected edit MAE=0 at the preregistered toy
threshold 0.5. This is a **tautological scorer smoke test**, not a model
performance result. The generated
fixture is not evidence of TandemX detection accuracy.
An additional independent test plan adds an eight-copy 8-bp locus to exercise
different unit lengths; it is not part of the committed 10-case toy receipt.

Input SHA-256 values:

- `toy_assembly.fa`: `20ba97916ccb680059a2c653a8d76b3250514c765c88649cf32b0e34cd2078eb`
- `toy_plan.json`: `d852f78a1397273f1cb1f96c909d5c064bce0f6b0aa8ed9749ba9acc1d0ea744`
- `controlled_collapse.py` generator: `f5679bc956f439e0ea5e72d64ca790eb817c871b2d035099e178037008018e93`

The full output SHA-256 values are stored in `evidence_v1/receipt.json`.
Scoring checks them before joining predictions.
