# Controlled array contraction benchmark B1 (development prototype)

This independent benchmark edits a copy of a supplied assembly. It does not alter
TandemX abundance inference, claim physical copy-number truth, or create a public
`tandemx` command. The first runnable fixture is synthetic and is a software
integrity check, **not** a real-data validation result. Large T7 inputs remain
paused after the documented I/O failure.

## Inputs and provenance gate

Run from the repository root in `tandemx-dev`:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.controlled_collapse generate \
  --source benchmarks/controlled_collapse/toy_assembly.fa \
  --plan benchmarks/controlled_collapse/toy_plan.json \
  --outdir /tmp/tandemx-controlled-collapse-toy
```

The source is a FASTA with unique record names. The JSON plan names the
assembly/read donor and haplotype IDs and states whether independent baseline
read consistency has been checked. It nominates nonoverlapping arrays as
`family_id`, `locus_id`, `contig`, `start`, `end`, `unit_bp`; coordinates are
zero-based, half-open and lengths must be multiples of four full units. The
`nonrepeat_interval` is a preselected region outside nominated arrays. Optional
`unit_labels` and `unit_monomers` specify per-copy truth and are checked
against source sequence before editing; the included fixture has a pure A
array and a second A/B mixed array sharing the A monomer. This exposes a
family-mixture/shared-sequence test condition without claiming that a tool
can resolve it. This
small prototype caps memory by buffering one nominated interval, but it does
not yet verify tandem identity, annotation quality, source assembly
completeness, or donor metadata against external records. Input hashes are
recorded; the user must perform source/provenance QC before using real data.

The editor produces ten cases from the same source: 0%, 25%, 50%, 75%, and
100% right-boundary contraction of **each** nominated locus; 25% internal
deletion; 25% left-boundary truncation; an untouched-array control with a
nonrepeat deletion equal to the combined 25% contraction size in the toy
fixture; 25% expansion; and a donor-swap pairing control. The last case leaves
the assembly untouched and flags the read/assembly pair as invalid. Multiple
homologous loci may share a family ID but must have distinct locus IDs; all
loci are kept in the ledger, including completely absent arrays. A control is
useful only when its interval is independently confirmed nonrepeat.

The toy fixture has two eight-copy loci and a separate control contig. It has
no simulated or measured reads, so its pairing state is `unverified_pair` and
its outputs cannot assess read–assembly abundance deficit.
The focused tests also run an independent three-locus plan with an 8-bp unit
beside the two 6-bp loci; that test is separate from the versioned toy receipt.

## Output schema

Each case writes `<case>.fa`, `<case>.ledger.tsv`, `<case>.chain.tsv`, and
`<case>.events.tsv`.
`receipt.json` records input paths and SHA-256 hashes, generator source hash,
case IDs, levels, and output hashes. FASTA headers retain source contig IDs.

`ledger.tsv` has one row per nominated family/locus, even at 100% deletion or
in the negative controls. `source_start/end` and `edited_start/end` are array
intervals; `unit_bp`, `source_full_units`, `retained_full_units`, and
`removed_full_units` describe complete annotated units. Residual columns are
currently zero because this prototype accepts full-unit arrays only.
`removed_bp` and `inserted_bp` are exact edit amounts. The two array sequence
hashes are SHA-256 of uppercase ASCII sequence. Breakpoints are source left
and right join coordinates plus edited join coordinate. At zero deletion, the
left and right breakpoint coordinates coincide and represent no break.
`source_unit_labels` and `edited_unit_labels` are JSON lists; they are empty
when a plan does not provide per-unit truth.
`truth_scope` is always `injected_edit_delta_only`. `pairing_status` is
`unverified_pair`, `declared_same_donor_baseline_consistent`, or
`invalid_donor_pair`; even the declared state is metadata, not independent
proof of donor identity or physical copy number.

`chain.tsv` is a bidirectional **piecewise coordinate relation**, not a
single-valued map across deletions or insertions. `match` rows map retained
source `[source_start,source_end)` to edited `[edited_start,edited_end)`.
`delete` rows map a source interval to a zero-width edited junction; `insert`
rows map a zero-width source junction to added edited sequence. Coordinates
are zero-based and half-open. Each case's chain includes all FASTA contigs.
The event ledger records every nominated interval and the nonrepeat deletion
as separate operations with source/edited sequence hashes, exact breakpoint
coordinates, edit size, unit labels, and pairing state. A nonrepeat deletion
is an event, not a positive array-family truth row. Donor swap is an invalid
pairing event; it does not fabricate a different read set.

## Scoring contract

Supply a tab-separated prediction file with header
`case_id family_id locus_id status score predicted_missing_bp`. `status` is
`ok`, `abstain`, `unmapped`, `failed`, or `not_run`. Only `ok` may carry a score
in `[0,1]` and a nonnegative integer missing-bp estimate. Its score must be
finite; every non-`ok` state must leave score and bp blank. The scorer joins
the **full** case × family/locus ledger, writes `scored.tsv`, and assigns
missing predictions `not_reported`; it never replaces them with zero. It checks
all generated output hashes before scoring and rejects duplicate or unknown
prediction keys. Its `absolute_error_bp` is descriptive error against the
**injected edit amount only**, shown only for `ok` predictions with a valid
read/assembly pairing label. The summary
retains all statuses, positive and 100%-deletion row counts. Pass a
preselected inclusive `--decision-threshold` in `[0,1]`. If **all** rows with
a valid pairing label are `ok`, the summary reports TP/FN/FP/TN,
sensitivity, FPR, and injected missing-bp MAE on that entire fixed subset.
Otherwise every metric is null with `status=blocked`; abstentions, failures,
and missing rows are never treated as true negatives. The invalid donor-swap
rows remain in the overall denominator and status counts but are excluded
from paired edit-recovery metrics. The toy fixture's complete-prediction test
copies ledger truth into predictions and is a tautological scorer smoke test,
not a model accuracy result or independent statistical replication. The
summary states both the total family-row denominator and valid-pair metric
denominator, so the donor-swap rows cannot silently disappear. `primary_auprc`
remains null until a prospective evaluation freezes the eligible denominator,
the primary same-input mapping baseline, and statistical protocol. The
donor-swap case is invalid for a paired read–assembly inference comparison.

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.controlled_collapse score \
  --generated /tmp/tandemx-controlled-collapse-toy \
  --predictions path/to/predictions.tsv \
  --decision-threshold 0.5 \
  --outdir /tmp/tandemx-controlled-collapse-scored
conda run -n tandemx-dev pytest -q tests/unit/test_controlled_collapse.py
```

The edit amount is not the original biological missing-copy truth. Evaluating
the induced *change* in read–assembly deficit requires same-donor/haplotype
reads, independent baseline consistency checks, frozen localization and
matching rules, and a fixed denominator. Natural collapse requires separate
flank-anchored molecular truth. No biological accuracy or method advantage
follows from passing this toy software test.
