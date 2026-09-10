# SRF development evidence inventory

## Scope

This is a read-only inventory of all located SRF development outputs.  It does
not rerun SRF, create seeds or datasets, or convert a development observation
into a formal comparison.  The authoritative run directories are under
`/Volumes/T7/Codex/TandemX/results`; small files under `paper/evidence/srf_*`
are repository copies or selected artifacts, not additional independent runs.

All three runs use the documented seven-stage native workflow in
`benchmarks/scripts/run_srf_pilot.py`: KMC k=151 counting, dump, SRF assembly,
`enlong`, minimap2 mapping, `paf2bed`, and `bed2abun`.  It uses one thread,
KMC `-m2`, count cutoffs ci100 or ci20, 120-s stage timeout, and observed FASTA
bases as the abundance denominator.  The runner retains native commands,
stderr, FASTA/BED/abundance products when produced, per-stage timing/RSS and
input/executable hashes in receipts.

The suite is development-only: all datasets use `s1101`, are approximately
500,000 input bases, and have no held-out split.  It is a sequence/array and
read-detection exercise.  Native mapped-bp abundance is retained, but there is
no paired TandemX-versus-SRF abundance error, calibrated copy-number endpoint,
confidence interval procedure, or predeclared superiority decision.

## Located executions

| Result directory | Input/split | Statuses | Useful observed endpoints | Resource receipt |
| --- | --- | --- | --- | --- |
| `results/srf_workflow_pilot_v1_20260906` | Two selected development scenarios, `clean_171_s1101` and `indel_01pct_s1101`, each ci100/ci20; 4 runs | 4 `ok` | Native catalogue/periods, cyclic recovery, array and read-detection metrics, negative-read calls, union bases, native abundance | `workflow_receipt.json` per run; wall 1.923--1.961 s and max sequential native-stage RSS 48.1--51.4 MiB |
| `results/srf_workflow_development_v1_20260906` | 16 scenario names × ci100/ci20; 32 runs | 17 `ok`, 4 `no_catalogue`, 11 `failed` | Same metrics for `ok`; `no_catalogue` is a completed native discovery with empty motif output; failed attempts retain command/stderr | Per-run receipts and 164 saved stage stderr logs; successful wall 1.863--2.026 s, RSS 45.2--51.5 MiB |
| `results/srf_workflow_empty_guard_v1_20260906` | Same 16 scenarios × ci100/ci20; 32 runs; rerun with the disclosed empty-dump guard | 17 `ok`, 4 `no_catalogue`, 11 `no_eligible_kmers`, 0 process failures | Same metrics for `ok`; empty-dump outcomes are explicit guarded skips, not SRF zero-output calls | Per-run receipts, input/executable/environment hashes and `validation.json`; successful wall 1.951--2.086 s, RSS 45.6--53.3 MiB |

The 16 named scenarios are `at_rich_control`, clean periods 60/171/421/729,
`dispersed_control`, `divergent_units`, indel 0.01/1/4%,
`low_complexity_control`, `related_families`, `short_arrays`, substitution
1/5%, and `two_arrays`.  The full names and every status are in each
`summary.tsv`.

## Status interpretation and retained evidence

The original 32-run development suite exposed 11 native SRF failures after an
empty KMC dump reached SRF's `ca_kmer_read` assertion.  The later guarded suite
does not relabel those as native SRF successes or catalogue negatives.  It
records `no_eligible_kmers`, skips SRF and downstream stages, and creates no
fake `srf.fa`.  The four `no_catalogue` rows are distinct: native SRF was
invoked successfully but emitted no motif, so their normalized prediction set
is empty with that status retained.

The guarded suite's `validation.json` reports `complete: true`, 32 runs and
zero process failures.  That records workflow completion after the guard; it
does not erase the prior 11 native assertion failures.  For example, ci20
recovers arrays in many conditions, while divergent units and severe indel or
substitution scenarios visibly retain no-eligible-kmer outcomes.  All negative
control rows in the guarded suite have zero negative-read calls under the
scoring scope, but this is a small constructed development result, not an SRF
specificity estimate.

For every `ok` row, the common scorer limits cross-tool array scoring to native
positive keep flags, motif periods 30--1000 bp and spans at least 100 bp.
Out-of-scope motifs, including possible HORs, remain in native output and are
counted rather than decomposed or scored as false calls.  This preserves SRF's
native result but leaves a family/HOR correspondence question unresolved.

## What is available now

The inventory supports these limited factual statements:

* SRF and its required KMC/minimap2/k8 stages were built and executed on the
  recorded macOS development inputs; executable hashes and source snapshots are
  saved in each `environment.json`.
* The workflow has observed success, completed-empty-catalogue and
  no-eligible-kmer statuses, plus the retained original native assertion
  failures.  It does not fabricate a catalogue for empty count input.
* Sequence recovery, array/read detection, union-base and negative-read
  endpoints are available for in-scope native outputs, and raw `bed2abun`
  output supplies a native mapped-bp abundance field.
* One-run resource observations exist.  They are sequential direct-child RSS
  maxima and wall times, not replicated performance estimates or peak scratch
  disk measurements.

## Formal gap before a fair comparator result

The existing evidence is not a final TandemX-versus-SRF comparison because it
lacks all of the following:

1. A fresh, prespecified held-out split distinct from all `s1101` scenarios.
2. A finite development calibration record for SRF k/count settings and the
   rule used to freeze any later condition; ci100/ci20 are documented probes,
   not an established optimum.
3. A locked family/HOR-aware correspondence table.  A valid SRF higher-order
   motif cannot be silently converted into a TandemX monomer miss or match.
4. A common, explicitly normalized abundance endpoint.  SRF mapped-bp fraction
   is not a calibrated physical copy-number estimate, and no paired error table
   exists.
5. Matched, repeated resource measurements and a formal failure-accounting rule
   for both methods on the same inputs.

Until those elements exist, cite the directory-specific observed outcomes and
their statuses rather than a comparator win/loss or an external-tool accuracy
claim.
