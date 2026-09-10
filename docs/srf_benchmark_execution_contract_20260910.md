# SRF benchmark execution contract — 2026-09-10

## Status and purpose

This is an implementation contract, not a benchmark request. It was checked
against `benchmarks/scripts/run_srf_pilot.py`,
`benchmarks/challenge/simulate.py`, existing SRF receipts, and the current
simulator outputs. It does not generate seeds, write datasets, invoke SRF, or
provide a command line to run.

The contract separates a finite **development parameter-selection set** from a
single later **hold-out set**. No SRF k/count choice, HOR equivalence rule, or
success margin is frozen for hold-out until the development outputs and their
mechanical scoring limitations have been reviewed.

## Current executable interface

`run_srf_pilot.py` accepts only four public inputs: a list of dataset
directories, a tools root, a new output directory, and one or more positive
`--minimum-counts`. Every dataset must provide `reads.fa`, `truth_arrays.tsv`,
and `truth_reads.tsv`. The runner itself fixes the native workflow to:

| Item | Actual current value | Consequence |
| --- | --- | --- |
| KMC input mode | `-fm` FASTA | It does not accept FASTQ through this runner. |
| k-mer size | `-k151`, hard-coded | `k=101` cannot be claimed as an executable comparison condition without a reviewed runner change. |
| threads | `-t1`, hard-coded | Resources are single-thread direct-child measurements. |
| KMC memory option | `-m2`, hard-coded | This is KMC's observed accepted minimum on this host, not measured peak RSS. |
| count cutoff | `-ci<minimum_count>` | The sole present parameter-selection interface. Positive integers only. |
| counter ceiling | `-cs100000` | Saturation must be checked from native outputs/receipts if a future input can approach it. |
| timeout | 120 seconds per native command | A timeout is a failure, never a zero-prediction result. |
| downstream stages | SRF, `enlong`, minimap2, `paf2bed`, `bed2abun` | End-to-end timing includes all seven native stages. |

The runner writes commands, stderr, native files, `workflow_receipt.json`,
normalized `predictions.tsv`, `metrics.json`, `summary.tsv`, source snapshots,
and `validation.json`. It refuses existing output directories.

## What 500,000 bases means here

The default challenge `Scenario` has 100 reads of 5,000 bp, hence 500,000
input bases. With defaults it places arrays in 70% of reads, one array per
positive read, 12 copies of a 171-bp period, selected among three families.
This is read-level array occupancy, **not a whole-genome coverage estimate**.
KMC `-ci` is an observed k-mer occurrence threshold, so it cannot be set from
“500 kb” alone. It depends on array occupancy, number of copies, read length,
error/divergence, family sharing, k, and the simulated sampling process.

The existing successful pilot demonstrates that `k151` with `ci20` and `ci100`
can produce catalogues for two particular 500-kb datasets. It does not establish
that `ci20` is a generally detectable lower bound or a valid future hold-out
setting. The SRF paper and official README also state that k and cutoff depend
on input/read coverage; the paper reports a maize CentC rescue at k=101 after
k=151 did not recover it. These facts motivate development calibration, not
post hoc choice of a winning setting.

## Existing development evidence and failure handling

The development receipt at
`/Volumes/T7/Codex/TandemX/results/srf_workflow_development_v1_20260906/`
contains 32 dataset/cutoff rows: 17 `ok`, 4 `no_catalogue`, and 11 `failed`.
The failed rows include both `ci20` and `ci100` under divergent, high-error,
short-array, AT-rich, and some substitution/indel conditions. The separate
guarded workflow records empty eligible-k-mer dumps as
`no_eligible_kmers` and skips native SRF graph assembly because upstream SRF
asserts on an empty dump. These are three distinct outcomes:

| State | Required treatment |
| --- | --- |
| `ok` | Score raw and normalized output; retain all native products. |
| `no_catalogue` | Native SRF completed but emitted no motif; score an empty normalized prediction set only with that status retained. |
| `no_eligible_kmers` | Guarded precondition: no eligible KMC k-mers; do not invoke native SRF and do not call it an SRF zero result. |
| `failed` or timed out | Preserve command/stderr/receipt; no accuracy number is imputed. |

Before any hold-out, the development selection report must list every tested
k/count configuration, native status, counter-saturation check, and resource
profile. It must fix one condition (or a predeclared reporting set, without
choosing the best result) and save its selection receipt. A k=101 experiment
requires a reviewed extension of the runner and tests; it cannot be silently
substituted for the hard-coded k=151 workflow.

## Available scoring today

For `ok`, `no_catalogue`, and guarded `no_eligible_kmers` states, the runner
uses `benchmarks.challenge.evaluate.score_arrays` and
`benchmarks.challenge.sequence_metrics.score_cyclic_recovery` after execution.
Existing fields include:

* array recall/precision/F1 at the scorer's current interval matching rule;
  read detection recall/precision and negative-read call rate;
* base-union recall/precision/F1, predicted/truth union bp, duplicate bp;
* matched period and boundary MAE;
* cyclic monomer recovery, distinct consensus count, homologous-consensus
  fraction, unmatched consensus count, and best cyclic edit similarity;
* native catalogue size, motif lengths, out-of-scope period/span counts, input
  and truth hashes; and all workflow resource/status fields.

These endpoints are suitable for raw satellite motif/interval behavior under
the current challenge's declared 30--1000 bp period and >=100-bp span scope.
They are not a complete SRF abundance comparison. `bed2abun` output is saved,
but no code currently matches each native motif/HOR to a truth family and then
computes per-correspondence mapped-bp abundance error. There is also no paired
TandemX-versus-SRF comparison table, CI procedure, or predeclared superiority/
noninferiority decision code for this endpoint.

## HOR/family equivalence: required review interface

SRF may validly emit a higher-order unit while TandemX may emit a constituent
monomer. The current `score_cyclic_recovery` compares each SRF consensus to a
single truth-family sequence and therefore cannot alone adjudicate that case.
Before development scoring is finalized, the main task must approve a
versioned correspondence interface with at least:

| Required field | Reason |
| --- | --- |
| `srf_motif_id`, sequence hash, length, native abundance row | Binds the adjudication to immutable native output. |
| one or more truth family IDs and their sequence hashes | Makes a many-to-one HOR relationship explicit. |
| relationship class (`monomer_match`, `hor_contains_monomer`, `related`, `unmatched`) | Prevents a valid HOR from becoming an automatic false negative. |
| cyclic/rotation and reverse-complement matching procedure and thresholds | Avoids hand-waved motif identity. |
| interval overlap/base-support evidence | Keeps sequence equivalence separate from read-array recovery. |
| adjudicator, rationale, and `predeclared` flag | Distinguishes a frozen rule from post-result manual rescue. |

Until this interface and scoring implementation are approved, report current
cyclic metrics as descriptive and preserve SRF HORs separately; do not convert
them to TandemX monomers or claim family-level abundance parity.

## One-time hold-out gate (future only)

The hold-out may begin only after all of the following are archived: selected
SRF parameter condition(s), runner revision/hash, scenario definitions,
HOR correspondence/scoring revision, resource/failure policy, and comparison
metrics. The later hold-out must use unobserved seed/family material generated
only after that archive, run both tools on identical observable `reads.fa`, and
keep truth inaccessible to external commands. A failure or resource limit ends
the attempt with a retained receipt; it does not trigger parameter expansion.

No hold-out scale, number of datasets, seed, k, cutoff, resource cap, or
command is specified here because each would currently be unvalidated or not
supported by the present runner. Those values belong in the reviewed development
selection receipt, not in a document that could accidentally become a post hoc
benchmark recipe.
