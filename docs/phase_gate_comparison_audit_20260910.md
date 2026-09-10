# Phase-gate comparison audit (development only)

## Scope and conclusion

This is a static audit of `benchmarks/challenge/local_phase.py` and
`benchmarks/scripts/run_phase_gate_development.py` at the working tree state
on 2026-09-10.  It did not execute the development runner or create or access a
holdout.  With the runner's explicit `haploid_depth=1.0`, A0 is in the same
sampled-read scale as A1--A3 before all are divided by observed coverage for
the synthetic source-bp classifier.  This is an intentional fixture convention,
not the genome-normalized public interpretation of a normal quantification
run.  The metric units are therefore comparable within this development
fixture, with the provenance and endpoint boundaries below.

The A1 implementation is a multi-target, per-base ambiguous mapping baseline.
It intentionally retains qualifying competing secondary PAF rows and excludes
only spans on which family assignments overlap.  Its assignment scope differs
from the repository's older conservative ONT occupancy audit, which excludes an
entire multi-family read.  Neither policy should be silently substituted for
the other.

No phase-gate value may be described as better or worse than the historical
0.363222 MARE.  That number is from the frozen three-genome, 1,485-family-row
quantify validation, whereas this fixture has two synthetic families, ten
development conditions, two development seeds, and sampled-read interval
endpoints.

## What each path actually measures

| Path | Observable input and operation | Value supplied to scoring | Unit at the point of scoring | Audit status |
| --- | --- | --- | --- | --- |
| A0 | Calls the public `quantify_toy_copy_number` code with the fixture reads and catalogue (`k=21`, Python backend, explicit haploid depth 1). | Reconstructed `estimated_copy_number * monomer_length` from `A0/copy_number.tsv`. | Sampled-read-scale diagnostic-k-mer depth times monomer length: an explicit depth of 1 prevents the normal depth normalization. | Comparable within this fixture after the runner's common observed-coverage conversion. |
| A1 | `minimap2` maps reads to each representative repeated 20 times; qualifying query intervals are grouped by target family. | Query/read intervals from PAF. | Sampled-read bp after the metric's family-wise interval union and cross-family overlap exclusion. | Same domain as A2/A3 interval metrics, with the assignment scope below. |
| A2 | Exact circular seed/phase-diagonal chains from `PeriodicContextPrototype`. | Read-coordinate intervals. | Sampled-read bp. | Comparable to A1/A3 for the per-base mask metrics. |
| A3 | `LocalPhasePrototype` performs a bounded predecessor DP using circular template advance, maximum gap/drift and phase coverage. | Read-coordinate seed intervals; seed-free spans are deliberately not filled. | Sampled-read bp. | Comparable to A1/A2 only as a development interval-occupancy rule. |

`score_predictions` appropriately computes precision, recall and relative error
for its intended input: predicted *sampled-read* unique bp versus truth
*sampled-read* bp.  It merges intervals within a family before sweeping events,
so overlapping mappings of the **same** family are not double counted.

## A0 unit check and output naming limitation

`quantify_toy_copy_number` computes `estimated_copy_number = median_depth /
haploid_depth`, then `estimated_bp = estimated_copy_number * monomer_length`.
Here `haploid_depth` is explicitly 1.0, so `a0_bp` is `median_depth *
monomer_length`, rather than a normal genome-scale copy-number estimate.  Under
the fixture's uniform sampled-read design this is a sampled-read-scale
quantity, like A1--A3's eligible read bp.  Dividing all of them by
`read_bases / genome_size` converts each to the fixture source-bp scale once;
there is no double normalization in this runner.

The runner nevertheless reconstructs the field from
`estimated_copy_number * len(catalogue[family])` instead of directly reading
the native `estimated_bp` column.  It should not imply that this fixture-only
depth-one quantity has the usual public `estimated_copy_number` biological
meaning.  Any reported result must say that A0 is a development endpoint with
explicit depth one, and must retain the source-bp conversion rule.  This audit
does not select an endpoint or alter parameters.

## A1 assignment scope

The A1 command maps to both 20x repeated targets with `-N 50`.  The parser
accepts every PAF row with block length at least 100 bp and identity at least
0.90.  It intentionally does not inspect `tp:A:P`, mapping quality, or a
best-alignment score: competing secondary mappings participate so ambiguity is
represented rather than forced to the highest-scoring family.

The subsequent per-base scorer does handle overlap correctly at a narrow level:

* intervals within one target family are merged, so repeated/supplementary
  same-family PAF intervals do not create double bp;
* segments concurrently covered by two or more families are recorded as
  ambiguous and excluded from all family unique/TP/FP totals; and
* total read bases are conserved across unique, ambiguous and unassigned bins.

Different-family mappings on non-overlapping portions of one read are retained
as two family contributions.  This is the stated per-base ambiguity policy, not
a double-counting defect; each accepted segment is attributed once unless it
overlaps a competing family.  It differs in scope from the repository's
existing ONT occupancy summarizer
(`evaluate_macadamia_ont_occupancy.py`), which accepts `tp:A:P` rows only,
unions query intervals within a family, and excludes a whole read when its
accepted mappings name more than one family.  A1 should be described as a
multi-target PAF-overlap, per-base ambiguity baseline; its results cannot be
equated with the whole-read-exclusion ONT endpoint.

## Baseline and provenance boundary

The receipt records `baseline: "81827c3"`, and A0 imports the public quantifier
function.  The runner does not execute a pinned checkout, however, and
`source_sha256` contains the runner, local phase, fixture and metric files only.
It omits `tandemx/quantify/mvp.py`, `tandemx/utils/kmers.py`, and the imported
`context_prototype.py`.  Commit `81827c3` itself changed `mvp.py` and
`kmers.py`.  The current static audit found no diff for `mvp.py` against that
commit, but the receipt does not itself establish the executed baseline bytes,
and it does not cover the other imported modules.  This is a reproducibility
limitation, not a claim that the current A0 path was a different implementation.

## Novelty boundary of A3

A3 is more specific than ordinary exact occurrence counting: it uses the
catalogue's circular monomer phase, requires forward template-coordinate
recurrence, applies bounded local drift, and gates terminal chains by phase
coverage and minimum units before counting their seed spans.  Those mechanics
can be falsified by comparing accepted/rejected occupancy against the fixed
decoys and by testing whether benefits persist after the baseline is genuinely
competitive.

They do not by themselves establish a new chaining algorithm.  The DP is still
a conventional bounded predecessor chain with a gap/drift penalty, and the
"reset" is a per-edge penalty rather than an explicit fitted collection of
local phase-offset segments or a change-point model.  Circular/periodic phase
constraints and seed chaining are established mechanisms.  Any future claim
must therefore rest on a predeclared, family-conditioned abundance eligibility
endpoint that outperforms a correctly implemented ordinary baseline on the
same endpoint, with adverse shared-fragment and indel outcomes retained.

## Reporting boundary

The only defensible present description is: a development-only mechanism
experiment with a common sampled-read/source-bp fixture endpoint, and a
deliberately per-base-ambiguous A1 assignment rule.  Do not compare any
phase-gate metric to the historical 0.363222 MARE, and do not use this audit to
change the frozen method, thresholds, parameters, or holdout design.
