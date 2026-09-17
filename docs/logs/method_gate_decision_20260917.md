# TandemX method gate, 2026-09-17

## Decision scope and provenance

This decision uses frozen M1 diagnostics (`993d441`), frozen B1 v2/v3
protocols, two independent M2 research prototypes, a simple paired-flank
length baseline, and B2 source/scorer commit `6c4df79`. B2 predictions were
made after that commit without changing either route. B2 is a **synthetic
structural hold-out**, separate in source/seed/monomer length from Col-CEN v3;
it is not a held-out biological donor or raw-read validation.

The full B2 denominator is 52 cases: 44 injected edits and 8 intact controls,
across four independently generated 100/200/300/400 bp founder templates.
Each source includes a nested BCBC block. Twenty-eight cases have three
synthetic read IDs; errorful cases perturb these independently, while exact
controls may repeat a sequence. Twenty-four cases have one read. Only array
interiors have simulated substitutions/indels; artificial flanks are protected.

The baseline initially emits `status=insufficient_read_support` for its 24
one-read cases. The precommitted B2 scorer accepts schema-v2 `abstain`, so a
truth-blind presentation adapter changed that **status string only** to
`abstain` in a separate archived file, retaining the raw predictions and
null scores. No event decision or threshold changed. The scorer's first
attempt rejected the raw string and wrote no result; the interface mismatch
is preserved in this log.

## Gate A: absolute family abundance

**NO-GO as a current Methods novelty claim.** M1 candidate estimators did
not stably beat the strong mapping baseline in the frozen tournament.
Signature rank analysis found exact non-identifiability for identical
families, weak/no family-exclusive markers in partially identifiable cases,
and family-specific rejection under representative divergence. The current
estimator remains a family-level evidence / assembly-audit score; it is not a
calibrated physical copy-number estimator. No additional abundance estimator
was opened.

## M2 route tournament and selection

The local edlib monomer-path route is the **selected research candidate among
the two M2 routes** because it issued decisions on native Col-CEN development
arrays where the exact-transition graph route rejected every case. That
selection is limited to research follow-up; it is not promotion to a public
method or a Methods-paper core contribution.

On Col-CEN B1 v3 development (39 edits/controls from three intervals of one
assembly lineage, exact synthetic full-span reads), alignment detected 24/33
events and missed 9; the paired-flank length baseline detected 18/33 and
missed 15; neither reported a false event among six intact controls. Alignment
rejected three cases and missed all three local rearrangements and all three
within-tile compression events. The transition graph rejected all 39 cases.
These repeated edits are not independent biological replicates.

On the prespecified B2 synthetic hold-out, intent-to-diagnose results are:

| Route | TP/44 | FN/44 | FP/8 | TN/8 | Unresolved intact | Decision coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Alignment M2 | 6 | 38 | 0 | 1 | 7 | 7/52 |
| Transition graph M2 | 0 | 44 | 0 | 0 | 8 | 0/52 |
| Paired-flank length baseline | 20 | 24 | 0 | 4 | 4 | 28/52 |

The reported FPR is 0/8 for all routes, but this must be read alongside
negative failure rates: alignment 7/8, graph 8/8, baseline 4/8. Eligible-only
alignment sensitivity of 6/6 is a selection-conditioned result from just 7
decisions and must not replace all-case 6/44. Alignment's six TPs all come
from the 100 bp founder; it returns no decisions for 200, 300, or 400 bp
founders. Its 45 abstentions comprise 38 cases over the 4,096 bp array-input
limit, six with insufficient resolved reads, and one 400 bp monomer beyond the
declared 300 bp monomer limit. In the 100 bp founder it detected one inversion
missed by the length baseline, but both missed low-support cases. The graph's
52 abstentions come from its frozen `split=development` input contract; the
held-out split was not relabeled to make it run.

Core per-case elapsed sums were approximately 12.22 s for alignment and
0.0123 s for the baseline; process peak RSS was 22.7 MB and 18.8 MB,
respectively. These are tiny synthetic workloads, not large-genome throughput
measurements. The graph's near-zero rejection time is not inference speed.

**Gate B: NO-GO for the frozen M2 architecture-audit route under the declared
100–400 bp synthetic structural operating range.** The selected alignment
candidate detects fewer edits and covers fewer cases than the simple length
baseline. The exact graph route is inapplicable. Real-read, donor-held-out
performance remains unmeasured, so this is a negative method-development
decision rather than a biological accuracy estimate. Do not use v3's more
favorable result to override the B2 stop-loss.

## Assembly-only competitor endpoint

One engineered 171 bp monomer control supports a TideCluster array-localization
score (base recall 1.0, precision approximately 0.99994) and native 855 bp
period output. CENdetectHOR ended at recorded upstream technical failures.
TandemX's supplementary `locate` run received truth-derived candidate
monomers and covered 98.73% of truth array bases; this different prior prevents
a fair cross-tool accuracy rank. A common no-prior monomer/HOR-order endpoint
is not completed. Assembly-only tools must not be marked false negatives on
the separate assembly-versus-read discordance task.

## Gate C and manuscript consequence

No real-material TandemX prediction here has independent raw-read/orthogonal
structural confirmation. **Gate C remains untested.** Biological utility and
superiority claims are not supported. With Gates A and B failing, there is no
current evidence basis for a Genome Research Methods core-method claim. The
evidence may support a bounded assembly QC/resource application after a
separate, appropriately validated scope decision. Formal release, Bioconda,
Zenodo, and manuscript promotion remain paused under project instructions.

## Reproduction and audit

See `benchmarks/controlled_collapse/b2/{protocol.json,metrics_frozen.json,
held_out_bundle/receipt.json,alignment_score/summary.json,graph_score/
summary.json,baseline_score/summary.json}`. Raw and normalized baseline
predictions are both retained. The independent hostile review is in
`docs/logs/hostile_method_gate_20260917.md`. The full local `tandemx-dev`
pytest suite passed 1,028 tests in 118.68 s on the final B2 code and scoring
state.
