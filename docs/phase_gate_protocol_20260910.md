# Bounded phase-coherence gate

Baseline is `81827c3`. Accuracy lives on `codex/phase-coherence-gate-20260910`;
efficiency lives on a separate worktree/branch. Production code, formal
manuscript and release remain unchanged. No new detector, HMM, EM or assembly
repair is permitted. At most two accuracy revisions are allowed.

## Development design and units

Two new explicitly development source seeds, ten conditions, three source
coverages (2/5/10x), two supplied families per source yield 60 datasets and
120 family-condition observations. They are dependent conditions of two source
seeds, not 120 independent biological replicates. The generator rejects holdout
labels. Unit variation/indels are introduced into source repeats; this first
development test does not cover independent empirical HiFi read errors.

A0 is the frozen Python k=21 estimator with explicit haploid depth one. It
therefore measures sampled repeat bp in this experiment; dividing once by
observed source coverage yields the synthetic source-bp estimate used only for
binary scoring. A1 maps reads competitively to 20-copy supplied templates using
minimap2 map-hifi, -N50, -f1000, -r100,100, one thread. Alignments >=100 bp and
identity >=0.9 contribute; secondary mappings remain eligible. All methods
union within family and preserve cross-family overlap as ambiguous rather than
forcing the highest-scoring family. This is not the old ONT whole-read exclusion
policy. A2 is the strict prototype at baseline. A3 uses bounded local phase
changes, weighted seeds and template recurrence without full alignment.

Sequence-only, phase, phase-coverage, recurrence/order and re-phasing ablations
are retained, including a specificity-weight ablation. Template, k and mapper
settings differ between native methods and are reported. A0 has no interval
precision/recall; those fields are unavailable, never zero.

MARE here uses sampled read truth and is **not comparable** to the historical
0.363222 genome-abundance validation mean. Synthetic binary assembly amounts
are 0.3 of source truth for F1 and full source truth for F2. Sampling variance
and low source coverage can therefore cause false binary calls even when read
intervals are accurate. This is an estimator stress test, not an assembly
localizer benchmark or a real-genome accuracy claim.

The synthetic FPR is specifically a false call among the 60 F2-role rows in
this fixture, whereas the historical classifier FPR is FP/(FP+TN) over a
conditional assembly/read matrix. They are not commensurate metrics. See
`fpr_scope_reconciliation_20260910.md` for their exact numerator, denominator,
truth and decision-rule reconciliation.

## Round 1 and one permitted modification

Round 1 connects seeds with at most 20-bp local drift, 45-bp seed-start gap,
reset cost 1 plus 0.05 per drift base, and at most 64 predecessor candidates.
Acceptance requires at least 70% distinct template start phases and two units
of template progression. Counts include the union of exact seed intervals.
Native commands, all scores and protocol/source hashes are retained under T7
`results/phase_coherence_gate_development_r1_20260910`.

The first observed development mean is A0 0.2070, A1 0.0568, A2 0.2455,
A3 0.0822 MARE. A3 restores indel recall from A2 0.1079 to 0.8982, but does
not outperform A1. The exact A0 table values must be used in final scoring,
rather than reconstructing bp from its separately rounded copy-number field.

The only second-round change is to include short gaps between accepted chain
seeds as part of the supported span. No gap, drift, phase, recurrence or weight
threshold changes. Long interruptions still split chains. This directly tests
whether exact-seed-only counting causes the remaining indel undercount; it can
also overcount short non-repeat insertions, which remains a limitation.
All old outcomes are preserved. There will be no third revision.

## Development decision criteria and withdrawn FPR gate

For a method-innovation claim, A3 needed to demonstrate a material benefit beyond
ordinary mapping, not merely the diagnostic-k-mer failure. The development
screen retained >=15% relative MARE reduction versus matched A0 **and A1**,
>=50% decrease of background-homology excess abundance versus A0, better indel
recall than A2, clean MARE worsening <=0.02 versus A0, and reporting of all
families including ambiguity/missingness. Runtime <=2x A0 is the preferred cost
boundary. The historical `binary_FPR_at_most_002` development gate is withdrawn:
the synthetic F2-role rate does not measure the same event or denominator as
the historical classifier FPR. Its raw receipt/code/value remain retained, but
it is not a pass/fail criterion. These development criteria do not establish
validation. If the final permitted revision does not show sufficient overall
advantage beyond A1, do not generate heldout data or tune against it.

Only if this screen passes may fresh source families, additional mixed/read-
error conditions and a one-time heldout execution be preregistered. Absolute
MARE <=0.30 on a different test distribution cannot substitute for the matched
comparisons above. SRF's existing development evidence stays labelled as such;
formal unified-endpoint comparison is deferred until the algorithm decision.

## Efficiency budget

The user sets approximately 30 minutes for this round. Finish started 100-Mb
tests with byte-identical outputs, wall/CPU time and RSS. A single 1-Gbase pair
may run only if >=2x speedup persists and the remaining budget accommodates it.
Do not launch three pairs now; if scaling drops below 1.5x or outputs differ,
stop expansion. GB-scale and paper-level performance remain unverified when
the budget does not permit the corresponding run.
