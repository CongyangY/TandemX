# Figure 7 legend

**Figure 7. Fresh held-out validation of a frozen depth-gated abundance
classifier.** **A,** Evidence sequence. The adverse 5701-5703 result was consumed
as development evidence; v3 kept single k=21 below estimated haploid depth 2 and
used the alpha-0.5/threshold-0.5 blend otherwise. Commit `62892a6` and hosted
Ubuntu/macOS CI preceded one-time evaluation of seeds 5801-5803 without fitting
or model selection. **B,** Sensitivity, false-positive-rate and precision deltas
relative to single k=21 for the six-genome development set and the fresh
three-genome held-out set. **C,** Held-out TP, FN, FP and TN counts across 2,430
paired family-conditions. **D,** The same metric deltas for each held-out seed.
Every seed met the predeclared sensitivity, FPR and precision guardrails.
**E,** Metric deltas by nominal read coverage. The rule is identical to baseline
for all 810 estimated-depth-below-2 rows, so the 1x stratum is unchanged; gains
occur at 5x and 20x without added false positives. **F,** Equal-weight mean
full-assembly base recall across nine family-condition rows per planted monomer
substitution/segment stratum. The dotted 0.95 line is the aggregate target shown
as a reference, not a requirement for every stratum. Across all positive
assemblies, mean base precision was 0.999617; no predicted bases occurred in 54
absent-family rows.

Sensitivity is TP/(TP+FN), false-positive rate is FP/(FP+TN), and precision is
TP/(TP+FP). The matrix contains three monomer lengths, three copy numbers, three
unit-divergence levels, one or three planted array segments, three nominal read
coverages, three read substitution rates and five assembly-retention levels.
The classifier endpoint uses a known monomer catalogue. No confidence interval
is drawn because only three new genome seeds were tested; all raw conditions and
per-seed counts are provided. This IID substitution simulation is not biological
collapse validation and does not support a universal performance claim.

Accepted files are in `figures_v2/`. Exact panel values are in
`figures_v2/panel_source.tsv`; input/output hashes are in
`figures_v2/figure_provenance.json`.
