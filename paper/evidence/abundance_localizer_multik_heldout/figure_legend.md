# Figure 5 legend

**Divergence-aware anchor bridging passes predeclared localization gates.**
(A) Development and held-out design. Development v1 and v2 use the same consumed
seeds; held-out seeds 5501--5503 were first used after the v2 rule, evidence
hashes and guard passed hosted CI. (B) Full-assembly mean base recall by unit
divergence, planted segment count and dataset. The dotted line is the aggregate
0.95 development/held-out gate. (C) Mean absolute difference between predicted
and planted fragment counts on full assemblies; the symmetric-log axis retains
zero while showing the large v1 errors. (D) Held-out base recall across positive
assembly retention levels, averaged across segment counts and repeat families.
(E) Held-out mean absolute relative error in localized repeat bp across the same
levels. (F) Sensitivity, false-positive rate and precision for the single-k21
baseline and the previously frozen multi-k/depth rule. The classifier gains
sensitivity while increasing false positives and lowering precision.

All localization panels use union-scored genomic bases. The experiment supplies
the planted founder catalogue and models independent length-preserving
substitutions plus 500-bp interruptions. It does not model unknown-catalogue
discovery, indels, chromosome-scale assembly fragmentation, empirical satellite
evolution or biological assembly collapse.
