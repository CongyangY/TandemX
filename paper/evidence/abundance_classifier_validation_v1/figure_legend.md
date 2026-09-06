# Figure 6 legend

**A seed-robust abundance classifier gains sensitivity but fails independent
false-positive and precision gates.** (A) Evaluation sequence. Pooled
development v1 failed parameter-stability gates; seed-robust development v2
passed, after which code, parameters and evidence hashes were frozen before
one-time evaluation of seeds 5701--5703. (B) Overall sensitivity, false-positive
rate and precision deltas from the single-k21 baseline in development v2 and
held-out data. Favorable directions are positive for sensitivity/precision and
nonpositive for false-positive rate. (C) Absolute held-out metrics for the
single-k21 baseline and the frozen alpha-0.5/threshold-0.5 blend. (D) Held-out
metric deltas by seed. Seed 5703 shows the largest false-positive and precision
regressions. (E) Sensitivity and false-positive-rate deltas by nominal read
coverage. All additional false positives occur at 1x; this stratification is a
post-hoc diagnosis and not a new validation. (F) Full-assembly localization
recall by unit divergence and planted segment count. The dotted line is the
predeclared aggregate mean-recall gate; the annotation reports precision across
positive assemblies and false-positive rows among absent families.

Classifier denominators comprise 1,458 positive and 972 negative family
conditions. Localization uses union-scored genomic bases. The known-catalogue
simulation models independent length-preserving substitutions and 500-bp
interruptions; it does not model unknown-catalogue discovery, indels, empirical
satellite evolution or biological assembly collapse.
