# Figure 6 TXF000695 QC source notes

This table is a descriptive manuscript-layer summary of the complete diagnostic
k-mer counts used in the frozen Macadamia Illumina analysis. It does not change
the estimator, family eligibility or final interpretation.

Inputs:

- `/Volumes/T7/Codex/TandemX/validation/orthogonal_abundance_v1/kmc_targets/macadamia_k21.tsv`
- `/Volumes/T7/Codex/TandemX/validation/orthogonal_abundance_v1/kmc_targets/macadamia_k31.tsv`
- `/Volumes/T7/Codex/TandemX/validation/orthogonal_abundance_v1/kmc/macadamia_k21_cs1e9_target_read_counts.tsv`
- `/Volumes/T7/Codex/TandemX/validation/orthogonal_abundance_v1/kmc/macadamia_k31_cs1e9_target_read_counts.tsv`
- `paper/evidence/orthogonal_abundance_validation_v1/results/macadamia_k21_family_metrics.tsv`
- `paper/evidence/orthogonal_abundance_validation_v1/results/macadamia_k31_family_metrics.tsv`

Counts absent from the KMC `-ci2` intersection were assigned 0.5 for the median
and distribution summaries, matching the frozen zero-or-one midpoint rule. The
cumulative columns `depth_below_100_count`, `depth_below_1000_count` and
`depth_below_10000_count` include all lower bins. The high-depth column covers
counts from 50,000 through 99,999; the observed maximum was below 100,000 at
both k values.

Derived cross-k values:

- median-depth fold, k21/k31: 763.2772277227723;
- abundance fold, k21/k31: 670.1605306966903.

The distribution is bimodal at both k values. The number of words in the
50,000–99,999 component decreases from 259 at k=21 to 202 at k=31, moving the
median from the high-depth to the low-depth component. This is reported as
median-boundary sensitivity, not as an independently validated classifier or a
resolved biological mechanism.
