# Held-out conditional assembly under-representation evidence

This compact archive records the first and only use of the predeclared held-out
seeds 5101–5103 from `benchmarks/configs/abundance_v1.json`. The unchanged
configuration supplies the planted monomer catalogue and crosses 1/5/20× read
coverage, 0/0.1/1% substitutions and assemblies retaining 100/75/50/25/0% of
each array. It isolates quantification, exact-copy localization and the fixed
0.6 assembly/read-ratio decision rule; it does not test de novo discovery,
biological sequence divergence, empirical HiFi errors or real assembly collapse.

All 177 commands completed under the source frozen at Git `f16596b`. The
environment warning is expected because the executed native extension is a
generated, untracked build; its hash and the complete source manifest identify
the executed code. The archiver rechecked the full source snapshot, benchmark
helpers, configuration, matrix dimensions and all 177 successful receipts.
`archive_manifest.json` verifies the 27 copied/generated evidence files.

Across 243 positive and 162 control family conditions, outcomes were TP=208,
FN=35, FP=8 and TN=154: sensitivity 0.855967, false-positive rate 0.049383 and
precision 0.962963. Every one of the 36 positive assembly/family localization
conditions had base recall 1.0; minimum base precision was 0.996732. This simple
localization result must not be generalized to divergent or fragmented arrays.

The aggregate conceals an important failure. At 20× with 1% substitutions,
50%-retained arrays were detected in only 2/9 family conditions. At 1×, the
complete-assembly control produced one false call in each error tier. These
patterns are consistent with error-dependent read copy-number underestimation
and finite-sampling variation, while they do not isolate a unique causal
mechanism. The held-out seeds are now consumed and must not be used for tuning.
Any revised decision model must be developed on development data and evaluated
once on newly predeclared independent seeds.

`copy_number_metrics.tsv`, `localization_metrics.tsv` and
`comparison_metrics.tsv` retain every family condition. Summary files keep the
requested strata. `resource_metrics.tsv` is generated from the 177 verified
receipts and records command JSON, wall time, direct-child RSS and CPU counters.
The full 263-MB run remains at the `source_snapshot` and receipt paths recorded
in `environment.json`; it is not required to read the compact scientific tables.
