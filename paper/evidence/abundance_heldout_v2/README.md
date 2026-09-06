# Fresh predeclared assembly-comparison baseline

This compact archive records the first and only execution of held-out seeds
5201–5203 from committed config `benchmarks/configs/abundance_v2.json`. Source
commit `1231743` was pushed and both hosted CI runs passed before execution.
All 177 commands succeeded, yielding 81 copy-number, 45 localization and 405
comparison-family rows.

The original k=21 rule produced TP/FN/FP/TN=198/45/14/148 across 243 positive
and 162 control conditions: sensitivity 0.814815, false-positive rate 0.086420
and precision 0.933962. At 20x with 1% substitutions and 50% retained copies,
sensitivity was 3/9. At 1x, fully retained assemblies produced five false calls
among 27 family conditions. These failures remain in the archive.

The genomes are only 199.1 kb and contain supplied, exact-copy monomers. They
isolate conditional quantification/localization behavior and are not discovery,
divergent-array or biological validation. `archive_manifest.json` covers 28
copied/generated evidence files and `resource_metrics.tsv` retains all 177
validated execution receipts. Seeds 5201–5203 are consumed.
