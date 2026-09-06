# Cross-tool source-query recovery on real plant reads

Actual consensus outputs from completed Mo17 111-Mb, Morex 115-Mb and
Nipponbare 114-Mb comparisons were scored against one source-backed query per
material. TandemX final family representatives and pre-clustering candidates
were kept as separate stages. TRF and TideHunter rows are per-array consensuses,
so their consensus counts and homologous fractions are not family precision.

The independent evaluator canonicalizes rotations/reverse complements and uses
global cyclic Levenshtein similarity. At the declared 0.90 threshold:

| Material / historical query | TandemX family | TandemX candidate | TRF | TideHunter |
| --- | ---: | ---: | ---: | ---: |
| Mo17 / CentC AF078922.1 | recovered, best 0.9872 | recovered, best 0.9872 | recovered, best 0.9872 | recovered, best 0.9872 |
| Nipponbare / Rice358 X55642.1 | recovered, best 0.9358 | recovered, best 0.9358 | recovered, best 0.9581 | recovered, best 0.9581 |
| Morex / HvT01 X16095.1:1–118 | not recovered | recovered, best 0.9068 | recovered, best 0.9068 | recovered, best 0.9068 |

The Morex boundary was traced to read `ERR4659246.235349`. TandemX detected
three threshold-supported 118-bp candidates; one belongs to `TXF000053`, whose
fixed observed representative has similarity 0.8983 to the historical query.
At threshold 0.89 all four output/stage categories recover the query, whereas
at 0.95 none does. This is evidence for threshold and family-representation
sensitivity, not justification to tune against this one query.

The query donors differ from the tested materials: CentC is from Seneca 60,
HvT01 from Donetsky A., and Rice358 from Cigalon; the latter also lacks validated
native monomer boundaries. Consequently these are selected source-query checks,
not complete known-family recall, false-negative rates or proof that the queries
occur unchanged in the tested donors. The exact output/input/script hashes and
stage examples are retained in each receipt. T7 versions 1 and 2 preserve the
initial evaluator and provenance-field development; version 3 is authoritative.

`archive_manifest.tsv` records byte-identical copies from
`/Volumes/T7/Codex/TandemX/results/known_query_recovery_v3_20260906`.
