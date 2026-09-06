# TRASH2 default assembly benchmark: one development genome

Pinned source f290a5e and immutable image ID are in execution.json. Actual input
is the independently generated 10-Mb s6301 genome, containing 55 families and
2,201,400 planted repeat bases, including a 1.026-Mb array. No truth/templates
were passed to TRASH2. Execution succeeded in 844.52 seconds; GNU time peak RSS
was 713.29 MiB and cgroup memory peak 1,528,770,560 bytes including cache.
This run predates whole-cgroup CPU counters. Its GNU time user/system totals
do not represent all R worker CPU. Concurrent-job timings are diagnostics.

| Endpoint | Observed result |
| --- | ---: |
| Cyclic sequence family recovery | 55/55 |
| Native approximate-array IoU + period matches | 49/55 |
| Native unit-membership-derived array matches | 55/55 |
| Unit-derived matched period MAE | 0 bp |
| Unit-derived matched boundary MAE | 3.1545 bp |
| Actual unit base-union precision | 0.9998987 |
| Actual unit base-union recall | 0.9999437 |
| Exact source-sequence coordinate agreement | 11,436/11,436 units |

The author README explicitly describes array-table boundaries as approximate.
`evaluation_v1` retains that original endpoint. `evaluation_v2` adds the main
repeat table's explicit array memberships and outer unit boundaries, while
retaining all original results. All 55 arrays have observed units and every
unit ID joins unambiguously; this normalization uses no planted truth. Native
peak and consensus-length policies give the same results here.

The six approximate-window IoU failures are not six missed repeat families.
This is a strong default comparator on this simple development genome. No
TandemX advantage is inferred from comparisons to a different HiFi input, from
coarse-window boundaries or from one synthetic seed. Additional independent
genomes, harder backgrounds, input scales and real assemblies remain required.
Original native products, logs, source/parameters and both evaluations are
preserved with checksums. Table/interval uncertainty is descriptive, because
arrays within a genome are not independent biological samples.
