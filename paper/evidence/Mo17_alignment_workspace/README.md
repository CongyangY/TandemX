# Exact rolling-workspace development profiles

Actual111,505,681-bp Mo17 input,8084 reads and39,438 native alignment calls.
All seven complete discovery products are byte-identical to the original
indexed-clustering baseline. No detection, alignment score, tie or clustering
threshold changed. Exact source snapshots/hashes are included because these
profiles ran during development on top of99c0c66, before the changes were committed.

| Workspace | Native alignment seconds | Full child seconds | Peak RSS MiB |
| --- | ---: | ---: | ---: |
| Original per-row allocations | 123.8034 | 151.2617 | 174.7500 |
| Two reused row pairs | 104.1338 | 139.0149 | 173.3594 |
| In-place score/peak pair | 93.3743 | 124.8605 | 184.2188 |

The original profile is archived in
`../Mo17_native_index/full_pipeline_111Mb_profile`. The final native stage is
24.58% shorter and the full profiled child17.45% shorter than that run, while
observed peak RSS is5.42% higher. These are single development profiles under
concurrent downloads/QC/container jobs, not isolated resource comparisons or
external-tool superiority. RSS fluctuation does not establish a memory gain.
Both intermediate and final profiles are retained, including unfavorable values.
Binary/profile hashes are recorded; the binaries and `.prof` files stay on T7.
