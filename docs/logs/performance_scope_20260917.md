# Performance evidence and its scope, 2026-09-17

This ledger separates measurements of the **production discovery workflow**
from research prototypes, competitor resource gates, and unmeasured scaling.
It is an expert-facing index to prior receipts, not a new benchmark or a
cross-platform rank. Historical runs below have not all been rerun in this
session; the cited evidence archives and dated `docs/current_status.md`
checkpoint retain their original commands, hashes, and conditions.

| Workload and comparison | TandemX wall / peak RSS | Other measured wall / peak RSS | What the observation supports |
| --- | ---: | ---: | --- |
| Mo17 11.681 Mb exact-output alignment optimization | 17.348 s / 81.703 MiB before; 12.361 s / 62.484 MiB after | same TandemX implementation before edit | One paired replay, 28.75% less wall and 23.52% less RSS; six core outputs byte-identical. Not a replicate distribution. |
| Mo17 111.506 Mb exact-output alignment optimization | 134.684 s / 178.438 MiB before; 100.039 s / 154.766 MiB after | same TandemX implementation before edit | One paired replay, 25.72% less wall and 13.27% less RSS; seven outputs byte-identical. |
| Chinese Spring 126.731 Mb real-read comparison | 139.529 s / 223.828 MiB | TRF 332.287 s / 311.938 MiB; TideHunter 48.230 s / 413.500 MiB | Descriptive same-input resource values, without a matched real-read discovery truth set. TandemX is faster than TRF and slower than TideHunter in this run. |
| Victoria oat 110.203 Mb real-read comparison | 164.738 s / 250.45 MiB | TRF 362.689 s / 176.72 MiB; TideHunter 81.982 s / 621.45 MiB | TandemX is faster than TRF, slower than TideHunter, and uses more RSS than TRF here. No universal memory advantage. |
| Morex 1.169928 Gb frozen three-tool diagnostic | 1,909.525 s / 913.891 MiB | TRF 2,409.145 s / 286.172 MiB; TideHunter 566.196 s / 580.375 MiB | Large-input completion. TandemX was 20.74% faster than TRF, 3.37 times TideHunter's wall time, and 3.19 times TRF's RSS. Jobs overlapped, so this is not an isolated final ranking. |
| Mo17 1.129 Gb complete discovery replay | 1,241.876 s / 697.109 MiB | Previous TandemX replay 1,625.652 s; reported RSS rose 0.81% | Byte-identical seven products and 23.61% less wall for that paired engineering replay. Different from a competitor comparison. |

Sources: `docs/current_status.md` sections “Exact-output elastic alignment
optimization”, “Published article draft and resolved figure export checkpoint”,
“Completed interval calibration and reference-concordance QC”, and “Joint-read
uncertainty, complete TRASH evaluation and seven-species file QC”. The compact
Mo17 parity receipts reside under `paper/evidence/discovery_packed_trace_batch_v1`;
other exact source receipts are identified in those dated checkpoints and, for
large datasets, on T7. The numeric values are historical archived observations,
not fresh measurements made by this log.

The fresh 2026-09-17 Col-CEN whole-reference audit used minimap2 on only six
selected reads, with minimap2 reporting 1.795 s real and 0.856 GB peak RSS
while indexing a 132-Mb reference. That is **source anchoring cost**, not a
TandemX runtime, discovery throughput, or assembly-audit scaling point. The
M1 full-read research timer, chunked mapping comparator timer, and production
timer cover different work and output stages; their raw wall/RSS values are
archived but deliberately not ranked. The Ey15 native span diagnostic uses
seven selected reads and nine small contexts and supplies no meaningful
large-genome speed or memory measurement.

TideCluster's Morex 100-Mb stage measurements were 331.52 s / 5,152,444 kB
for its TideHunter stage and 60.99 s / 7,770,936 kB for clustering. The
frozen 1-Gb preflight refused execution because the 100-Mb clustering peak
nearly filled the same 8.217-GB Docker limit. This is a **host/container
resource gate**; it is not evidence that TideCluster cannot run at 1 Gb on a
machine with more memory. The separate small HOR-control profiler sums live
process-tree RSS and can double-count shared pages, so it cannot be compared
to the direct-child or GNU-time peaks above.

The current evidence does **not** establish a universal speed or memory
advantage, multithread scaling efficiency, end-to-end performance on an
isolated several-Gb wheat/rye genome, or accuracy per unit of compute. Current
T7 access is through a 480-Mb/s USB2 path with prior I/O failures; the full
large-input rebenchmark must wait for a verified reliable storage path and a
method/endpoint freeze. Small development runs and successful CI do not close
that gate.
