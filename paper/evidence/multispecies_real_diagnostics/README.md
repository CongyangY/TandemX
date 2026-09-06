# Additional whole-library real-input diagnostics

Committed native code96fd5dc; fixed whole-file seed6101 subsets, after complete
source QC. Same FASTA, period30–1000 bp, span100 bp and one thread for all three
tools. All runs listed here exited successfully and normalized without error.

| Material / observed input | TandemX seconds / MiB | TRF seconds / MiB | TideHunter seconds / MiB |
| --- | ---: | ---: | ---: |
| Morex, 12,542,086 bp / 572 reads | 22.143 / 127.81 | 27.316 / 171.00 | 6.456 / 257.64 |
| Morex, 115,272,341 bp / 5294 reads | 159.600 / 173.03 | 274.546 / 179.00 | 60.409 / 489.81 |
| Nipponbare, 11,418,016 bp / 626 reads | 10.910 / 106.53 | 16.629 / 201.17 | 4.372 / 168.23 |
| Nipponbare, 113,627,173 bp / 6141 reads | 104.784 / 212.75 | 140.957 / 333.86 | 46.330 / 450.64 |

Call counts and union bases are in each original summary.tsv. More calls or more
covered bases are not accuracy or recall without independent truth. TandemX is
faster than TRF and slower than TideHunter in these observations; its measured
peak RSS is lower than both here. Concurrent downloads, QC and other jobs mean
these are development diagnostics, not final isolated performance rankings.
No all-metric superiority is claimed. Larger subsets are nested technical
samples, not additional plants or independent biological replicates.

Each folder retains source/parameter/tool hashes, raw execution resources,
logs and output hashes. Full FASTA, native predictions, normalized intervals and
source snapshots remain at the recorded T7 paths. Public ENA source manifests,
download/QC/sampling commands and comparator source recipes permit regeneration.
