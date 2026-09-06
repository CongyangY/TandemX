# Additional whole-library real-input diagnostics

Committed source 96fd5dc for smaller Morex/Nipponbare runs, d2aa4d6 for
Victoria, and 438d078 for Chinese Spring/Lo7 plus the 1.170-Gb Morex run;
fixed whole-file seed6101 subsets, after complete
source QC. Same FASTA, period 30–1000 bp, span 100 bp and one thread for all three
tools. All runs listed here exited successfully and normalized without error.

| Material / observed input | TandemX seconds / MiB | TRF seconds / MiB | TideHunter seconds / MiB |
| --- | ---: | ---: | ---: |
| Morex, 12,542,086 bp / 572 reads | 22.143 / 127.81 | 27.316 / 171.00 | 6.456 / 257.64 |
| Morex, 115,272,341 bp / 5294 reads | 159.600 / 173.03 | 274.546 / 179.00 | 60.409 / 489.81 |
| Morex, 1,169,928,427 bp / 53,715 reads | 1909.525 / 913.89 | 2409.145 / 286.17 | 566.196 / 580.38 |
| Nipponbare, 11,418,016 bp / 626 reads | 10.910 / 106.53 | 16.629 / 201.17 | 4.372 / 168.23 |
| Nipponbare, 113,627,173 bp / 6141 reads | 104.784 / 212.75 | 140.957 / 333.86 | 46.330 / 450.64 |
| Victoria, 11,361,028 bp / 615 reads | 16.827 / 122.97 | 44.997 / 165.78 | 11.097 / 257.20 |
| Victoria, 110,202,741 bp / 5972 reads | 164.738 / 250.45 | 362.689 / 176.72 | 81.982 / 621.45 |
| Chinese Spring, 12,698,239 bp / 765 reads | 13.410 / 97.08 | 32.879 / 196.31 | 4.537 / 150.95 |
| Chinese Spring, 126,730,796 bp / 7633 reads | 139.529 / 223.83 | 332.287 / 311.94 | 48.230 / 413.50 |
| Lo7, 11,639,976 bp / 693 reads | 13.159 / 106.64 | 32.334 / 167.05 | 4.498 / 277.14 |

Call counts and union bases are in each original summary.tsv. More calls or more
covered bases are not accuracy or recall without independent truth. TandemX is
faster than TRF and slower than TideHunter in these observations. Peak RSS is
input-dependent: TandemX exceeds TRF by 41.72% in the 110-Mb Victoria run, and
the 1.170-Gb Morex run uses 913.89 MiB versus 286.17 MiB for TRF and 580.38 MiB
for TideHunter. In that Morex run TandemX is 20.74% faster than TRF but takes
3.37-fold the TideHunter wall time. Concurrent downloads, QC and other jobs mean
these are development diagnostics, not final isolated performance rankings.
No all-metric superiority is claimed. Larger subsets are nested technical
samples, not additional plants or independent biological replicates.

Each folder retains source/parameter/tool hashes, raw execution resources,
logs and output hashes. Full FASTA, native predictions, normalized intervals and
source snapshots remain at the recorded T7 paths. Public ENA source manifests,
download/QC/sampling commands and comparator source recipes permit regeneration.

## Cross-run diagnostic figure contract

The six-panel renderer `benchmarks/scripts/plot_real_comparator_diagnostics.py`
compares wall time, throughput, peak RSS, paired wall-time and RSS ratios to TRF,
and called-base fraction across completed compact archives. Its analytical
question is whether the present one-thread observations reveal stable scaling
and resource trade-offs that should guide optimization. The supported takeaway
is descriptive: ranking varies by endpoint and input, and current TandemX does
not dominate TideHunter in speed. Called-base fraction measures output extent,
not precision or recall.

Every point is linked to a hash-validated run archive and retained in
`panel_source.tsv`. Color identifies tools and marker shape identifies material;
nested sample sizes are connected only within the same material/tool. The
PDF/SVG/PNG must preserve editable SVG text and vector marks. Because acquisition
overlapped these single executions, this is a supplementary development
diagnostic rather than the final isolated repeated performance figure.

`figures_v1/` retains the first render because its tool/material legends overlap
the top panel titles. It is rejected for manuscript use. `figures_v2/` adds a
dedicated legend row and passed visual inspection. Version 2 has 78 editable SVG
text nodes, no raster image node, matching output/provenance hashes and a panel
source TSV byte-identical to version 1. Version 2 is the manuscript figure.
