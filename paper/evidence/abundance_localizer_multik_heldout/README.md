# Divergence-aware localizer and frozen multi-k held-out validation

This archive applies the previously frozen exact-copy multi-k/depth rule to the
completed 5501--5503 held-out baseline. It reused all reads and localization
outputs, evaluated 162 multi-k read conditions and did not rerun the 1,062
baseline commands or fit a threshold on held-out data. All 7,290 method rows are
present; unlike the earlier domain-shift validation, no multi-k estimate was
unavailable in this matrix.

The single-k baseline had TP/FN/FP/TN=875/583/12/960: sensitivity 0.600137,
false-positive rate 0.012346 and precision 0.986471. The frozen hybrid had
1339/119/83/889: sensitivity 0.918381, false-positive rate 0.085391 and precision
0.941632. The sensitivity gain of 0.318244 therefore accompanied 71 additional
false positives and a precision decrease of 0.044840. Multi-k alone reached
sensitivity 0.940329, false-positive rate 0.128601 and precision 0.916444. These
results do not support universal classifier superiority.

`archive_manifest.json` verifies 19 files and rechecks frozen calibration,
baseline hashes, source snapshots, complete paired keys and confusion counts.
`figures_v3/` is the accepted six-panel render. Its SVG has 113 editable text
nodes and no raster image node; PNG and rendered PDF were visually inspected.
`figures_v1/` and `figures_v2/` are retained as rejected layouts because panel-A
text columns overlapped. Every version keeps its panel source and output hashes.

The successful localization result is limited to a known-catalogue simulation
with independent length-preserving substitutions. Exact-k-mer IID identity is
not alignment identity, and neither the planted ratios nor the classifier calls
constitute biological evidence of assembly collapse.
