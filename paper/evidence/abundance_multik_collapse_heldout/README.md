# Frozen multi-k assembly-comparison validation

This archive applies the rule calibrated only on development seeds 4101–4103
to fresh held-out seeds 5201–5203. The rule and calibration hashes were committed
in `1231743` before either held-out command ran. Held-out mode performed no
threshold selection: it used k=15/21/27/31 log-linear copy estimates when
available, otherwise k=21, with a ratio threshold of 0.5 below observed haploid
depth 2 and 0.6 otherwise.

Across all 405 conditions, baseline TP/FN/FP/TN=198/45/14/148 and frozen-rule
TP/FN/FP/TN=208/35/12/150. Sensitivity improved 0.814815 to 0.855967,
false-positive rate 0.086420 to 0.074074 and precision 0.933962 to 0.945455.
Every held-out seed gained true positives without gaining false positives, but
seed 5202 retained 12 false positives and a 0.222222 false-positive rate. The
20x/1%-substitution/50%-retention stratum improved from 3/9 to 9/9; 1x complete-
assembly false calls decreased from five to three.

The aggregate gain is not uniform. At 1x with 0% or 0.1% substitutions and 50%
retention, sensitivity decreased from 6/9 to 4/9 in each stratum. The source
heatmap retains these adverse cells; future work must improve the decision rule
without tuning on the consumed held-out genomes.

The multi-k-only arm had 15 unavailable rows and is not a complete-method
comparison. These exact-copy 199.1-kb simulations do not establish performance
on divergent/fragmented arrays, unknown catalogues, empirical HiFi errors or
biological under-representation. `archive_manifest.json` covers 19 files and
rechecks the frozen calibration hashes, baseline receipts, paired metric rows,
method confusion counts and selected source snapshot. Calibration files and
baseline receipts are copied into the archive for durable review.

`figures_v2/` is the inspected six-panel manuscript render. It shows the
predeclared design, overall metrics, confusion counts, per-seed sensitivity and
false-positive rates, and all nine 50%-retention coverage/error cells. The
source-backed heatmap therefore displays both gains and the two adverse 1x
strata. Its SVG contains 111 editable text nodes and no raster image node; the
single-page PDF was rendered to PNG and visually inspected without clipping or
overlap. `figures_v1/` is retained as a rejected layout because the panel-B
legend obscured low-value labels.
