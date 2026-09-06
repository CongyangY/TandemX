# Conditional copy-number evaluation on three 10-Mb source genomes

Development diagnostic, not a final comparative publication figure. Each source
genome has 54 factorial founder arrays (three periods, copy counts and target-GC
fractions crossed with two unit-divergence rates) and one 1.026-Mb array. Read
lengths follow the complete Mo 17 library's histogram; only lengths are empirical.
Three source coverages and three independent-error models yield 27 paired read
conditions. High error means substitution 1%, insertion 0.5%, deletion 0.5%; low
error means 0.1% each. Biological unit substitutions are independent of read errors.

**A**, Mean signed relative error at 20x for the 54 ordinary factorial families.
Each divergence/error group contains 81 family conditions (27 per source genome).
Large connected markers are means; small markers show the three genome-specific
means. Grey native median k 21 and orange mean k 21 frequently overlap. Dashed
line denotes zero bias; points are not independent biological samples.

**B**, Mean absolute relative error for 2%-divergent families in high-error reads.
All three methods are restricted to the same available pairs:73/81 at 1x and 81/81
at 5x/20x. Unavailable fits are not assigned zero error; all-input outcomes are in C.

**C**, Per-family comparison with native median k 21 over all 495 conditions at each
coverage, including the megabase arrays and all error settings. Lower/higher error
uses a 1e-12 absolute relative-error difference tolerance. Zero-support fits remain
unavailable (63 at 1x; three at 5x; none at 20x). Two positive-slope fits remain in
the available categories with model-violation flags in the raw table. Counts,
including categories too small to label on the plot, are in `panel_source.tsv`.

**D**, Empirical truth coverage of the original diagnostic-k-mer 10th–90th percentile
spread for 2%-divergent ordinary arrays. Each cell contains 81 conditions. These
bounds are not a nominal sampling confidence interval and are not recalibrated
by relabelling them. Across the complete baseline,77/1,485 spreads contain truth.

**E**, Estimated/planted copy ratio for the separate 171-bp×6,000 arrays (1% unit
divergence) in high-error reads. Small points represent three source genomes;
connected points are their means. The dashed line is exact planted copy count.

**F**, At 20x, maximum absolute log-fit residual versus absolute estimator-minus-
source-sampling-oracle error, normalized by planted copies. All 495 fitted family
conditions are shown. The oracle uses sampled source repeat bases and actual
source coverage; it is evaluation-only. Colour indicates biological divergence;
1% denotes the megabase arrays. The residual axis is symmetric-log with a 1e-4
linear threshold to display zeros. Fit residuals are not sampling confidence.

All measured rows, scripts, provenance and model failures are retained. No formal
independence-based significance, real-data calibration or external superiority
is inferred from this controlled three-genome development experiment.
