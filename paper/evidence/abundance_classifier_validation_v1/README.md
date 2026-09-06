# Frozen abundance-classifier development and held-out validation

This compact archive preserves the complete evidence chain for the transparent
single/multi-k classifier experiment. Development v1 used seeds 5601--5603 and
failed its predeclared leave-one-seed-out stability gate. A seed-robust v2
selector reused only those development rows, selected alpha 0.5 and decision
threshold 0.5, and passed its development gates. Commit `300e48d` and both
hosted Ubuntu/macOS checks preceded the first use of held-out seeds 5701--5703.

The frozen blend did **not** pass held-out validation. Across 2,430 paired family
conditions, single k=21 produced TP/FN/FP/TN=916/542/29/943 (sensitivity
0.628258, false-positive rate 0.029835 and precision 0.969312). The frozen blend
produced 1044/414/45/927 (0.716049, 0.046296 and 0.958678). Sensitivity increased
by 0.087791, but false-positive rate increased by 0.016461 and precision decreased
by 0.010634. Seed 5703 had the largest adverse changes: false-positive rate
increased by 0.040123 and precision decreased by 0.025985.

All 16 additional false positives occurred in the nominal 1x stratum. At 5x and
20x, both methods had zero false positives while the blend increased sensitivity
by 0.094650 and 0.082305, respectively. This coverage split is a post-hoc
diagnosis on consumed held-out data. It cannot convert the failed model into a
validated one or support refitting and retesting on seeds 5701--5703.

The same held-out baseline independently retained useful localization behavior:
full-assembly mean base recall was 0.973858, positive-assembly mean precision was
0.999631, and no predicted bases occurred in 54 absent-family rows. The weakest
full-assembly stratum was 5% unit divergence with three array segments (mean
recall 0.937682). Thus the classifier failure is retained separately from the
aggregate localizer result.

`archive_manifest.json` covers the copied development, baseline, paired multi-k
and held-out files plus their frozen benchmark-source snapshots. The archive
script verifies every chain hash, exact method pairing, seed separation and
reported confusion count. The baseline's 1,062 command receipts are retained in
`baseline/resource_metrics.tsv`; their summed direct-child wall times were
19.434 s for 90 localization, 26.848 s for 162 quantification and 95.292 s for
810 comparison commands. Maximum direct-child peak RSS was 32.719 MiB. These
small synthetic command measurements are reproducibility diagnostics, not
production performance or external-tool comparisons.

`figures_v2/` is the accepted six-panel render. Its SVG contains 123 editable
text nodes and no raster image node. The PNG and PDF rendering were inspected.
`figures_v1/` is retained as a rejected layout because the panel-A arrow crossed
the stage labels; its numerical panel source is unchanged.

The benchmark supplies a known planted catalogue and uses synthetic IID,
length-preserving substitutions. It is neither de novo family recovery nor
biological evidence of assembly collapse, and it does not establish universal
TandemX superiority.
