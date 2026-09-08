# TideCluster factorial-validation figure contract

Initial status: frozen before any TideCluster output from
`tidecluster_factorial_validation_v1` was inspected. That runner stopped after
the first external-process failure. The v2 continuation and missing-value
rendering rules below were frozen after that failure but before any previously
unattempted seed-setting cell was executed.

## Claim boundary

The figure evaluates pinned TideCluster 1.21.2 on three deterministic 10-Mb
factorial assemblies with planted truth. Seeds 6401--6403 were previously used
for TandemX quantification validation, but both TideCluster settings were frozen
before TideCluster output inspection. These are technical comparisons on
same-process simulations, not biological replicates or a whole-genome plant
benchmark.

The unchanged TideCluster default is primary. The second setting changes only
the accepted period range from 40--3000 bp to the benchmark-matched 30--1000 bp;
it is a declared sensitivity setting, not result-dependent tuning.

## Panels

- **A, Array detection:** recall and precision for every seed and setting.
- **B, Base-union agreement:** base-union recall and precision for every seed
  and setting.
- **C, Boundary and period error:** matched boundary MAE and period MAE for
  every seed and setting; unmatched cases remain unavailable, not zero.
- **D, Monomer recovery:** cyclic-Levenshtein monomer recall and homologous
  distinct-consensus fraction for every seed and setting.
- **E, Stage elapsed time:** TideHunter and clustering wall time for every seed
  and setting on a log scale.
- **F, Stage memory:** TideHunter and clustering maximum RSS for every seed and
  setting on a log scale.

## Visual and source-data rules

- Preserve all six frozen cell fates; do not omit an unfavorable seed or
  setting. Accuracy panels show only independently verified successful cells
  and state the successful denominator out of six.
- Failed external stages are labelled as failed/missing and are not plotted as
  zero-valued accuracy. Their measured failed-stage time/RSS may be shown with
  an `x`, while dependent unstarted stages remain unavailable.
- Use a restrained blue/orange palette for the two settings, with point shapes
  and stage line styles as redundant non-colour encodings.
- Export SVG, PDF and PNG plus a panel-level source TSV and provenance JSON.
- SVG text must remain editable and the SVG must contain no raster image nodes.
- Independently verify all plotted values before rendering, inspect the direct
  PNG and a separately rendered PDF page, and reject clipping or overlap.
