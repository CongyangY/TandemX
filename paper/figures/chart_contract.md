# Figure contracts

Publication graphics use Matplotlib, editable SVG text and embedded PDF fonts.
No AI-generated quantitative panels. Raw data and executable plotting code are
the authority; these specifications do not imply that a planned panel exists.

## Figure 6: frozen abundance-classifier validation

- Analytical question: does the seed-robust single/multi-k blend transfer to
  untouched simulated genomes while preserving false-positive rate and precision?
- Takeaway: held-out sensitivity increases by 0.087791, but false-positive rate
  increases by 0.016461 and precision decreases by 0.010634; the frozen model
  fails its predeclared gate, with all additional false positives at nominal 1x.
- Data sufficiency: 2,430 exactly paired family conditions from three held-out
  genomes, plus 270 localization rows; panels retain aggregate, seed, coverage
  and divergence/fragmentation views.
- Forms: process/status diagram; grouped bars; signed delta bars; coverage line
  comparison; two-series localization line chart.
- Palette: grey baseline, blue selected/development, orange adverse/held-out and
  muted rose precision delta; marker shape and direct labels supplement color.
- Surface: reproducible Matplotlib SVG/PDF/PNG, with editable SVG text and
  `panel_source.tsv`; accepted export is
  `evidence/abundance_classifier_validation_v1/figures_v2`.
- QA: exact input hashes and pairing validated before rendering; 123 SVG text
  nodes, zero raster nodes; PNG and rendered PDF visually inspected.

## Challenge diagnostic figure (development evidence)

- Question: which repeat architectures expose a weakness in the current method?
- Supported takeaway: report only the completed benchmark values; the first
  development baseline exposes indel-sensitive boundaries and one-array-per-read
  detection. It is not a final method-performance claim.
- Panels a-c: annotated heatmaps of array recall, array precision and strictly
  sequence-supported family recovery on the same 13 positive scenarios x 3 tools.
  Repeated chart family is deliberate: each panel compares the same experimental
  matrix with a different explicitly labelled denominator.
- Panel d: interval diagram from a deterministic first indel-case mismatch, with
  real truth and normalized tool intervals; use the same 0-based bp axis.
- Missing or failed data: grey cells labelled NA; never coerce to zero.
- Palette: blue sequential heatmap, grey missing cells; interval tracks use
  blue/gold/olive and neutral truth, with hatches and direct tool labels.
- Delivery: standalone PDF/SVG/PNG plus source TSV and source SHA-256 manifest.
- QA: inspect full PNG; verify matrix values against TSV, consistent [0,1]
  limits, complete axis labels, panel lettering and no raster image in SVG.
- Timing from this development run is exploratory because other development
  checks and a separate real-read pilot ran on the host. Rerun publication
  timing in isolation after freezing the source and settings.

## Final manuscript structure (pending evidence)

1. Read-to-assembly method, uncertainty model and output interfaces (a-f).
2. Simulated detection, family recovery, false calls and resource scaling (a-f).
3. Copy-number calibration, error, interval coverage and collapse detection (a-f).
4. Independent plant read/assembly applications and known-repeat agreement (a-f).
5. Probe prioritization and independently supported biological examples (a-d),
   only if sufficient evidence is obtained.
6. Supplement: seed sensitivity, comparator parameters, ablations, bias checks,
   cross-platform checks and source-data provenance; each figure multi-panel.
