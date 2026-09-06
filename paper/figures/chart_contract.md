# Figure contracts

Publication graphics use Matplotlib, editable SVG text and embedded PDF fonts.
No AI-generated quantitative panels. Raw data and executable plotting code are
the authority; these specifications do not imply that a planned panel exists.

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
