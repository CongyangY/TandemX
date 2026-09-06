# Domain-shift assembly-comparison validation

This archive applies the previously frozen exact-copy multi-k/depth rule to
fresh seeds 5401–5403 after crossing 1%, 3% and 5% founder-to-unit substitutions
with one or three same-family array segments. Three-segment arrays contain
500-bp independent interruptions. The configuration was committed in `c15dad7`
and passed hosted Ubuntu/macOS CI before the seeds were consumed once. The
baseline completed 1,062/1,062 commands. The multi-k evaluation reused those
inputs and outputs without rerunning baseline commands or fitting held-out data.

Across 2,430 paired family conditions, the single-k baseline had
TP/FN/FP/TN=1230/228/625/347: sensitivity 0.843621, false-positive rate
0.643004 and precision 0.663073. The frozen hybrid had 1263/195/663/309:
sensitivity 0.866255, false-positive rate 0.682099 and precision 0.655763.
The sensitivity gain therefore accompanied 38 additional false positives and
lower precision. Multi-k alone had 270 unavailable rows and is not a complete
method comparison.

The failure is strongly associated with localization. At 1% unit divergence,
full-assembly mean base recall was 0.754905 for one segment and 0.692814 for
three segments. At 3%, both means were 0.000793; at 5%, both were zero. Under
3–5% divergence, false-positive rates were 0.858025–0.888889 for the baseline
scenario strata and 0.888889 for the frozen hybrid strata. This exposes the
current exact diagnostic-k-mer support filter as a domain-shift boundary. It
does not justify tuning on consumed seeds or claiming robust collapse detection.

`archive_manifest.json` covers 19 files and rechecks the frozen calibration,
baseline receipts, source snapshots, complete scenario matrix, exact paired
keys and confusion counts. `figures_v3/` is the accepted six-panel render; its
SVG has 104 editable text nodes and no raster image node, and its PDF was
rendered to PNG and visually inspected without overlap or clipping. `figures_v1/`
and `figures_v2/` are retained as rejected layouts because bottom labels
overlapped the footer. All versions retain source rows and output hashes.

This is a small known-catalogue simulation with independent length-preserving
substitutions and same-chromosome interruptions. It is not an empirical model of
satellite evolution, graph/contig fragmentation, unknown-catalogue discovery or
biological assembly collapse.
