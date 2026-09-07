# Figure version audit

- `figures_v1`: rejected after visual QA because the panel-D legend overlapped
  the stacked outcome bar.
- `figures_v2`: visual layout accepted, but rejected for source-data reuse
  because panel-F runtime and RSS rows used coverage alone rather than a unique
  execution key.
- `figures_v3`: accepted. The legend is clear, the six panels contain no raster
  SVG elements, and all 54 runtime/RSS pairs have unique
  `seed:condition:coverage` keys in `panel_source.tsv`.

All versions retain the same scientific values. Use `figures_v3` for the paper.
