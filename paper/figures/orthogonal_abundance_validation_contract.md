# Orthogonal abundance validation figure contract

## Analytical question and takeaway

Do independent Illumina and ONT observations support TandemX abundance above a
newer assembly for the six frozen deficit candidates, or do they instead expose
HiFi quantification bias?

The figure must show the answer and its boundary together: three families have
stable orthogonal support for residual under-representation, three are
unresolved, no family has stable cross-k support for quantification bias, and
continuous abundance/deficit magnitude is not physical missing-base truth.

## Data and grain

- Six candidate rows, fixed before orthogonal result inspection: three Ey15-2
  and three Macadamia families.
- Assembly, frozen HiFi, Illumina k=21 and Illumina k=31 abundance per family.
- Direction-only ONT abundance ranges for the three Macadamia candidates across
  the frozen 616--780-Mb genome-size sensitivity.
- One binary evidence state for each k/platform/family and one final
  interpretation per family.
- No aggregation of the selected 3/6 result into an accuracy rate.

## Three-panel plan

1. Absolute abundance on a log scale: newer-assembly localized bp, frozen HiFi,
   Illumina k=21/k=31 and Macadamia ONT range. Directly label family/species.
2. Newer-assembly/read ratio on a log scale with the frozen 0.6 threshold. Show
   Illumina k=21/k=31 and Macadamia ONT sensitivity ranges. Ratios below 0.6 are
   directionally consistent with under-representation, not proof of physical
   missing bp.
3. Evidence-state matrix for HiFi, Illumina k=21, Illumina k=31, ONT and final
   interpretation. Use symbol shape/text as well as fill so color is not the
   only distinction. Preserve not-available and unresolved states.

## Rendering and palette

Use reproducible static Matplotlib and export editable-text SVG plus PDF and
PNG. The palette has a hard two-root cap: blue for residual-support evidence,
orange for bias-direction evidence, charcoal for assembly/HiFi/reference text,
and grey/open marks for unresolved or unavailable states. Use quiet grids, no
gradients and consistent log scales.

## QA and outputs

Write `panel_source.tsv`, `figure_legend.md` and `receipt.json` with all source
and output hashes. Verify six candidate rows, the 3 supported/3 unresolved
final split, editable SVG text, zero SVG raster images and no non-finite plotted
values. Render the exported PDF independently to PNG and inspect it together
with the direct PNG for clipping, overlaps, detached labels and honest scales.
