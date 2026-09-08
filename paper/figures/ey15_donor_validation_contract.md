# Ey15 donor-matched validation figure contract

## Analytical question and takeaway

Does a read-based TandemX estimate identify repeat families that are
under-represented in the older Ey15-2 CLR-Canu assembly when the donor-matched
HiFi-Hifiasm assembly is used as a retrospective high-quality reference proxy?

The figure must show both the primary classification result and its limits: the
15-kb denominator is small, missing-bp agreement is not perfect, the 5-kb
sensitivity denominator contains errors, normalization can change predictions,
and orthogonal annotation/alignment context is not independent copy truth.

## Data and grain

- One point per source-eligible family for the primary 15-kb panels; expected
  count is the frozen evaluator's full eligible denominator (currently 19).
- All families at each predeclared 5-, 15- and 50-kb new-assembly threshold for
  the threshold panel; no result-dependent family filtering.
- The same primary family universe for total-bases versus explicit-107x
  normalization.
- All eligible families in the post hoc author-annotation and frozen
  assembly-alignment context panels.

## Six-panel plan

1. New- versus old-assembly localized bp per family, with equality and frozen
   0.6 collapse lines. Preserve old=0 points with a labelled symmetric-log axis.
2. Read-estimated versus old-assembly bp, with equality and frozen 0.6
   prediction lines; use the same family labels and state encoding as panel A.
3. Predicted missing bp versus observed old-to-new gain with an equality line,
   all 19 points and Pearson/Spearman values in the annotation.
4. TP/FN/FP/TN counts at the predeclared 5-, 15- and 50-kb thresholds. Show the
   denominator above each threshold and do not hide adverse 5-kb outcomes.
5. Primary total-bases versus explicit-107x outcome/metric comparison after the
   sensitivity run completes. State that these are two normalization choices on
   the same reads, not biological replicates.
6. Orthogonal context by reference state: author-annotation overlap fraction and
   same-chromosome primary aligned-query coverage. Identify the annotation test
   as post hoc and the assembly alignment as a frozen explanatory audit.

## Rendering and palette

Use reproducible static Matplotlib and export editable-text SVG plus PDF and PNG.
Palette policy is a hard two-root cap: blue for reference-collapse/focal values,
orange for discordant or adverse outcomes, with charcoal/grey for retained and
reference elements. Point shape, open/filled markers and line style must repeat
the distinctions so color is not the only encoding. Use a white background,
quiet grids, visible axis anchors and no gradients.

## QA and outputs

Write a panel-source TSV that retains family identifiers, denominators, all
metrics and warnings. Write a provenance JSON with SHA-256 for every source and
output. Verify the SVG contains no raster image elements and retains editable
text. Render the PDF independently to PNG and inspect both it and the direct PNG
for clipping, overlap, detached labels and honest scales. Do not promote the
figure to the manuscript if either independent evaluator receipt fails.

