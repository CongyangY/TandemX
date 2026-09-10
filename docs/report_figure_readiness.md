# Report-figure readiness audit

Audit date: 2026-09-10. This records the boundary between the reusable
software report renderer and the manuscript figure package. It does not claim
that the manuscript figures are complete.

## Software report exports

The current `write_report_figures(outdir, data)` renderer provides six figure
classes:

| Class | SVG/PDF/PNG | Source TSV | Receipt | Assessment |
| --- | --- | --- | --- | --- |
| `summary` | yes | yes | yes | usable as a compact run overview |
| `family_abundance_vs_assembly` | yes | yes | yes | usable for family-level relationship review |
| `top_underrepresented` | yes | yes | yes | usable for ranked positive deficits; reports paired read/assembly values |
| `family_landscape` | yes | yes | yes | usable for length versus abundance/support context |
| `family_hierarchy` | yes | yes | yes | usable as a selected heuristic relationship view |
| `family_evidence_cards` | yes | yes | yes | usable for up to six family evidence summaries |

Each receipt records output hashes, byte sizes, and SVG text/raster node counts.
The accepted real-data QA used the frozen Ey15-2 discovery catalogue with the
primary quantification and old-assembly localization inputs. It generated
2,133 families and 1,040 architecture edges; the hierarchy view disclosed
6/1,040 selected edges. Direct PNG and independently rasterized PDF checks
passed after the evidence-card layout fix. Toy, empty, and no-assembly paths
also rendered without layout warnings under `python -W error`.

The accepted real-data card review is retained at
`/tmp/tandemx-ey15-report-qa-20260910-cards-final/figures/family_evidence_cards.png`
and its independent PDF raster at
`/tmp/tandemx-ey15-report-qa-20260910-cards-final/pdf_render/family_evidence_cards-1.png`.
Warnings are displayed as human-readable summaries with explicit further-warning
counts; complete warning strings remain in the source TSV.

The report figures are suitable for offline HTML embedding and software-level
evidence review. They preserve unresolved values and warnings, and describe
architecture links as heuristic candidates rather than validated HORs.

Final parent visual review additionally rejected the first T7 report layout:
the deficit legend overlapped its last row, same-column hierarchy edges
occluded each other, and landscape symlog padding showed meaningless negative
abundance ticks. The corrected report uses an external deficit legend,
six separate pairwise relationship rows, explicit confidence groups in the
scatter, and read/assembly landscape panels with nonnegative ticks, shared GC
and source-consistent support-size encoding. The final report directory is
`/Volumes/T7/Codex/TandemX/results/ey15_family_report_recovery_20260910_v2/`.
The original T7 report remains retained as layout-rejected evidence. No
scientific stage was rerun to make these display corrections.

## Manuscript gap assessment

The current manuscript construction plan requires six biological figures with
different evidence packages:

1. An editable-vector biological problem and workflow schematic. The report
   figures do not provide this schematic.
2. An integrated controlled-benchmark figure covering recovery, abundance
   error, localization, binary metrics, comparator results, runtime, and
   memory. The report figures do not assemble these benchmark panels.
3. A cross-species, nested-sampling, provenance-aware plant landscape figure.
   The report landscape is family-level and single-run; it is insufficient.
4. A two-species historical/newer assembly reference-proxy figure with
   denominator sensitivities and explanatory alignment context. The report
   scatter is only one family-level view and is insufficient.
5. An orthogonal-read validation figure with the six preselected families,
   Illumina k values, ONT direction, and final evidence matrix. The report
   figures do not contain these validation panels.
6. A focused TXF000695 cross-k/cross-platform confidence-boundary figure. The
   report figures do not contain diagnostic-word distributions or this QC
   analysis.

The manuscript also calls for exact panel source tables, supplementary figure
migration, editable core graphics, legend cross-references, and final
reference/metadata checks. Those remain manuscript-package tasks.

## Decision boundary

Software report-figure completion is not equivalent to manuscript-figure
completion. The six report classes are reusable diagnostic views with source
receipts; they do not replace the six planned evidence-specific manuscript
figures, their source-data tables, or the required biological and benchmark
panel assembly. A manuscript figure should be marked ready only after its
panel-level evidence, source table, legend, editable-artwork requirements, and
direct/PDF visual QA have each passed independently.
