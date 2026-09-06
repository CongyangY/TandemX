# Figure 7 chart contract

- **Question:** Does the post-failure depth-gated rule reproduce its development
  benefit on untouched genome seeds while preserving baseline FPR and precision?
- **Takeaway:** On seeds 5801-5803, sensitivity rises by 0.058985 with unchanged
  FPR and a 0.001608 precision gain; every seed passes the frozen guardrails.
- **Grain:** 2,430 paired family-conditions from three genome seeds; 810 rows have
  estimated haploid depth below 2 and 1,620 have depth at least 2.
- **Panels:** sequence/status cards; grouped delta bars; grouped confusion-count
  bars; seed-level grouped delta bars; coverage-level grouped delta bars; grouped
  localization-recall bars.
- **Palette:** gray baseline, blue v3, green independent pass, orange adverse
  evidence/FPR, magenta precision. Labels and bar position duplicate color.
- **Scale:** count panels start at zero. Delta panels include zero. Panel F uses
  a labelled 0.90-1.012 focused scale because all six values exceed 0.95.
- **Outputs:** editable SVG, vector PDF, 240-dpi PNG, source TSV and provenance
  JSON in `figures_v2/`.
- **QA:** inspect PNG and rendered PDF; require no clipping or overlap, zero SVG
  raster image nodes, editable text and hashes for every input/output artifact.
