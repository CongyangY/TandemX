# TandemX manuscript v5 figure and table plan

This plan implements the five-Result structure in `manuscript_v5.md`. It maps
every panel to existing evidence and distinguishes reusable source material from
a finished manuscript composite. No sixth Result or Figure 6 is included. A
future newly discovered tandem repeat and its independent sequence or FISH
validation will be evaluated before any sixth Result is written.

## Main-figure sequence

### Figure 1. TandemX links raw-read tandem-repeat families to assembly representation

- **Purpose:** make the distinct problem understandable before any benchmark.
- **Panels:** (A) a long read-supported satellite array versus a shorter array
  in a continuous assembly; (B) periodic intervals, family construction,
  diagnostic-k-mer abundance and localization of the same representative; (C)
  paired read/assembly estimates and supported or unresolved evidence states.
- **Status:** the editable Figure 1 drafts were rejected and retained only for
  provenance. The accepted biology-first AI raster is
  `paper/0910/figures/figure1_ai_biology_first/figure1_ai_biology_first_v1.png`.
  The author may add labels or terminology later; do not restore the rejected
  editable drafts.
- **Sources:** `docs/algorithms.md`, `docs/file_formats.md` and the public-command
  behavior documented in `README.md`.
- **Boundary:** orthogonal Illumina/ONT tests are validation analyses, not claimed
  as additional native TandemX input modes.

### Figure 2. Endpoint-matched comparisons define analytical scope

- **Purpose:** compare the final production workflow with external methods at
  matched endpoints, without using internal development history as the main
  narrative.
- **Panels:** (A) production TandemX, TRF and TideHunter family recovery and
  interval/base-union performance; (B) production TandemX, SRF k=151/k=101 and
  competitive-mapping family recovery and positive-family MARE; (C)
  shared-fragment background attribution.
- **Reusable assets:** `paper/evidence/cascade_gap_free_validation_v1/figures_v2/`,
  `paper/evidence/factorial_multik_replay/figures/`,
  `paper/evidence/cascade_native_screen_heldout_v1/figures_v2/` and
  `paper/evidence/srf_formal_unified_v1/condition_summary.tsv`.
- **Status:** earlier composites containing internal versions or incompatible
  binary endpoints were rejected and retained for provenance only. The current
  production-only candidate is
  `paper/0910/figures/figure2_v3_1_production_only_preview/`; it is regenerated
  from authoritative TSVs and remains pending author approval.
- **Supplementary move:** detailed runtime/RSS, negative controls, rejected
  development variants and tool-specific endpoint tables.
- **Coverage-saturation conclusion:** retain the concise 5x truth-recovery and
  20x operational-saturation result in this Results section. The complete
  four-panel curve is Supplementary Figure S1 so the main benchmark figure does
  not expand at the expense of the plant results. Source:
  `paper/evidence/discovery_saturation_validation_v1/figures/`.

### Figure 3. Operating range across plant long-read datasets

- **Purpose:** demonstrate applicability across contrasting plant genomes without
  presenting the cohort as replicated cross-species biology.
- **Panels:** (A) ten libraries/eight species and total bases; (B) read-length
  distributions; (C) matched-size call density and positive-read fraction in
  wheat, barley, rye and oat; (D) repeat-base fraction across nested scales.
- **Reusable assets:** `paper/evidence/multispecies_input_qc/figures_v2/` and
  `paper/evidence/multispecies_real_diagnostics/figures_v2/`.
- **Status:** the prior local composite remains a rejected editorial draft
  (`paper/0910/figures/figure3_operating_range/`). A source-linked replacement
  candidate is `paper/0910/figures/figure3_v2_operating_range/`; it separates
  call density and positive-read fraction onto aligned axes and is pending
  author approval.

### Figure 4. Improved assemblies preferentially recover read-prioritized families

- **Purpose:** make the historical-to-modern biological result the visual center
  of the manuscript.
- **Panels:** (A) two-system design; (B-C) Ey15-2 family gains and deficit/gain
  relationship; (D-E) Macadamia equivalents; (F) primary 8/8 and 2/2 summary.
- **Reusable assets:**
  `paper/evidence/ey15_donor_matched_collapse_v1/results/figures_v2/`,
  `paper/evidence/ey15_donor_matched_collapse_v1/results/evaluation_primary_total_bases/family_metrics.tsv`,
  `paper/evidence/macadamia_jansenii_donor_matched_collapse_v1/results/evaluation_primary_total_bases/family_metrics.tsv`
  and `paper/evidence/macadamia_jansenii_donor_matched_collapse_v1/headline_summary.json`.
- **Status:** the former local composite remains rejected provenance only. A
  source-linked two-species candidate is
  `paper/0910/figures/figure4_v2_assembly_recovery/`; it is regenerated from
  the frozen primary family tables and remains pending author approval.

### Figure 5. Orthogonal reads support three residual deficits and define three unresolved cases

- **Purpose:** end the present Results with independent biological evidence and a
  compact confidence boundary.
- **Panels:** (A) six-candidate design; (B) Ey15-2 supported families; (C)
  Macadamia TXF000496; (D) TXF001517/TXF000563 discordance; (E) TXF000695
  cross-k depth distributions; (F) supported versus unresolved evidence summary.
- **Reusable assets:** `paper/evidence/orthogonal_abundance_validation_v1/figures_v2/`,
  `paper/0910/source_data/figure6_txf000695_qc.tsv` and
  `paper/0910/source_data/figure6_txf000695_qc_notes.md`.
- **Status:** orthogonal source figure and TXF000695 source data exist. A
  source-linked Figure 5 candidate at
  `paper/0910/figures/figure5_v2_orthogonal_abundance/` places the complete
  TXF000695 cross-k diagnostic in panel E and the six-family evidence state in
  panel F; it is pending author approval.

## Main tables

### Table 1. Native analytical scope of relevant tools

- Replaces the former main-text public-library table.
- Presents inputs and native outputs rather than a simplistic winner/loser matrix.
- Capability descriptions are based on the cited primary papers for TideHunter,
  SRF, TRASH, TandemTools, Merqury, KAT and AniAnn's.
- TandemX orthogonal validation is not listed as a native software input.

### Table 2. Family-level biological evidence summary

- Summarizes the primary Ey15-2 and Macadamia denominators, historical/newer
  prioritization and orthogonal outcomes.
- Values are 19 eligible/8 recovered/8 prioritized/3 tested/2 supported/1
  unresolved for Ey15-2 and 43/2/2/3/1/2 for Macadamia.

## Supplementary tables and figures

- **Table S1:** public HiFi cohort, moved from former main Table 1; actual table
  in `supplementary_tables_v5.md`, full source in `paper/tables/input_cohort.tsv`.
- **Table S2:** unified SRF abundance matrix, moved from former main Table 2;
  actual table in `supplementary_tables_v5.md`.
- **Figure S1:** pre-specified discovery-saturation curve across three simulated
  validation genomes.
- **Figure S2:** specified simulated production endpoint and negative controls,
  from `supplementary_figures/figure_s2_production_validation/`.
- **Figure S3:** production-only simulation abundance and assigned-base endpoint,
  from `supplementary_figures/figure_s3_production_quantification/`.
- **Figure S4:** evaluable TideCluster endpoint cells and the external-tool
  completion record, from `supplementary_figures/figure_s4_tidecluster_endpoint/`.
- **Figure S5:** held-out simulated endpoint and negative controls for TandemX,
  TideHunter and TRF, from
  `supplementary_figures/figure_s5_heldout_external_endpoint/`.
- **Figure S6:** formal SRF endpoint record, including the shared-catalogue
  boundary for competitive mapping, from
  `supplementary_figures/figure_s6_formal_srf_endpoint/`.
- **Figure S7:** full plant-input QC and descriptive TandemX calls, from
  `supplementary_figures/figure_s7_multispecies_inputs/`.
- **Figure S8:** Ey15-2 donor-matched ratio, read-depth sensitivity and alignment
  context, from `supplementary_figures/figure_s8_ey15_reference_context/`.
- **Figure S9:** Macadamia donor-matched ratio, read-depth sensitivity and
  alignment context, from
  `supplementary_figures/figure_s9_macadamia_reference_context/`.
- **Figure S10:** complete six-candidate orthogonal abundance evidence and final
  evidence state, from `supplementary_figures/figure_s10_orthogonal_candidates/`.
- **Figure S11:** TXF000695 diagnostic-k-mer QC and unresolved abundance context,
  from `supplementary_figures/figure_s11_txf000695_qc/`.
- **Tables S3-S4:** depth-wise discovery/recall and adjacent-depth marginal-yield/
  Jaccard summaries; complete per-seed values remain in the compact evidence
  archive.
- **Supplementary Result:** the bounded Ey15-2 recovery experiment remains a
  complete negative result and is not returned to the main narrative.

The source-artwork review and required redraws are recorded in
`paper/0910/supplementary_figure_visual_qa_v5.md`. S1 remains the original
directly reusable figure. S2--S11 now have source-linked v5 artwork under
`paper/0910/supplementary_figures/`; each directory includes `panel_source.tsv`,
PNG, single-page PDF and editable SVG.
