# TandemX manuscript v5 figure and table plan

This plan implements the five-Result structure in `manuscript_v5.md`. It maps
every panel to existing evidence and distinguishes reusable source material from
a finished manuscript composite. No sixth Result or Figure 6 is included. A
future newly discovered tandem repeat and its independent sequence or FISH
validation will be evaluated before any sixth Result is written.

## Main-figure sequence

### Figure 1. TandemX links raw-read tandem-repeat families to assembly representation

- **Purpose:** make the distinct problem understandable before any benchmark.
- **Panels:** (A) continuity versus family representation; (B) raw reads to
  operational family catalogue; (C) read abundance *R*, assembly representation
  *A* and deficit; (D) family-level output and evidence state.
- **Status:** new editable schematic required.
- **Sources:** `docs/algorithms.md`, `docs/file_formats.md` and the public-command
  behavior documented in `README.md`.
- **Boundary:** orthogonal Illumina/ONT tests are validation analyses, not claimed
  as additional native TandemX input modes.

### Figure 2. Controlled validation defines analytical scope

- **Purpose:** present one compact validation figure rather than a sequence of
  development-history sections.
- **Panels:** (A) truth/endpoints; (B) discovery and read-interval recovery;
  (C) abundance error; (D) assembly localization; (E) unified SRF and competitive-
  mapping comparison.
- **Reusable assets:** `paper/evidence/cascade_gap_free_validation_v1/figures_v2/`,
  `paper/evidence/factorial_multik_replay/figures/`,
  `paper/evidence/cascade_native_screen_heldout_v1/figures_v2/` and
  `paper/evidence/srf_formal_unified_v1/condition_summary.tsv`.
- **Status:** evidence complete; final five-panel composite still required.
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
- **Status:** source panels exist; final composite and consistent styling required.

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
- **Status:** Ey15 source figure and both family tables exist; a balanced two-
  species composite remains to be built.

### Figure 5. Orthogonal reads support three residual deficits and define three unresolved cases

- **Purpose:** end the present Results with independent biological evidence and a
  compact confidence boundary.
- **Panels:** (A) six-candidate design; (B) Ey15-2 supported families; (C)
  Macadamia TXF000496; (D) TXF001517/TXF000563 discordance; (E) TXF000695
  cross-k depth distributions; (F) supported versus unresolved evidence summary.
- **Reusable assets:** `paper/evidence/orthogonal_abundance_validation_v1/figures_v2/`,
  `paper/0910/source_data/figure6_txf000695_qc.tsv` and
  `paper/0910/source_data/figure6_txf000695_qc_notes.md`.
- **Status:** orthogonal source figure and TXF000695 source data exist; the old
  Figure 6 material must be condensed into Figure 5E-F.

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
- **Figures S2-S6:** detailed discovery, quantification, localization, comparator
  and performance evidence contributing to main Figure 2.
- **Figure S7:** full plant-input QC and unmatched real-data diagnostics supporting
  Figure 3.
- **Figures S8-S9:** assembly-alignment context and eligibility sensitivities for
  Ey15-2 and Macadamia.
- **Figure S10:** complete candidate/context-family orthogonal analysis.
- **Figure S11:** full TXF000695 diagnostic distributions and mapping sensitivities.
- **Tables S3-S4:** depth-wise discovery/recall and adjacent-depth marginal-yield/
  Jaccard summaries; complete per-seed values remain in the compact evidence
  archive.
- **Supplementary Result:** the bounded Ey15-2 recovery experiment remains a
  complete negative result and is not returned to the main narrative.
