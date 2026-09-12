# Manuscript v5 supplementary-figure visual QA

Reviewed on 2026-09-12 against `manuscript_v5.md` and
`manuscript_v5_figure_table_plan.md`. This is a review of currently available
evidence artwork. Source figures remain provenance; they are not silently
relabeled as finished v5 supplements. The completed v5 redraws are listed
explicitly below and live in `paper/0910/supplementary_figures/`.

| Planned figure | Existing source artwork | Visual/content review | Resolution for v5 |
| --- | --- | --- | --- |
| S1, discovery saturation | `paper/evidence/discovery_saturation_validation_v1/figures/discovery_saturation.pdf` | Four panels are legible, have coherent depth axes and explicitly separate planted-family recall from catalogue stability. Threshold and 20x marker are visible. | **Pass.** Retain as S1 with its current source table and v5 legend. |
| S2, detailed discovery | `paper/evidence/cascade_gap_free_validation_v1/figures_v2/cascade_validation.pdf` | Legible but is an older cascade-development audit. Its “all 14 preregistered gates pass” title and runtime/RSS panels do not describe the final external-endpoint narrative. | **Rebuilt.** `figure_s2_production_validation/` retains only specified production endpoints and negative controls. |
| S3, quantification | `paper/evidence/factorial_multik_replay/figures/factorial_quantification.pdf` | Layout is readable, but it contains a retired prototype, development data and calibration limitations. | **Rebuilt.** `figure_s3_production_quantification/` uses frozen production estimator rows only. |
| S4, localization/comparator detail | `paper/evidence/tidecluster_factorial_validation_v3/results/figures_v2/tidecluster_factorial_validation.pdf` | Six-panel rendering is legible and retains two failed cells as NA. It mixes native TideCluster stages, resource traces and conditional assembly validation. | **Rebuilt.** `figure_s4_tidecluster_endpoint/` shows evaluable endpoint cells plus the explicit completion record. |
| S5, held-out/rejected development variants | `paper/evidence/cascade_native_screen_heldout_v1/figures_v2/cascade_heldout.pdf` | The rendering exposes historical promotion gates rather than the final external-comparator endpoint. | **Rebuilt.** `figure_s5_heldout_external_endpoint/` uses held-out endpoint and negative-control values only. |
| S6, formal SRF comparison | `paper/evidence/srf_formal_unified_v1/condition_summary.tsv` and `summary.tsv` | The former source was tabular only. | **Built.** `figure_s6_formal_srf_endpoint/` keeps the shared-catalogue boundary for competitive mapping explicit. |
| S7, full plant-input QC and diagnostics | `paper/evidence/multispecies_input_qc/figures_v2/multispecies_input_qc.pdf`; `paper/evidence/multispecies_real_diagnostics/figures_v2/real_comparator_diagnostics.pdf` | Both PDFs are legible and preserve the non-ranking resource boundary. | **Built.** `figure_s7_multispecies_inputs/` contains input QC and descriptive TandemX output only. |
| S8, Ey15-2 alignment context and eligibility sensitivity | `paper/evidence/ey15_donor_matched_collapse_v1/results/figures_v2/ey15_donor_validation.pdf` | The older panel wording conflicts with v5’s abundance boundary. | **Rebuilt.** `figure_s8_ey15_reference_context/` contains ratio, depth sensitivity and alignment context without physical-sequence inference. |
| S9, Macadamia alignment context and eligibility sensitivity | primary/reported-depth and old-new alignment tables | The former source was tabular only. | **Built.** `figure_s9_macadamia_reference_context/` retains only source-linked ratio and alignment context. |
| S10, complete orthogonal candidate/context-family analysis | `paper/evidence/orthogonal_abundance_validation_v1/figures_v2/orthogonal_abundance_validation.pdf` | The source retains all six candidates and the family-level evidence boundary. | **Rebuilt.** `figure_s10_orthogonal_candidates/` retains six candidates and three supported/three unresolved states. |
| S11, TXF000695 diagnostic QC and unresolved context | `paper/0910/source_data/figure6_txf000695_qc.tsv`; `paper/0910/source_data/figure6_txf000695_qc_notes.md` | The source identifies the cross-k boundary and leaves the family unresolved. | **Built.** `figure_s11_txf000695_qc/` presents the two frozen settings separately and retains the unresolved interpretation. |

## Completed repairs in this figure pass

- Main Figure 4 was replaced with the source-linked candidate in
  `figures/figure4_v2_assembly_recovery/`. Its continuous labels now read
  `read--assembly abundance deficit` and `assembly gain`, rather than missing
  physical bases. The old Figure 4 remains rejected provenance.
- Main Figure 5 now absorbs the former Figure 6 TXF000695 material as panels E
  and F. The full diagnostic expansion remains an explicit S11 deliverable.
- S2--S11 were rebuilt as source-linked PNG, single-page PDF and editable SVG
  assets, each with a figure-local `panel_source.tsv` and README. Visual review
  found no clipped values, overlapping panel text or raster embedding in the
  rebuilt SVG files.
- No method, comparator, estimator, Result section or biological conclusion was
  added during this artwork review.
