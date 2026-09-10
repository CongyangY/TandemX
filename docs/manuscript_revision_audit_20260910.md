# Manuscript revision audit — 2026-09-10

This is an editorial audit of `paper/0910/manuscript_v2.md`; it does not alter
the manuscript, algorithms, thresholds or result tables. It compares the
manuscript with the three supplied editorial analyses:

* `/Users/ycy/.codex/attachments/790959cf-d401-4958-bbd4-05f3e61122ed/pasted-text.txt`
* `/Users/ycy/.codex/attachments/954dc7ad-973e-438c-a87b-2f4329447b27/pasted-text.txt`
* `/Users/ycy/.codex/attachments/9e1d1a77-4719-4b26-ac6e-81d19a05318a/pasted-text.txt`

The three analyses converge on a focused article: present TandemX as a
read-first, assembly-aware measurement framework, retain the negative and
unresolved evidence, and move implementation history and recovery mechanics
out of the main narrative. This is an editorial recommendation, not a new
scientific result.

## Recommended edits by manuscript line

Line numbers refer to the current `manuscript_v2.md` and should be refreshed
after edits.

| Location | Finding | Suggested action / replacement |
| --- | --- | --- |
| 145–164, 580–593, 649–654 | The definition of *R*, *A*, deficit, the 0.6 rule and unresolved fates is stated three times. | Keep the compact biological definition at 145–164 and the exact operational rule at 649–654. In the Methods definition block (580–593), retain only terminology needed to interpret tables, for example: “We report read-derived abundance, assembly representation and their positive difference; a binary candidate flag used the frozen *A/R* < 0.6 rule, while discordant evidence remained `unresolved`.” Remove repeated explanatory sentences elsewhere. |
| 149–154, 283–285, 308–309, 355–364, 551–569 | “Deficit is not physical missing sequence / assembly is not truth” is repeated in Results, Discussion and Limitations. | Preserve the full evidence boundary once at 149–154, a short reminder at the first old/new conclusion (283–285), and the consolidated limitations paragraph (551–569). At 308–309 and 355–364 use a short result sentence: “The magnitude remains a read–assembly abundance estimate, not a physical missing-base measurement.” |
| 156–164, 368–375, 525–538, 735–742 | The same unresolved logic is explained generically and then repeated for each discordant family. | Keep the family-specific reasons at 368–375 and the interpretation rule in Methods 735–742. Shorten the generic Discussion passage to: “Discordant cross-k or platform evidence was retained as `unresolved`; heuristic QC labels were not promoted to a classifier.” Keep TXF000695’s mechanism uncertainty (525–538), because it is a principal result. |
| 174–181, 799–833 | Report directory, configuration interface and “interface additions” read as project-status text. | Move command/interface details to Software Availability or Supplementary Methods. Replace 174–181 with: “The workflow emits an offline report linking family abundance, assembly representation, confidence labels and warnings to source tables.” Keep the reproducibility claim only if the corresponding artifacts are included in the final package. |
| 237–247, 493–498 | Comparator and failure-history language is useful, but the development-history framing is more prominent than the biological result. | Keep unavailable cells as unavailable (237–238) and the task-matched comparator conclusion (240–247, 493–498). Move parameter-selection history, resource profiles and failed-run chronology to Supplementary Results; do not replace unavailable values with zero. |
| 445–471, 759–797 | The bounded recovery test occupies a complete Results subsection plus detailed Methods, although it did not recover additional repeat representation and does not establish a method. | Move the detailed recovery subsection to Supplementary Results/Methods. In the main text retain a short negative-result paragraph: “A bounded Ey15-2 flank-and-read test produced one unpolished span, but it did not increase measured repeat representation or improve agreement with the newer reference proxy. The test therefore does not establish a recovery method; the complete loci, filters and adverse outcomes are provided in Supplementary Results.” Preserve the non-generalizability statement (469–471) in the supplement. |
| 467–471 | The “does not establish…” caveat is repeated in both the result and the surrounding project-scope language. | Keep the first sentence as the scientific conclusion. Move “do not justify further algorithm expansion for this release” to an internal release note or Supplementary Methods; it is project management rather than a biological result. |
| 262–274, 681–690 | Denominator sensitivity (5, 15 and 50 kb) and old/new proxy rules are clearly described, but they can be mistaken for threshold sensitivity. | Label these explicitly as “eligibility-denominator sensitivity.” Do not call them 0.6-threshold sensitivity. If threshold sensitivity is added later, report it in a separate table with threshold, eligible families, TP/FP/FN/NA and the frozen analysis snapshot. |
| 835–842 | “Versioned release, Bioconda recipe and Zenodo deposition remain in preparation” is internal delivery status. | Move this sentence to the release checklist. In the manuscript, state only currently available source/evidence locations and do not imply a DOI or release identifier. |
| 854–949 | The section is titled “Figure legends and construction plan” and repeatedly says figures “require,” “will be assembled,” “will be rebuilt” or “must be generated.” | Keep final legends in the manuscript only after the corresponding artwork and source data are accepted. Move construction instructions, migration decisions and accepted/rejected layout history to an internal figure manifest or Supplementary package notes. The scientific panel descriptions can remain as legends after rewriting in past/present tense. |
| 950–974 | The migration table is valuable project bookkeeping, but it interrupts the manuscript narrative. | Move the whole migration map to the package/repository documentation. Retain one manuscript sentence pointing to Supplementary Figures and Source Data, without commit IDs, evidence IDs or transfer history. |
| 1029–1031 | Reference-manager and author-information status is an unfinished internal checklist. | Remove from the manuscript. Resolve reference 12 and dataset metadata before submission and track the work in the editorial checklist. Do not infer missing author, funding or competing-interest information. |

## Caveats that must remain visible

The following are evidence limits rather than stylistic repetition and should
survive condensation:

1. The newer assembly is a reference proxy, not absolute copy-number truth
   (283–285, 317–319, 688–690).
2. The continuous read–assembly deficit is not a physical missing-base count
   (149–154). Binary state and continuous magnitude are separate quantities.
3. The 0.6 rule is frozen and analytical, not a universal biological boundary
   (649–654). Discordant k values, platforms and applicable genome-size
   sensitivities remain `unresolved`.
4. The orthogonal set contains six preselected candidates; “three supported and
   three unresolved” is not accuracy, prevalence or a genome-wide collapse rate
   (355–364, 551–558).
5. Macadamia metadata do not prove identical DNA extraction for all libraries
   (551–556, 692–699); donor/material matching must remain explicit.
6. ONT supplies direction-level occupancy evidence, not exact copy-number truth
   (346–353, 520–523).
7. Failed, unstarted or unavailable comparator cells remain unavailable rather
   than zero (237–238, 664–667).
8. TXF000695’s mechanism is unresolved; its QC labels are descriptive and not a
   validated classifier (525–538, 744–757).
9. The recovery test is bounded and negative; moving it to the supplement must
   preserve the candidate span, no-increase result and explicit scope (445–471).

## What can be consolidated in Methods or Limitations

The following generic wording need not recur in every Results paragraph:

* **Methods:** the definitions of abundance, assembly representation, deficit,
  *A/R* < 0.6, eligibility fates, unavailable evidence, and the rule that
  cross-k/platform disagreement is `unresolved` (580–593, 649–654, 735–742).
* **Limitations:** sampling and genome-size normalization, family specificity,
  donor/material matching, platform/library effects, lack of physical repeat
  truth, simulated-data simplifications, and the boundary against claiming
  megabase-array reconstruction (551–569).
* **Supplementary Methods/Results:** complete command settings, seeds, hashes,
  resource profiles, adverse rows, operator-stopped runs, recovery filters and
  migration history (245–247, 759–797, 950–974).

The main Results should then use one short evidence qualifier at the point of
each claim. This keeps uncertainty attached to the relevant result without
repeating the same disclaimer in every section.

## 0.6 threshold sensitivity audit

The current evidence archive shows that 0.6 is exercised as the frozen
baseline rule and that some held-out tables contain a 0.5 comparison. For
example, these source tables contain a `decision_threshold` column:

* `paper/evidence/abundance_multik_collapse_heldout/comparison_metrics.tsv`
* `paper/evidence/abundance_localizer_multik_heldout/comparison_metrics.tsv`
* `paper/evidence/abundance_classifier_depth_gated_validation_v1/heldout/comparison_metrics.tsv`

The observed values in these tables are 0.5 and 0.6, with different row counts
and outcomes. The development code also records a decision-threshold grid in
`benchmarks/abundance/run.py` and frozen-selection checks in
`benchmarks/abundance/select_robust_blend.py`. This is evidence of threshold
evaluation in the controlled classifier-development branch, not proof that a
manuscript-level sensitivity analysis of the final old/new and orthogonal
catalogues has been completed.

I did not find a final, catalogue-level 0.6 sensitivity table reporting the
family calls and evidence fates across a declared threshold grid. The current
manuscript’s 5/15/50-kb results (262–274, 681–690) are denominator sensitivity,
and the genome-size and cross-k results are different perturbations. Therefore
the audit status is **partial threshold evidence; final 0.6 threshold
sensitivity is not verified as complete**. Do not write that the 0.6 result is
threshold-robust. If the parent task later elects to complete this bounded
analysis, reuse the existing comparison-metrics schema and add at least:
`threshold`, `family_id`, `R`, `A`, `A/R`, continuous deficit, binary state,
confidence/evidence fate, and the analysis snapshot. No rerun was performed for
this audit.

## Software report figures versus manuscript figures

The six offline software report figure classes (summary, abundance-versus-
assembly, top under-represented families, landscape, hierarchy and evidence
cards) are useful provenance-backed report outputs. Their SVG/PDF/PNG exports,
source TSVs and receipts are tracked in `docs/report_figure_readiness.md`.
That software-report completion does not complete the manuscript figure plan:
the manuscript still requires its stated workflow schematic, controlled
benchmark composite, cross-species/assembly comparisons, orthogonal validation,
and TXF000695 boundary figure with final legends and source-data mappings.
Report plots may support those panels, but they should not be described as the
manuscript figures until panel selection, statistical annotations, artwork and
source-data links are finalized.
