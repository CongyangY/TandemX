# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Thirteen inspected six-panel figures, including exact-copy validation, adverse domain shift, failed/passed classifier validations, ten-library QC, failed cascade promotion, quantify development, frozen quantify validation and passed guarded-cascade validation; two four-panel single-library QC figures | Isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Ten-library/eight-species input table; complete per-condition and native comparison tables; exact-output 11/111-Mb discovery optimization table | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Six independent factorial quantify genomes/54 conditions; three other predeclared three-genome held-out matrices; the consumed failed cascade split and fresh seed-2201 guarded-cascade validation; ten full library QCs across eight species; several real comparator diagnostics | Family/process-held-out evaluation, broader species/material replication, true biological recall and multiple isolated large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained. TideCluster 1.21.2 now has a pinned runnable image and successful planted-truth assembly smoke with all dependency/resource evidence archived | Plant-scale TideCluster and further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | The 108-run development ablation and exact-output speed replay completed. On frozen seeds 6401–6403, 54/54 public commands and 9/9 gates passed: MARE 0.408767→0.363222 with positive improvement in every genome and 661/495/329 improved/equal/worse family pairs. Ungated controls were slightly lower at 0.359443. The rule is an opt-in public CLI mode | Real control specificity, discovered-catalogue propagation, calibrated sampling intervals, catalogue uncertainty and ploidy remain |
| Assembly interpretation | The failed 5701–5703 result was retained. A subsequently frozen depth gate passed fresh 5801–5803, increasing sensitivity 0.594650→0.653635 with unchanged FPR 0.016461 and precision 0.981880→0.983488; localizer recall/precision were 0.980059/0.999617 | Obtain biological under-representation and donor-matched validation; do not generalize the known-catalogue IID result |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim; the adverse split and a later successful frozen-rule split are both retained | Establish reproducible added value against task-matched external tools and biological truth; AI only if independently useful |
| Software quality | Frozen cascade-validation commit 054b935 passed 531 Python tests and hosted checks before seed 2201 was used. Published archive commit 9b09bfd passes 533 Python tests and hosted runs 34097112749/34097146588 with Ubuntu/macOS Rust checks and wheel builds. The guarded cascade completed 96/96 frozen validation commands and 14/14 gates; its compact archive has 13 verified entries and an accepted editable six-panel figure | Complete production defaults, broader clean-install portability, repeated isolated timing and real-data validation |
| GitHub/reuse | `9b09bfd` is on main and the working branch; `docs/current_status.md` is the explicit new-window restart record and the compact validation evidence is published | Complete release portability and archive deposition |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
