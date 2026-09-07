# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Twelve inspected six-panel figures, including exact-copy validation, adverse domain shift, failed/passed classifier validations, ten-library QC, failed cascade promotion, quantify development and frozen quantify validation; two four-panel single-library QC figures | Isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Ten-library/eight-species input table; complete per-condition and native comparison tables; exact-output 11/111-Mb discovery optimization table | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Six independent factorial quantify genomes/54 conditions; three other predeclared three-genome held-out matrices including the consumed cascade split; ten full library QCs across eight species; several real comparator diagnostics | Family/process-held-out evaluation, broader species/material replication, true biological recall and multiple isolated large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained. TideCluster 1.21.2 now has a pinned runnable image and successful planted-truth assembly smoke with all dependency/resource evidence archived | Plant-scale TideCluster and further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | The 108-run development ablation and exact-output speed replay completed. On frozen seeds 6401–6403, 54/54 public commands and 9/9 gates passed: MARE 0.408767→0.363222 with positive improvement in every genome and 661/495/329 improved/equal/worse family pairs. Ungated controls were slightly lower at 0.359443. The rule is an opt-in public CLI mode | Real control specificity, discovered-catalogue propagation, calibrated sampling intervals, catalogue uncertainty and ploidy remain |
| Assembly interpretation | The failed 5701–5703 result was retained. A subsequently frozen depth gate passed fresh 5801–5803, increasing sensitivity 0.594650→0.653635 with unchanged FPR 0.016461 and precision 0.981880→0.983488; localizer recall/precision were 0.980059/0.999617 | Obtain biological under-representation and donor-matched validation; do not generalize the known-catalogue IID result |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim; the adverse split and a later successful frozen-rule split are both retained | Establish reproducible added value against task-matched external tools and biological truth; AI only if independently useful |
| Software quality | Frozen commit 9eda196 passed 514 Python tests and hosted runs 34088623796/34088640260. Current source passes 520 Python tests, Rust format/Clippy, wheel build and an isolated full toy workflow. The validation has a 93-payload input audit, 67-entry compact archive and 54 execution-artifact rows. Local Rust tests compile but hosted CI must execute them because conda lacks dynamic libpython | Publish the opt-in CLI and evidence, then complete production defaults, broader clean-install portability and repeated isolated timing |
| GitHub/reuse | `9eda196` is on main and the working branch; the complete validation and opt-in CLI update are prepared locally. `docs/current_status.md` is the explicit new-window restart record | Commit/push the validation evidence and CLI, then complete release portability and archive deposition |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
