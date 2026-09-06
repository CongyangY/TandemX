# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Nine inspected six-panel figures, including exact-copy validation, adverse domain shift, fresh localizer validation, failed and passed classifier validations, and ten-library QC; two four-panel single-library QC figures | Isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Ten-library/eight-species input table; complete per-condition and native comparison tables; exact-output 11/111-Mb discovery optimization table | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; two predeclared three-genome conditional held-out matrices; ten full library QCs across eight species; several real comparator diagnostics | Factorial/family-held-out evaluation, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained | Further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | The failed 5701–5703 result was retained. A subsequently frozen depth gate passed fresh 5801–5803, increasing sensitivity 0.594650→0.653635 with unchanged FPR 0.016461 and precision 0.981880→0.983488; localizer recall/precision were 0.980059/0.999617 | Obtain biological under-representation and donor-matched validation; do not generalize the known-catalogue IID result |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim; the adverse split and a later successful frozen-rule split are both retained | Establish reproducible added value against task-matched external tools and biological truth; AI only if independently useful |
| Software quality | Commit 62892a6 passed hosted checks before held-out execution; evidence commit 3b21825 passed hosted Ubuntu/macOS Python, Rust and wheel runs 34047548628/34047562996. A local packed-trace/batched-period change passes 472 Python tests, Rust formatting and Clippy and preserves six/seven real discovery products exactly while reducing one-replay runtime and RSS | Complete hosted executable Rust checks; production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse | `f734004` is on main and the working branch; the failed and passed held-out chains are published with compact hash-checked archives and editable six-panel figures. The exact-output discovery optimization is local | Publish the optimization checkpoint, then complete a clean release, portability checks and archive deposition |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
