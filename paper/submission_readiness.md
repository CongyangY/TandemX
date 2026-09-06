# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Seven inspected six-panel figures, including exact-copy validation, adverse domain shift, fresh localizer validation and ten-library QC; two four-panel single-library QC figures | Isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Ten-library/eight-species input table; complete per-condition and native comparison tables | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; two predeclared three-genome conditional held-out matrices; ten full library QCs across eight species; several real comparator diagnostics | Factorial/family-held-out evaluation, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained | Further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | Bounded bridging passed fresh 5501–5503 localization gates. Pooled classifier development failed seed stability; seed-robust refinement selected alpha 0.5/threshold 0.5 and improved sensitivity, FPR and precision within every 5601–5603 seed | Commit and CI-check the double-frozen held-out guard before touching 5701–5703; biological under-representation and donor-matched validation remain unresolved |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim | Establish reproducible added value and ablations; AI only if independently useful |
| Software quality | Current source passes 456 Python tests; unchanged Rust source passed 15 tests, Clippy and formatting; ebfaaf5 passed hosted Ubuntu/macOS CI | Hosted CI for the double-frozen held-out guard; production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse | ebfaaf5 pushed to main and working branch; seed-robust development selection, raw multi-k checkpointing, held-out localizer validation, YSD56 full QC and figure evidence are published; both hosted CI runs passed | Double-frozen held-out config and evaluator need the next verified commit/push; clean release, portability and archive deposition remain pending |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
