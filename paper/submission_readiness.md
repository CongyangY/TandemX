# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Six inspected six-panel figures, including exact-copy validation, adverse domain shift and ten-library QC; two four-panel single-library QC figures | Isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Ten-library/eight-species input table; complete per-condition and native comparison tables | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; two predeclared three-genome conditional held-out matrices; ten full library QCs across eight species; several real comparator diagnostics | Factorial/family-held-out evaluation, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained | Further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | Frozen transparent rule improved fresh exact-copy metrics, but predeclared 5401–5403 domain shift increased sensitivity while worsening FPR 64.30%→68.21% and precision 66.31%→65.58%; localization recall collapsed at 3–5% unit divergence | Develop a divergence-tolerant localizer on independent seeds, then freeze and validate again; biological under-representation and donor-matched validation remain unresolved |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim | Establish reproducible added value and ablations; AI only if independently useful |
| Software quality | Current source passes 442 Python tests; unchanged Rust source passed 15 tests, Clippy and formatting; 0c24b24 passed hosted Ubuntu/macOS CI | Hosted CI for the newest evidence checkpoints; production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse | b85a266 pushed to main and working branch; held-out validation, YSD56 full QC and figure v2 are published; both hosted CI runs passed | YSD56 sampling evidence needs the next verified commit/push; clean release, portability and archive deposition remain pending |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
