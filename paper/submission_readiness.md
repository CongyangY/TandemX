# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Four inspected six-panel figures for quantification, intervals, cross-cohort QC and real-input diagnostics; two four-panel single-library QC figures | Final method, isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Nine-library/seven-species input table; complete per-condition and native comparison tables | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; one predeclared three-genome conditional held-out matrix; nine full library QCs across seven species; several real comparator diagnostics | Factorial/family-held-out evaluation, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained | Further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | Engineered retained-copy matrix completed once on seeds 5101–5103: sensitivity 85.60%, false-positive rate 4.94%, precision 96.30%; exact-copy localization recall 36/36; a transparent rule improved all three confusion metrics during development | Frozen rule still requires one-time validation on untouched seeds 5201–5203; divergent/fragmented-array, biological under-representation and donor-matched validation remain pending |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim | Establish reproducible added value and ablations; AI only if independently useful |
| Software quality | Current source passes 426 Python tests; unchanged f16596b Rust source passed 15 tests, Clippy and formatting | Production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse | f16596b pushed to main and working branch; both hosted CI runs passed | Current clustering replay archiver and resulting evidence need the next verified commit/push |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
