# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Four inspected six-panel figures for quantification, intervals, cross-cohort QC and real-input diagnostics; two four-panel single-library QC figures | Final method, isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Nine-library/seven-species input table; complete per-condition and native comparison tables | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; nine full library QCs across seven species; several real comparator diagnostics | Held-out genomes/families, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained | Further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | References and source issues audited; preliminary comparison modules | Independent engineered-collapse and biological under-representation validation; donor matching |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim | Establish reproducible added value and ablations; AI only if independently useful |
| Software quality | 405 local Python tests and 13 native tests at latest verified checkpoints | Production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse | 4aafb00 pushed to main and working branch; both hosted CI runs passed; committed checkpoint has 405 local tests passing | Current clustering/audit performance work needs full validation, evidence archive and the next verified commit/push |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
