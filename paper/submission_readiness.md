# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Existing six-panel quantification and interval figures; two four-panel QC figures | Fix interval-panel-A legend; final method, comparative scaling and biological figures; visual/source QA of final package |
| Tables and source data | Seven-species input table; complete per-condition and native comparison tables | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; seven full library QCs; several real comparator diagnostics | Held-out genomes/families, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; native alternatives/failed parsing retained | Further applicable tools, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | References and source issues audited; preliminary comparison modules | Independent engineered-collapse and biological under-representation validation; donor matching |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim | Establish reproducible added value and ablations; AI only if independently useful |
| Software quality |389 local Python tests and13 native tests at latest verified checkpoints | Production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse |d2aa4d6 pushed to main and working branch; earlier CI checkpoints verified | New local mapping/manuscript/evidence changes need commit/push; latest hosted CI status needs confirmation |

On 6 September 2026 automatic approval rejected the T7 figure-v2 render because
the review service reached its usage limit. The render did not execute. Existing
jobs continue; read-only checks and local source/evidence preparation are
unaffected. Do not bypass that refusal. User authorization for project Git/T7
work already exists, but approval-service availability must be restored before
dependent actions proceed. This temporary service issue is separate from the
substantial scientific and release work that still remains.
