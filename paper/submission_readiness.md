# Submission readiness audit

The manuscript is an evidence-backed development draft. The requested mature
software, biological benchmark, complete paper and reusable public release
have **not** been achieved.

| Requirement | Evidence now | Remaining work |
| --- | --- | --- |
| Complete article sections | `manuscript.md` contains abstract, background, results, discussion, methods, legends and references | Scientific completion, editorial revision, full bibliography, author/funding/declaration metadata |
| Multi-panel figures | Eight inspected six-panel figures, including exact-copy validation, adverse domain shift, fresh localizer validation, failed classifier validation and ten-library QC; two four-panel single-library QC figures | Isolated repeated comparative scaling and biological-validation figures; visual/source QA of final package |
| Tables and source data | Ten-library/eight-species input table; complete per-condition and native comparison tables | Final numbered standalone tables, workbook/package and manuscript cross-reference audit |
| Real and simulated testing | Three independent factorial genomes/27 conditions; two predeclared three-genome conditional held-out matrices; ten full library QCs across eight species; several real comparator diagnostics | Factorial/family-held-out evaluation, broader species/material replication, hard negatives, true biological recall and multiple large-input runs |
| Fair comparators | TRF/TideHunter executed; TRASH1/2 task-specific evaluations; real source-query recovery and threshold sensitivity; native alternatives/failed parsing retained | Further applicable tools, donor-matched truth, controlled sensitivity settings and isolated same-platform resource comparisons |
| Quantification | Improved conditional point estimation; joint sampling coverage/missingness | Nuclear-depth normalization, real background specificity, catalogue uncertainty, ploidy and held-out calibration |
| Assembly interpretation | Bounded bridging retained strong aggregate localization on fresh 5701–5703. The double-frozen seed-robust classifier increased sensitivity 0.628258→0.716049 but failed held-out FPR and precision gates; all 16 added false positives were at 1× | Treat 5701–5703 as consumed development evidence, freeze any low-depth fallback before new seeds, and obtain biological under-representation and donor-matched validation |
| Probe performance | Computational prioritization framework | Specificity calibration and independent published/new experimental concordance; no fabricated FISH rate |
| Novelty/AI | Explicit prior-art comparison; no unsupported AI claim; independent classifier failure retained | Establish reproducible added value and ablations on a new frozen split; AI only if independently useful |
| Software quality | Current source passes 462 Python tests; unchanged Rust source passed 15 tests, Clippy and formatting; 300e48d passed hosted Ubuntu/macOS CI before held-out execution | Verify the new compact archiver/figure/depth-gate tests in hosted CI; production defaults, clean install/release verification, workflow portability and complete command validation |
| GitHub/reuse | 300e48d is pushed to main and the working branch; double-freeze guards passed hosted CI; the failed held-out result now has a compact hash-checked archive and editable six-panel figure in the working tree | Publish this evidence checkpoint, then complete a clean release, portability checks and archive deposition |

On 6 September 2026 an initial T7 figure-v2 render was rejected because the
automatic review service reached its usage limit. Later Git operations passed
review. After a read-only diff established the layout-only change and verified
that the destination was new, a retry through the same approval channel passed
and version2 was rendered and inspected. The transient service issue is resolved;
the substantial scientific and release work above still remains.
