# TandemX evidence and release programme

Started 2026-09-06 at the user's request. This expands the original toy MVP
scope; it does not retrospectively validate production or biological claims.
Target outlets are Genome Research, Genome Biology and Plant Communications;
journal suitability depends on the resulting evidence, not an acceptance promise.

The 2026-09-06 comparator and optimization expansion is superseded by the
2026-09-09 feature freeze. Do not add algorithms, comparators, convenience
features or low-yield optimizations, and do not retry `unitFinder`. The final
method-science gate is the frozen Ey15-2/Macadamia orthogonal abundance
validation in [orthogonal_abundance_validation.md](orthogonal_abundance_validation.md).
After it passes, work moves to manuscript, figure, Bioconda, Zenodo and release
packaging. If it demonstrates quantification bias, only that bounded bias may be
repaired before the gate is repeated.

## Storage and provenance

Source checkout: `/Users/ycy/Codex/Sofw/TandemX`. Data and result root on this
machine: `/Volumes/T7/Codex/TandemX`. Keep code and environment on the internal
APFS volume; keep large data, scratch and releases on T7. A new user can choose
any writable data root. Never use an assembly-alignment BAM as raw HiFi input.
Retain originals until copied files pass checksum verification. Keep downloaded
data, compiled binaries and caches out of Git. Archive source revisions and
machine-readable commands alongside results.

## Acceptance gates

| Gate | Required evidence | Current evidence |
| --- | --- | --- |
| Engineering baseline | Full pytest; Rust tests and release build; output validation | Validation freeze 054b935 passed 531 local Python tests and hosted checks before seed 2201 was used. Published archive commit 9b09bfd passes 533 local Python tests, compileall and hosted runs 34097112749/34097146588 on Ubuntu/macOS with Rust checks and wheel builds. The quantified same-machine replay retained byte-identical 5,940-row metrics and all 108 copy-number files while reducing full-matrix time 610.507→349.217 s; this remains an engineering check rather than publication timing |
| Challenge benchmark | Independent seeds; indels, unit divergence, mixed families, multiple arrays, short arrays, hard negatives; strict normalization | Once-only seeds 3101–3103 are complete and archived: the first promotion failed because nine TRF control runs timed out and the TandemX/TideHunter wall-time ratio was 2.457943 (>2.0). On fresh development seed 1201, a guarded gap-free candidate reduced the paired ratio 2.357698→1.989148 with minimum base-union F1 0.997597 and maximum boundary MAE 2.921 bp; a faster unguarded attempt was rejected for a 14.221-bp indel-boundary regression. After source/config freeze and CI, validation seed 2201 ran once: 96/96 commands and 14/14 gates passed, with runtime ratio 1.977877, RSS ratio 0.396864, minimum array recall/precision 1.0, minimum base-union F1 0.997445, maximum boundary MAE 2.35 bp and zero TandemX negative calls. This remains synthetic distribution-level evidence; some scenario runtime ratios exceeded 4 and real-data superiority is pending |
| Quantification | Orthogonal read support for frozen real families; finite-read and normalization bias; independent exact-count reproduction | Simulations passed their frozen gates. Ey15 PCR-free Illumina `ERR8666067`, Macadamia Illumina `SRR11191912` and Macadamia PromethION `SRR11191910` are source-frozen but not yet downloaded or analyzed. The gate separates read--assembly deficit from binary collapse and includes the Macadamia 616/653/738/780-Mb normalization sensitivity |
| Localization/comparison | Known array coordinates; engineered collapse levels; false collapse calls and boundary errors | Simulated gates passed. Ey15 and Macadamia old/new comparisons are high-quality assembly-reference proxies, not copy-number truth. Their binary concordance is retained separately from continuous magnitude; orthogonal reads must now distinguish predictor bias from residual assembly collapse |
| Real data | Public accession, source checksums, matched sample/assembly, bounded pilot then scaling; known-family recovery | Ten complete libraries across eight species and 262.731 Gb passed prior QC. The final source audit found same-tissue PCR-free Illumina for Ey15 and same-old-study BioSample Illumina plus PromethION for Macadamia. Complete QC and family-level abundance results remain the final biological-method blocker |
| Comparators | Frozen task-matched methods in comparator_matrix.md | TRF, TideHunter, ULTRA, SRF, TRASH1/2 and TideCluster evidence is final for this version. No new comparator is permitted. `unitFinder` is permanently stopped; its incomplete smoke remains provenance rather than runtime or accuracy evidence |
| Innovation | Audited prior art; method ablations; held-out evaluation; AI retained only with reproducible benefit | The transparent multi-k/depth rule failed domain transfer, pooled selection failed seed stability and the seed-robust refinement failed independent FPR/precision gates. A transparent depth gate then passed a second frozen split without refitting. The catalogue now emits a pairwise candidate period-multiple graph without forcing 171/342/684-style alternatives into one family or claiming validated HOR structure. These are simulation/software advances; biological validation and external-tool superiority remain pending |
| Reuse | LICENSE; CI; clean source install; toy and real-data workflows; versioned release and hashes | Source install/toy workflows and compact evidence exist. A strict `tandemx import tidehunter` route validates native calls against source reads, retains provenance and produces standard downstream catalogues; clean release, portability matrix and archive deposition remain pending |
| Paper | Abstract, background, results, methods, discussion; multi-panel figures; tables, supplement, source data and references | Evidence-backed draft and figures exist. Language now treats the newer assembly as a measurement, continuous results as read--assembly deficit/estimated under-representation, and binary assembly-proxy classification separately from magnitude. Orthogonal abundance validation, then final bibliography/metadata and release packaging, remain pending |

## Benchmark discipline

Distinguish **read detection**, **array recovery**, **period recovery** and
**sequence-supported family recovery**. They have different denominators.
Do not call read-level period recall family recall. Match arrays one-to-one
by read, interval overlap and period; unmatched duplicates count against
array precision. Score every tool on the same observable input. A crash,
missing output, malformed row or incompatible schema is a failed run, not a
zero-recall measurement. Publish failures and exclusion reasons.

Use separate development and held-out seeds/families. Tuning on a benchmark
makes it development data thereafter. Repeated timing runs are not independent
biological replicates. Report uncertainty across independently generated
datasets or biological samples, retaining per-dataset metrics. Fix threads,
record versions, executable/source hashes, wall time and direct-child peak RSS.
Do not equate RSS with total pipeline memory for tools spawning children.

## Scientific differentiation to test

SRF already supports satellite discovery and abundance estimation from accurate
long reads. Local periodic k-mers alone therefore do not establish novelty.
Test whether catalog ambiguity, calibrated read evidence and assembly evidence
can be combined into better-supported abundance and under-representation calls.
Learned confidence is an optional experiment: compare a transparent baseline,
calibration model and no-calibration ablation on family-held-out data, with
calibration error and domain-shift reporting. Synthetic-only improvement is
insufficient to claim a general plant model. FISH remains prioritization until
independent experimental or published probe evidence supports validation.

Primary starting references (consulted 2026-09-06):

- [TideHunter](https://doi.org/10.1093/bioinformatics/btz376)
- [TRASH](https://doi.org/10.1093/bioinformatics/btad308)
- [SRF source and workflow](https://github.com/lh3/srf)
- [Arabidopsis HiFi/assembly data study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9757041/)
- [Col-CEN HiFi study accession PRJEB46164](https://www.ebi.ac.uk/ena/browser/view/PRJEB46164)

Each gate must be updated with actual artifacts and checks before being marked
complete. Intermediate manuscripts are research drafts, not submission-ready
papers. Published FISH observations and new experiments must be distinguished.
