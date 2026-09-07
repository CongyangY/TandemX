# TandemX evidence and release programme

Started 2026-09-06 at the user's request. This expands the original toy MVP
scope; it does not retrospectively validate production or biological claims.
Target outlets are Genome Research, Genome Biology and Plant Communications;
journal suitability depends on the resulting evidence, not an acceptance promise.

User clarification, 2026-09-06: expand comparators beyond the initial three and
actively improve speed, memory, accuracy and additional relevant metrics.
Multi-metric superiority is a development objective; it is not an assumed result.
A deficit triggers diagnosis, method improvement and a new frozen evaluation.
Do not redefine the benchmark or omit a strong comparator to hide a deficit.
See [comparator scope and measurement contract](comparator_matrix.md).

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
| Engineering baseline | Full pytest; Rust tests and release build; output validation | Published commit aabaa29 passes 522 local Python tests and hosted runs 34092483185/34092507958. The current guarded cascade candidate passes all 531 Python tests, compileall, Rust formatting and release Clippy; its hosted checks are required before frozen validation may run. The quantified same-machine replay retained byte-identical 5,940-row metrics and all 108 copy-number files while reducing full-matrix time 610.507→349.217 s; this remains an engineering check rather than publication timing |
| Challenge benchmark | Independent seeds; indels, unit divergence, mixed families, multiple arrays, short arrays, hard negatives; strict normalization | Once-only seeds 3101–3103 are complete and archived: TandemX minimum positive array recall/precision 0.985714, zero negative-read calls and RSS ratio 0.455130 versus TideHunter. Promotion failed because nine TRF control runs timed out and TandemX/TideHunter wall-time ratio was 2.457943 (>2.0). On fresh development seed 1201, a guarded gap-free candidate reduced the paired ratio 2.357698→1.989148 with minimum base-union F1 0.997597 and maximum boundary MAE 2.921 bp; a faster unguarded attempt was rejected for a 14.221-bp indel-boundary regression. Frozen validation seed 2201 is untouched and must run once only after the source/config commit passes both hosted workflows |
| Quantification | Known copy number and coverage; finite-read and error bias; catalog ambiguity; calibration/interval coverage | The frozen 6401–6403 run completed 54/54 public commands and passed 9/9 gates: candidate MARE 0.363222 versus 0.408767 baseline, with positive improvement in all three genomes and 661/495/329 improved/equal/worse pairs. Ungated controls were slightly lower at 0.359443. The opt-in CLI gate is implemented; real controls, discovered-catalogue propagation, a calibrated interval and ploidy remain pending |
| Localization/comparison | Known array coordinates; engineered collapse levels; false collapse calls and boundary errors | After the failed alpha-0.5 classifier, a depth-gated v3 rule passed fresh 5801–5803: sensitivity +0.058985, unchanged FPR, precision +0.001608 and all seed guardrails met. Localizer recall/precision were 0.980059/0.999617 with 0/54 absent-family FP. This remains known-catalogue IID simulation; biological truth is absent |
| Real data | Public accession, source checksums, matched sample/assembly, bounded pilot then scaling; known-family recovery | Ten complete libraries across eight species and 262.731 Gb passed full-file QC and deterministic nested sampling. Exact ENA metadata for four retrospective candidates and 12 PacBio WGS runs is hash-archived, but none yet establishes identical donor material across reads and old/new assemblies |
| Comparators | Expanded task-matched methods in comparator_matrix.md | TRF, TideHunter, ULTRA, SRF and TRASH1/2 executed at task-specific stages. A pinned TideCluster 1.21.2 image passed dependency QA and a planted-assembly smoke test with 3/3 families and arrays; plant-scale TideCluster, unitFinder and other applicable assembly/genotyping tests remain pending |
| Innovation | Audited prior art; method ablations; held-out evaluation; AI retained only with reproducible benefit | The transparent multi-k/depth rule failed domain transfer, pooled selection failed seed stability and the seed-robust refinement failed independent FPR/precision gates. A transparent depth gate then passed a second frozen split without refitting. This is a simulation advance; biological validation and external-tool superiority remain pending |
| Reuse | LICENSE; CI; clean source install; toy and real-data workflows; versioned release and hashes | Source install/toy workflows and compact evidence exist; clean release, portability matrix and archive deposition remain pending |
| Paper | Abstract, background, results, methods, discussion; multi-panel figures; tables, supplement, source data and references | Evidence-backed draft plus ten six-panel main figures and source tables exist; Figures 8–10 retain the failed cascade promotion, development trade-offs and passed frozen quantify validation. Scientific completion, final bibliography/metadata and biological validation remain pending |

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
