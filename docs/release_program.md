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
| Engineering baseline | Full pytest; Rust tests and release build; output validation | Current source passes all 431 Python tests; unchanged Rust source passed 15 tests, fmt and clippy; both b85a266 hosted CI runs passed on Ubuntu and macOS |
| Challenge benchmark | Independent seeds; indels, unit divergence, mixed families, multiple arrays, short arrays, hard negatives; strict normalization | Development and separate validation seeds executed; related-family merge failure retained; held-out publication seeds remain unused |
| Quantification | Known copy number and coverage; finite-read and error bias; catalog ambiguity; calibration/interval coverage | Three-genome/27-condition multi-k and joint-read experiments complete; 20× conditional interval coverage 94.58% with explicit missingness; held-out/real calibration pending |
| Localization/comparison | Known array coordinates; engineered collapse levels; false collapse calls and boundary errors | Fresh 5201–5203 validation improved TP/FN/FP/TN from 198/45/14/148 to 208/35/12/150 with a frozen rule; divergent/fragmented arrays and biological truth remain pending |
| Real data | Public accession, source checksums, matched sample/assembly, bounded pilot then scaling; known-family recovery | Ten complete libraries across eight species and 262.731 Gb passed full-file QC and deterministic nested sampling; donor-matched truth remains pending |
| Comparators | Expanded task-matched methods in comparator_matrix.md | TRF, TideHunter, ULTRA, SRF and TRASH1/2 executed at task-specific stages; unitFinder and other applicable assembly/genotyping tools remain pending |
| Innovation | Audited prior art; method ablations; held-out evaluation; AI retained only with reproducible benefit | Frozen transparent multi-k/depth rule improved sensitivity, false-positive rate and precision on all three fresh held-out genomes; current evidence is narrow exact-copy simulation, with domain-shift and biological validation pending |
| Reuse | LICENSE; CI; clean source install; toy and real-data workflows; versioned release and hashes | Source install/toy workflows and compact evidence exist; clean release, portability matrix and archive deposition remain pending |
| Paper | Abstract, background, results, methods, discussion; multi-panel figures; tables, supplement, source data and references | Evidence-backed draft plus four six-panel figures and source tables exist; scientific completion, final bibliography/metadata and biological validation remain pending |

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
