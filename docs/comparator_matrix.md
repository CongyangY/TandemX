# Comparator scope and measurement contract

Updated 2026-09-07. This is an evaluation contract with explicit execution
status, not a claim that every listed comparison is complete. The aim is to establish reproducible advantages against
strong, task-matched methods over several metrics. A method is not declared
superior simply because its combined workflow has more commands.

## Review-relevant methods

| Method | Relevant endpoint | Planned role and boundary | Actual state |
| --- | --- | --- | --- |
| TRF | De novo intervals and repeat period | Established comparator; raw duplicate calls and standardized interval-union metrics both required | Development benchmark executed; first cascade held-out retained nine low-complexity timeouts |
| TideHunter | Read-local repeat consensus and intervals | Strong long-read comparator; retain recommended settings and matched period range | Development and first held-out comparisons executed; fresh guarded-cascade validation completed 96/96 commands and 14/14 gates |
| TideCluster | Assembly tandem-repeat detection and family clustering | Tier-1 assembly comparator; run its TideHunter detection and MMseqs2/BLAST clustering as one end-to-end method, then score intervals and sequence-supported families independently | Reproducible 1.21.2 linux/amd64 image and planted-truth smoke passed. Nested 10/100-Mb MorexV3 reference-window runs completed with external-stage resources and family-aware interval provenance. The frozen 1-Gb same-host gate formally refused execution after the 100-Mb clustering peak used 96.836% of the 8.217-GB Docker limit; accuracy remains under frozen factorial evaluation, and whole-genome context remains unavailable |
| ULTRA | Genomic tandem-repeat annotation | Indel-aware probabilistic baseline; default and documented tuning configurations, including tuning cost | v1.2.2 built; default and tuned ten-read pilots completed |
| TRASH | Assembly tandem arrays and monomers | Plant/large-genome relevance; matched assembly or synthetic sequence endpoints | Pinned Linux ARM64 offline author positive-control completed;355 monomer sequences require a recorded -1 bp coordinate offset for exact source agreement; no plant-scale comparison yet |
| TRASH2 | Assembly tandem arrays and monomers | Author-designated early-development successor, evaluated separately from TRASH1 | Default 10-Mb synthetic run completed: 55/55 families and unit-derived arrays; native coarse-window and unit endpoints retained. Two author-control runs differ in unit/consensus outputs; no real plant assembly benchmark yet |
| SRF | De novo satellite families, units/HORs and abundance | Essential read-first novelty comparator; include k-mer construction and abundance mapping in pipeline cost | Four-run pilot and 32-run development suite complete; native failures and explicit empty-count guard outcomes both retained |
| RepeatExplorer2/TAREAN | Satellite clustering and abundance from short reads | Separate matched short-read evaluation if suitable reads exist; no artificial failure on unsupported HiFi input | Not yet run |
| HiCAT/HiCAT-human | Monomer/HOR organization | Prior-art and specialized analysis; human pretrained classifiers do not supply a fair plant de novo baseline | Applicability review pending |
| StringDecomposer/NCRF | Known-motif decomposition | Candidate targeted annotation baselines, with identical supplied motif catalogues; not de novo family discovery | Applicability review pending |
| TandemTools | Long-read mapping, polishing and quality assessment for assembled extra-long tandem repeats | Apply to engineered assembly errors or candidate arrays after target coordinates exist; it is not a de novo monomer-discovery baseline | Applicability review pending |
| RaMA | Pairwise alignment of centromere assemblies and higher-order-repeat structure | Relevant to assembly-structure concordance and alignment resource use; it does not discover read-level families or estimate whole-sample copy number | Applicability review pending |
| TRsv | Reference-anchored tandem-repeat CNV plus SV/indel calling from long-read alignments | Non-human mode requires a user-supplied repeat BED/unit catalogue; compare only on known-locus CNV tests, not de novo satellite discovery | Applicability review pending |
| Chorus2 | Genome-based oligo-FISH probe specificity | Review specificity and oligo constraints; unique chromosome-painting oligos and repeated satellite probes have different objectives | Official source reviewed; task-matched experiment pending |
| unitFinder | Plant centromeric monomer discovery/decomposition in assemblies | New Genome Biology soybean study; iterative TRF/nucmer workflow, supplied chromosome and de novo/reference-assisted modes must be matched separately | Primary methods reviewed; implementation/dependency audit and experiment pending |
| CentIER | Assembly centromere-region prediction using repeats, retrotransposons and k-mer features | Plant Communications method; optional annotation/Hi-C inputs require matched evidence. Region prediction is distinct from read-first monomer recovery | Official source reviewed; dependency/build and matched-region experiment pending |

Add or retire a comparator only with a documented scientific or reproducibility
reason. Missing dependencies, installation failures, resource limits and version
incompatibilities stay in the audit table. They are not zero accuracy values.
Do not use human-trained centromere labels as plant ground truth.

The SRF primary paper explicitly reports recovery of maize CentC with k=101
after failure at k=151. Its HiFi high-count filter also depends on coverage.
The existing k=151/ci20/ci100 toy workflow is therefore a development condition,
not a sufficient tuned comparator. Include frozen k/count sensitivity and HOR
decomposition endpoints before attributing misses to the method generally.

The current unitFinder README defines a chromosome-by-chromosome assembly
workflow built around iterative TRF calls, nucmer grouping (`-c 10 -l 10`),
cross-chromosome merging and a separate reference-assisted `--cen` pass. Its
published example also includes manual concatenation, script-format edits and
Clustal Omega steps. A reproducible comparison must therefore freeze the exact
repository state, containerize every external dependency and scripted edit, and
measure de novo and supplied-reference modes separately. Dependency and manual
post-processing time cannot be omitted from end-to-end resource results.

The [real/simulated cohort and QC programme](cohort_and_qc.md) governs enrollment,
data scale, true independent units, reference quality and sampling.

## Metrics to freeze before the publication test

1. Detection: array and family recall/precision/F1, read-level false-positive
   rate, base-level interval precision/recall, duplicate rate, fragmentation,
   harmonics, and failures. Ground-truth-family recall requires sequence evidence.
2. Structure: circular sequence identity/edit distance, monomer length error,
   start/end error, interval IoU, multiple-array recovery and partial-array behavior.
3. Quantification: absolute and relative copy-number/abundance error, systematic
   bias by coverage and error rate, uncertainty interval coverage/width, family
   ambiguity and contamination sensitivity.
4. Assembly comparison: collapse sensitivity/specificity, precision-recall curves,
   false collapse calls under normal assembly representation, and calibration.
5. Probe candidates: target coverage, off-target genomic abundance, family
   cross-reactivity, GC/Tm/secondary-structure constraints, and independently
   documented FISH validation where available. In-silico scores alone are not
   successful experiments.
6. Resources: end-to-end wall time, CPU time, aggregate peak memory including
   child processes, temporary disk usage, throughput, thread scaling and input
   size scaling. Also report per-stage costs; include index/model preparation,
   downloaded model sizes and inference costs for any AI component.
7. Reliability: deterministic outputs, graceful empty results, truncated or
   invalid-input handling, resumability, install portability and reproducibility.

## Fairness and claims

- Freeze development/validation/test families and biological samples. Do not
  tune on held-out results. Related synthetic families must not leak across splits.
- Run all methods on the same observable input, hardware and resource allocation.
  Use single-thread and comparable multi-thread experiments. Record default and
  equally budgeted tuned settings separately; count tuning time explicitly.
- Pair comparisons by dataset. Timing repeats quantify machine variability;
  they are not biological replicates. Use confidence intervals over independent
  datasets/samples; publish all per-dataset rows and failed runs.
- Report raw outputs and a tool-independent normalization (e.g. interval union)
  to reveal whether apparent precision gains are only duplicate suppression.
- Predeclare meaningful superiority and noninferiority margins after development
  variability is measured, before test. Report trade-offs and Pareto fronts;
  a composite score cannot hide a substantial accuracy or resource deficit.
- An all-metric advantage is claimed only over the named methods, dataset domain,
  settings and measured uncertainty actually tested. Synthetic development gains
  do not establish production readiness or universal superiority.

## Primary sources

- [ULTRA official implementation and tuning](https://github.com/TravisWheelerLab/ULTRA)
- [ULTRA 2024 paper](https://doi.org/10.1093/bioadv/vbae149)
- [SRF official implementation](https://github.com/lh3/srf)
- [SRF Genome Research paper (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760446/)
- [TRASH official implementation](https://github.com/vlothec/TRASH)
- [TideHunter official implementation](https://github.com/Xinglab/TideHunter)
- [TideCluster official implementation](https://github.com/kavonrtep/TideCluster)
- [Chorus2 official implementation](https://github.com/zhangtaolab/Chorus2)
- [TandemTools paper](https://doi.org/10.1093/bioinformatics/btaa440) and [official implementation](https://github.com/ablab/TandemTools)
- [RaMA Genome Research paper](https://doi.org/10.1101/gr.279763.124) and [official implementation](https://github.com/pinglu-zhang/RaMA)
- [TRsv Genome Biology paper](https://doi.org/10.1186/s13059-025-03718-z) and [official implementation](https://github.com/stat-lab/TRsv)
- [unitFinder primary Genome Biology study](https://doi.org/10.1186/s13059-025-03924-9)
- [unitFinder official implementation](https://github.com/HuangYicheng-Bio/unitFinder)
- [CentIER official implementation](https://github.com/simon19891216/CentIER)
