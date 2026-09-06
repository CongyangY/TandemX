# Comparator scope and measurement contract

Updated 2026-09-06. This is a prospective evaluation contract, not a table of
completed comparisons. The aim is to establish reproducible advantages against
strong, task-matched methods over several metrics. A method is not declared
superior simply because its combined workflow has more commands.

## Review-relevant methods

| Method | Relevant endpoint | Planned role and boundary | Actual state |
| --- | --- | --- | --- |
| TRF | De novo intervals and repeat period | Established comparator; raw duplicate calls and standardized interval-union metrics both required | Development benchmark executed |
| TideHunter | Read-local repeat consensus and intervals | Strong long-read comparator; retain recommended settings and matched period range | Development benchmark executed |
| ULTRA | Genomic tandem-repeat annotation | Indel-aware probabilistic baseline; default and documented tuning configurations, including tuning cost | v1.2.2 built; default and tuned ten-read pilots completed |
| TRASH | Assembly tandem arrays and monomers | Plant/large-genome relevance; matched assembly or synthetic sequence endpoints | Source reviewed; not yet run |
| SRF | De novo satellite families, units/HORs and abundance | Essential read-first novelty comparator; include k-mer construction and abundance mapping in pipeline cost | Source reviewed; not yet run |
| RepeatExplorer2/TAREAN | Satellite clustering and abundance from short reads | Separate matched short-read evaluation if suitable reads exist; no artificial failure on unsupported HiFi input | Not yet run |
| HiCAT/HiCAT-human | Monomer/HOR organization | Prior-art and specialized analysis; human pretrained classifiers do not supply a fair plant de novo baseline | Applicability review pending |
| StringDecomposer/NCRF | Known-motif decomposition | Candidate targeted annotation baselines, with identical supplied motif catalogues; not de novo family discovery | Applicability review pending |
| Chorus2 | Genome-based oligo-FISH probe specificity | Review specificity and oligo constraints; unique chromosome-painting oligos and repeated satellite probes have different objectives | Official source reviewed; task-matched experiment pending |

Add or retire a comparator only with a documented scientific or reproducibility
reason. Missing dependencies, installation failures, resource limits and version
incompatibilities stay in the audit table. They are not zero accuracy values.
Do not use human-trained centromere labels as plant ground truth.

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
- [TRASH official implementation](https://github.com/vlothec/TRASH)
- [TideHunter official implementation](https://github.com/Xinglab/TideHunter)
- [Chorus2 official implementation](https://github.com/zhangtaolab/Chorus2)
