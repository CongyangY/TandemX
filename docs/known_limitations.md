# TandemX Known Limitations

## Current Toy-scale Limits

The current MVP:

1. supports FASTA, FASTQ and gzip-compressed sequence inputs; engineering benchmarks include bounded real-read subsets, but biological validation remains toy-scale;
2. performs de novo discovery from reads in `tandemx discover`;
3. assumes small toy data;
4. uses a bounded k-mer spacing prefilter and local period refinement with Python and optional Rust backends;
5. does not use external tandem repeat finders;
6. does not perform read mapping;
7. flags possible higher-order or partial family relationships but does not infer higher-order repeat structure;
8. does not model ploidy or subgenomes;
9. does not provide experimentally calibrated FISH probe prediction.

The downstream `--catalog` input reuses the de novo repeat catalog produced by `tandemx discover`. It is not a requirement that users already know the repeat sequence before running TandemX.

Known-repeat annotation is post hoc. `tandemx annotate-repeats` can compare
discovered monomers to a known-repeat FASTA after discovery, but that library is
not an input to `tandemx discover` and does not make discovery template-guided.

Anti-hardcoding and randomized toy workflow tests now cover non-default repeat lengths and fixed random seeds. This checks that the toy MVP is not narrowly tied to simulator defaults and that output schemas remain valid across small controlled cases. It does not establish performance on real plant repeats, related repeat families, noisy real reads, polyploid genomes, or chromosome-scale assemblies.

## Not Suitable for Real Large Genomes Yet

The current code should not be used for production analysis of wheat, rye, barley, oat, maize or other 7-20 Gb plant genomes.

Reasons:

1. read processing is streaming at the parser layer but the algorithms are not optimized for large inputs;
2. read-local seed spacing can use Rust, but global k-mer counting and downstream large-scale algorithms are not production backends;
3. assembly scanning is simple and not indexed;
4. monomer clustering is period-based and not robust to related repeat families;
5. copy-number calibration depends on user-provided or rough haploid depth;
6. localization uses k-mer evidence, not alignments;
7. probe specificity is a heuristic, not a validated hybridization model.

Optional family collapse is disabled by default. `--collapse-redundant-families`
collapses only `likely_redundant` relationships and keeps an audit table.
`possible_higher_order_or_partial` relationships are retained for user review
and should not be described as definitively redundant or definitively
higher-order without assembly, read-support, known-repeat annotation or
experimental evidence.

Known-repeat annotation and optional collapse improve interpretability, but they
do not replace assembly localization, FISH validation, or literature-supported
repeat biology interpretation.

An optional guided mode using user-supplied known-repeat FASTA files may be added later, but it is not the default MVP workflow.

The synthetic benchmark harness can measure tiny and manual scale runs, but it does not make the MVP suitable for real 7-20 Gb production analysis. Real data should be limited to pilot subsets until chunking, resumable execution, memory reporting and external benchmarks are implemented.

Discover now provides incremental candidates and progress logs, but it remains single-process. `--chunk-size` does not yet provide checkpoint/resume behavior. The Python backend is suitable for toy and staged subset tests; the Rust backend makes larger real-read pilots practical but does not address clustering memory, restartability, full-dataset resource reporting or downstream scaling.

`tandemx run --resume` operates between commands only. It validates expected output files before skipping a step, but it does not resume a partially processed FASTQ, fingerprint inputs, or detect that an upstream file changed after a downstream result was written. `--threads` is recorded and currently must be 1.

Quantify's in-process counters retain only catalogue-derived diagnostic k-mers, which avoids all-read distinct-k-mer memory growth. This is appropriate for toy and bounded pilot runs; it is not a substitute for KMC, meryl, Jellyfish, or production-scale coverage calibration.

## Claims Not Supported

Do not claim that TandemX can:

1. fully resolve megabase-scale satellite arrays;
2. precisely locate every repeat copy;
3. completely assemble satellite arrays from reads;
4. guarantee FISH probe success;
5. outperform TRF, TideHunter, TRASH or RepeatExplorer2/TAREAN on real genomes.

## Current Deferred Command

`tandemx compare` remains deferred. The implemented assembly-vs-read comparison currently lives in `tandemx locate`.
