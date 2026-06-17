# TandemX Known Limitations

## Current Toy-scale Limits

The current MVP:

1. supports FASTA, FASTQ and gzip-compressed sequence inputs, but only for toy-scale datasets;
2. assumes small toy data;
3. uses simple k-mer and shifted-periodicity heuristics;
4. does not use external tandem repeat finders;
5. does not perform read mapping;
6. does not infer higher-order repeat structure;
7. does not model ploidy or subgenomes;
8. does not provide experimentally calibrated FISH probe prediction.

Anti-hardcoding and randomized toy workflow tests now cover non-default repeat lengths and fixed random seeds. This checks that the toy MVP is not narrowly tied to simulator defaults and that output schemas remain valid across small controlled cases. It does not establish performance on real plant repeats, related repeat families, noisy real reads, polyploid genomes, or chromosome-scale assemblies.

## Not Suitable for Real Large Genomes Yet

The current code should not be used for production analysis of wheat, rye, barley, oat, maize or other 7-20 Gb plant genomes.

Reasons:

1. read processing is streaming at the parser layer but the algorithms are not optimized for large inputs;
2. k-mer counting is in-memory Python;
3. assembly scanning is simple and not indexed;
4. monomer clustering is period-based and not robust to related repeat families;
5. copy-number calibration depends on user-provided or rough haploid depth;
6. localization uses k-mer evidence, not alignments;
7. probe specificity is a heuristic, not a validated hybridization model.

## Claims Not Supported

Do not claim that TandemX can:

1. fully resolve megabase-scale satellite arrays;
2. precisely locate every repeat copy;
3. completely assemble satellite arrays from reads;
4. guarantee FISH probe success;
5. outperform TRF, TideHunter, TRASH or RepeatExplorer2/TAREAN on real genomes.

## Current Deferred Command

`tandemx compare` remains deferred. The implemented assembly-vs-read comparison currently lives in `tandemx locate`.
