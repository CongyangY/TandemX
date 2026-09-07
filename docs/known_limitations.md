# TandemX Known Limitations

## Current Toy-scale Limits

The current MVP:

1. supports FASTA, FASTQ and gzip-compressed sequence inputs; engineering benchmarks include bounded real-read subsets, but biological validation remains toy-scale;
2. performs de novo discovery from reads in `tandemx discover`;
3. assumes small toy data;
4. uses a bounded k-mer spacing prefilter and local period refinement with Python and optional Rust backends;
5. keeps native discovery independent of external tandem repeat finders, while
   an explicit provenance-labelled importer can ingest TideHunter `-f 2` calls;
6. does not perform read mapping;
7. emits a directional candidate period-multiple graph for higher-order/partial
   review, but does not infer or validate higher-order repeat structure;
8. does not model ploidy or subgenomes;
9. does not provide experimentally calibrated FISH probe prediction.

The downstream `--catalog` input normally reuses the de novo repeat catalog
produced by `tandemx discover`. It can also use a catalogue from
`tandemx import tidehunter`, whose external detector provenance remains in the
candidate, family and import-audit outputs. It is not a requirement that users
already know the repeat sequence before running TandemX.

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
4. sequence-aware monomer clustering remains insufficiently validated across related natural repeat families;
5. copy-number calibration depends on user-provided or rough haploid depth;
6. localization uses k-mer evidence, not alignments;
7. probe specificity is a heuristic, not a validated hybridization model.

Optional family collapse is disabled by default. `--collapse-redundant-families`
collapses only `likely_redundant` relationships and keeps an audit table.
`possible_higher_order_or_partial` relationships are retained for user review
and should not be described as definitively redundant or definitively
higher-order without assembly, read-support, known-repeat annotation or
experimental evidence.
`family_hierarchy.tsv` makes those pairwise alternatives machine-readable but
does not select a rooted hierarchy or validate a HOR model.

Known-repeat annotation and optional collapse improve interpretability, but they
do not replace assembly localization, FISH validation, or literature-supported
repeat biology interpretation.

An optional guided mode using user-supplied known-repeat FASTA files may be added later, but it is not the default MVP workflow.

The synthetic benchmark harness can measure tiny and manual scale runs, but it does not make the MVP suitable for real 7-20 Gb production analysis. Real data should be limited to pilot subsets until chunking, resumable execution, memory reporting and external benchmarks are implemented.

Discover provides incremental candidates, bounded read batches and progress logs.
Rust read-local scanning can use multiple threads; the Python scan remains GIL
limited. `--chunk-size` does not provide intra-step checkpoint/resume. Candidate
clustering still retains state that can grow with the number of candidates.
Neither the Rust backend nor parser streaming establishes full-workflow scaling.

`tandemx run --resume` operates between commands only. It validates outputs and
SHA-256 fingerprints of input files and effective commands before skipping a
step. It does not resume a partially processed FASTQ. Discover accepts multiple
Rust scan threads; this does not make every downstream step parallel.

## Development challenge findings (2026-09-06)

An opt-in `--discovery-method elastic` now supports indel-aware local alignment,
observed-unit consensus and multiple arrays per read. Legacy remains the default
and the baseline below remains archived. In the corrected elastic development
run (seed 1101), all 13 positive scenarios recovered all planted arrays and all
three families; three 100-read negative scenarios yielded no calls. This follows
an initial elastic attempt with three AT-rich false-positive reads and one
incorrect consensus length; those failed development results are retained too.
The corrected metrics are development performance, not held-out evidence.

Separate validation seed 2101 (48/48 successful runs, three repetitions per
scenario) retained perfect array recall/precision and zero negative-control
calls, but recovered only two of three sequence families in `related_families`.
The existing sketch-based clustering merged two distinct related monomers.
This failure is retained and requires method improvement; the final held-out
seeds 3101/3102/3103 have not been used.

The composition filter is heuristic, boundaries can include chance-matching
flank bases, consensus samples at most 32 observed units, and local traces are
capped at 32 million cells. Natural family divergence, long arrays, close array
transitions and genome-wide false discovery require broader validation. Neither
mode currently supplies calibrated detection probabilities. Retaining legacy
provides a reproducible algorithm ablation, not a recommended workaround for
indel-rich data.

The independent planted-array benchmark at development seed 1101 found full
array recall on clean and substitution-only scenarios but serious loss of
interval recall under indels (0.70 at 0.1% total indels, 0.0286 at 1%, and 0 at
4%, requiring IoU >=0.5 and the declared period tolerance). The one-candidate
model recovered half of the arrays in two-array reads. These are **array**
metrics: sequence-supported family recovery can remain high when array
boundaries are incomplete. The baseline is retained for algorithm development.

The old discovery command also exited with an error on valid negative controls.
This behavior has been repaired: completed zero results now have a validated
receipt and pipeline dependencies are explicitly skipped. This repair does not
resolve the indel or multiple-array algorithm limitations.

Quantify's in-process counters retain only catalogue-derived diagnostic k-mers, which avoids all-read distinct-k-mer memory growth. This is appropriate for toy and bounded pilot runs; it is not a substitute for KMC, meryl, Jellyfish, or production-scale coverage calibration.

## Claims Not Supported

Do not claim that TandemX can:

1. fully resolve megabase-scale satellite arrays;
2. precisely locate every repeat copy;
3. completely assemble satellite arrays from reads;
4. guarantee FISH probe success;
5. outperform TRF, TideHunter, TRASH or RepeatExplorer2/TAREAN on real genomes.

## Compare Limitations

`tandemx compare` is implemented only as a toy-scale single-run comparison between read-based repeat abundance and assembly-based array abundance. It does not compare multiple samples or populations, and it does not prove array collapse. `repeat_density.bedgraph` is not a family-level compare input because it has no `family_id`.
