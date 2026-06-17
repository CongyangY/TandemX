# TandemX MVP Status

## Implemented

Current implemented commands:

1. `tandemx simulate toy`
2. `tandemx discover`
3. `tandemx quantify`
4. `tandemx locate`
5. `tandemx probe`
6. `tandemx visualize`
7. `tandemx validate`

`tandemx compare` remains deferred. Assembly-vs-read comparison is currently produced inside `tandemx locate` as `assembly_vs_read_cn.tsv`.

Implemented sequence input support:

1. `.fa`
2. `.fasta`
3. `.fq`
4. `.fastq`
5. `.fa.gz`
6. `.fasta.gz`
7. `.fq.gz`
8. `.fastq.gz`

## Implemented Outputs

Toy simulator:

1. `reads.fa`
2. `assembly.fa`
3. `truth_monomers.fa`
4. `truth_arrays.bed`
5. `truth_copy_number.tsv`
6. `simulation_config.yaml`

Discover:

1. `candidate_reads.tsv`
2. `monomers.fa`
3. `families.tsv`

Quantify:

1. `copy_number.tsv`

Locate:

1. `repeat_density.bedgraph`
2. `arrays.bed`
3. `assembly_vs_read_cn.tsv`

Probe:

1. `probes.fa`
2. `probes.rank.tsv`
3. `in_silico_fish.tsv`

Visualize:

1. `catalogue_summary.svg/pdf`
2. `assembly_vs_read.svg/pdf`
3. `in_silico_fish.svg/pdf`

Validate:

1. checks recognized core TSV files for required fields and numeric values;
2. checks `arrays.bed` and `repeat_density.bedgraph` coordinates;
3. checks TandemX FASTA headers for `monomers.fa` and `probes.fa`;
4. reports clear validation errors for empty or malformed recognized outputs.

## Test Status

The MVP includes unit and integration tests for toy simulation, sequence I/O, discovery, quantification, localization, probe ranking, visualization, output validation, negative CLI behavior, and the end-to-end toy workflow.

The test suite includes an anti-hardcoding integration test using non-default simulated repeat lengths, currently 421 bp and 729 bp. The test runs simulate, discover, quantify, locate, and probe without passing truth files to the analysis commands. Truth files are used only for test assertions.

Randomized toy workflow tests use fixed seeds `1`, `7`, `13`, `42`, and `99`. Each seed generates two monomer lengths in the 180-900 bp range plus controlled read and assembly copy counts, then runs simulate, discover, quantify, locate, probe, and validate. These tests are intended to catch brittle assumptions in the toy MVP, not to validate real-genome performance.

All tests should be run in the `tandemx-dev` conda environment.

## Scope Status

The MVP is strictly toy-scale. It is intended to establish architecture, file formats, command behavior, and reproducible tests.

The MVP is not validated for real 7-20 Gb plant genomes. Before real large-genome use, TandemX needs streaming algorithm optimization, parallelization, external benchmarks, and validation on real reads and assemblies.
