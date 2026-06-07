# TandemX MVP Status

## Implemented

Current implemented commands:

1. `tandemx simulate toy`
2. `tandemx discover`
3. `tandemx quantify`
4. `tandemx locate`
5. `tandemx probe`
6. `tandemx visualize`

`tandemx compare` remains a placeholder. Assembly-vs-read comparison is currently produced inside `tandemx locate` as `assembly_vs_read_cn.tsv`.

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

## Test Status

The MVP includes unit and integration tests for toy simulation, discovery, quantification, localization, probe ranking, visualization, and the end-to-end toy workflow.

The test suite includes an anti-hardcoding integration test using non-default simulated repeat lengths, currently 421 bp and 729 bp. The test runs simulate, discover, quantify, locate, and probe without passing truth files to the analysis commands. Truth files are used only for test assertions.

All tests should be run in the `tandemx-dev` conda environment.

## Scope Status

The MVP is strictly toy-scale. It is intended to establish architecture, file formats, command behavior, and reproducible tests.

The MVP is not validated for real 7-20 Gb plant genomes.
