# TandemX Toy Workflow

This example runs the current toy-scale TandemX MVP from simulation through probe ranking.

The workflow uses only simulated sequences. It does not analyze real plant genome data and should not be used to claim performance on 7-20 Gb genomes.

## Run

From the repository root:

```bash
bash examples/toy/run_toy_workflow.sh
```

By default, outputs are written to:

```text
examples/toy/results/
```

To choose another output directory:

```bash
bash examples/toy/run_toy_workflow.sh /tmp/tandemx-toy-results
```

## Steps

The workflow runs:

```bash
tandemx simulate toy
tandemx discover
tandemx quantify
tandemx locate
tandemx compare
tandemx probe
tandemx visualize
tandemx validate
```

Expected key outputs:

```text
results/simulated/reads.fa
results/simulated/assembly.fa
results/discover/families.tsv
results/discover/monomers.fa
results/quantify/copy_number.tsv
results/locate/repeat_density.bedgraph
results/locate/arrays.bed
results/locate/assembly_vs_read_cn.tsv
results/compare/assembly_vs_read_cn.tsv
results/probe/probes.fa
results/probe/probes.rank.tsv
results/probe/in_silico_fish.tsv
results/visualize/catalogue_summary.svg
results/visualize/assembly_vs_read.svg
results/visualize/in_silico_fish.svg
```

The final validation step checks recognized TandemX output schemas, numeric fields, BED-style coordinates and TandemX FASTA headers.

Do not commit generated `results/` directories.

## Two-sample cohort example

The cohort example generates two samples with the same seeded monomers and
different planted copy counts, runs the required single-sample commands, then
builds and validates a pan-repeat catalogue:

```bash
bash examples/toy/run_toy_cohort.sh /tmp/tandemx-toy-cohort
```

Key cohort outputs include `cohort_input_qc.tsv`, `pan_families.tsv`,
`family_membership.tsv`, three abundance matrices (point, lower endpoint and
upper endpoint), the assembly representation matrix, and the four-panel
`cohort_overview.svg`/`cohort_overview.pdf`. `cohort_plot_source.tsv` is the
machine-readable source for the displayed top families and
`cohort_figure_receipt.json` records input/output hashes and editable SVG-node
counts. This example verifies the data contract and rendering path; it is not
biological cohort validation.

## TideHunter import example

With a TideHunter executable available, this workflow generates new toy reads,
runs TideHunter `-f 2`, imports the actual native calls, validates the catalogue
and uses it for TandemX quantification:

```bash
bash examples/toy/run_tidehunter_import.sh \
  /tmp/tandemx-tidehunter-import \
  /path/to/TideHunter
```

The importer validates TideHunter IDs, lengths and coordinates against the
generated reads. Its audit preserves full native IDs separately from TandemX's
normalized IDs. This is an executable interoperability example, not a detector
accuracy benchmark.
