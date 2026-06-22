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
