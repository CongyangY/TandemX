# TandemX

TandemX is an open-source computational genomics project for tandem repeat and satellite repeat analysis in large plant genomes.

The intended long-term design is a read-first, assembly-aware workflow for:

1. discovering candidate tandem repeat monomers from HiFi-like reads;
2. estimating repeat family copy number with diagnostic k-mer depth;
3. localizing repeat evidence on assemblies when assemblies are available;
4. comparing read-based and assembly-based abundance to flag possible under-representation;
5. prioritizing FISH probe candidates; and
6. generating documented, publication-oriented visual outputs.

## Current Status

This repository currently contains the project skeleton, CLI skeleton, a toy dataset simulator, toy-scale `discover`, `quantify`, `locate`, `probe`, and `visualize` MVPs.

No production-scale tandem repeat discovery, copy-number estimation, assembly localization, probe scoring, comparison, or visualization algorithm is implemented yet.

The first implementation target is a toy-scale MVP. It should run on small simulated data and should not claim support for real 7-20 Gb plant genomes until benchmarked.

## Development Environment

Use the project-local conda environment definition. Do not install TandemX dependencies into `base`.

```bash
conda env create -f environment.yml
conda activate tandemx-dev
pip install -e .
pytest
```

The current `environment.yml` intentionally contains only lightweight dependencies needed for the MVP and tests. `matplotlib` is included for static SVG/PDF visualization. Larger scientific dependencies should be added only when the corresponding module actually uses them.

## Install

After activating `tandemx-dev`, install the package in editable mode:

```bash
pip install -e .
```

## CLI

Available commands:

```bash
tandemx simulate toy --help
tandemx discover --help
tandemx quantify --help
tandemx locate --help
tandemx probe --help
tandemx compare --help
tandemx visualize --help
```

## Quick Start: Toy Workflow

After activating `tandemx-dev` and installing TandemX in editable mode, run:

```bash
bash examples/toy/run_toy_workflow.sh
```

This writes simulated input and MVP outputs under:

```text
examples/toy/results/
```

Generated results are ignored by git and should not be committed.

`tandemx simulate toy` generates a reproducible simulated toy dataset. `tandemx discover` currently implements a toy-scale FASTA-only MVP that writes `candidate_reads.tsv`, `monomers.fa`, and `families.tsv`. `tandemx quantify` estimates toy read-based copy number from diagnostic k-mer depth and writes `copy_number.tsv`. `tandemx locate` scans a toy assembly with monomer k-mers and writes `repeat_density.bedgraph`, `arrays.bed`, and `assembly_vs_read_cn.tsv`. `tandemx probe` ranks toy FISH probe candidates and writes `probes.fa`, `probes.rank.tsv`, and `in_silico_fish.tsv`. `tandemx visualize` writes basic SVG/PDF static plots.

The `compare` command is still a placeholder. Placeholder commands currently:

1. parse arguments;
2. check that required input files exist;
3. create the output directory;
4. write `run_config.yaml`;
5. write `run.log`; and
6. report that the algorithm is not implemented yet.

Generate a toy dataset and run discover:

```bash
mkdir -p results
tandemx simulate toy --outdir results/toy
tandemx discover \
  --reads results/toy/reads.fa \
  --outdir results/discover
tandemx quantify \
  --reads results/toy/reads.fa \
  --catalog results/discover/monomers.fa \
  --genome-size 7744 \
  --outdir results/quantify
tandemx locate \
  --assembly results/toy/assembly.fa \
  --catalog results/discover/monomers.fa \
  --copy-number results/quantify/copy_number.tsv \
  --window-size 500 \
  --step-size 250 \
  --outdir results/locate
tandemx probe \
  --catalog results/discover/monomers.fa \
  --assembly results/toy/assembly.fa \
  --copy-number results/quantify/copy_number.tsv \
  --arrays results/locate/arrays.bed \
  --outdir results/probe
tandemx visualize \
  --catalog results/discover/monomers.fa \
  --copy-number results/quantify/copy_number.tsv \
  --comparison results/locate/assembly_vs_read_cn.tsv \
  --probes results/probe/probes.rank.tsv \
  --fish results/probe/in_silico_fish.tsv \
  --outdir results/visualize
```

## Tests

```bash
pytest
```

Formal development and test runs should use Python 3.10 or newer.

For this repository, use the `tandemx-dev` conda environment:

```bash
conda activate tandemx-dev
pytest
```

The tests currently validate CLI help, missing-input errors, simulator reproducibility, toy-scale discover output behavior, toy-scale quantify behavior, toy-scale assembly localization behavior, toy-scale probe ranking behavior, and basic static visualization output.

The test suite also includes non-default toy repeat lengths, currently 421 bp and 729 bp, to check that the MVP workflow is not hardcoded to the simulator defaults. This is an anti-hardcoding check for toy-scale behavior, not validation for real 7-20 Gb plant genomes.

## MVP Documentation

Current status and limitations:

```text
docs/mvp_status.md
docs/known_limitations.md
docs/roadmap.md
```
