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

This repository currently contains the project skeleton, CLI skeleton, a toy dataset simulator, and a toy-scale `discover` MVP.

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

The current `environment.yml` intentionally contains only lightweight dependencies needed for the skeleton and tests. Larger scientific dependencies should be added only when the corresponding module actually uses them.

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

`tandemx simulate toy` generates a reproducible simulated toy dataset. `tandemx discover` currently implements a toy-scale FASTA-only MVP that writes `candidate_reads.tsv`, `monomers.fa`, and `families.tsv`.

The `quantify`, `locate`, `probe`, `compare`, and `visualize` commands are still placeholders. Placeholder commands currently:

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

The tests currently validate CLI help, missing-input errors, simulator reproducibility, and toy-scale discover output behavior.
