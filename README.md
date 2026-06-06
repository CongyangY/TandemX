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

This repository currently contains only the project skeleton and CLI skeleton.

No real tandem repeat discovery, copy-number estimation, assembly localization, probe scoring, comparison, or visualization algorithm is implemented yet.

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

Available placeholder commands:

```bash
tandemx discover --help
tandemx quantify --help
tandemx locate --help
tandemx probe --help
tandemx compare --help
tandemx visualize --help
```

Each command currently:

1. parses arguments;
2. checks that required input files exist;
3. creates the output directory;
4. writes `run_config.yaml`;
5. writes `run.log`; and
6. reports that the algorithm is not implemented yet.

The toy dataset has not been added yet. It will be added in the next development phase.

Temporary skeleton run with a minimal local FASTA:

```bash
mkdir -p tmp
printf ">read1\nACGTACGTACGT\n" > tmp/reads.fa
tandemx discover \
  --reads tmp/reads.fa \
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

The tests currently validate CLI help, missing-input errors, and placeholder output-directory/config behavior.
