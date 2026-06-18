# TandemX

TandemX is an open-source computational genomics project for tandem repeat and satellite repeat analysis in large plant genomes.

The intended long-term design is a read-first, assembly-aware workflow for:

1. discovering candidate tandem repeat monomers from HiFi-like reads;
2. estimating repeat family copy number with diagnostic k-mer depth;
3. localizing repeat evidence on assemblies when assemblies are available;
4. comparing read-based and assembly-based abundance to flag possible under-representation;
5. prioritizing FISH probe candidates; and
6. generating documented, publication-oriented visual outputs.

The default TandemX workflow is de novo:

```text
raw reads -> tandemx discover -> de novo repeat catalog -> quantify/locate/probe/visualize
```

`tandemx discover` starts from reads only. It does not require pre-existing repeat sequences, a user library or simulator truth file. Its `monomers.fa`, `families.tsv` and `candidate_reads.tsv` outputs are the de novo repeat catalog used by downstream commands.

## Current Status

This repository currently contains a toy dataset simulator and toy-scale `discover`, `quantify`, `locate`, `probe`, and `visualize` MVPs.

No production-scale tandem repeat discovery, copy-number estimation, assembly localization, probe scoring, comparison, or visualization algorithm is available yet.

The first implementation target is a toy-scale MVP. It should run on small simulated data and should not claim support for real 7-20 Gb plant genomes until benchmarked.

## Development Environment

Use the project-local conda environment definition. Do not install TandemX dependencies into `base`.

```bash
conda env create -f environment.yml
conda activate tandemx-dev
pip install -e .
pytest
```

`environment.yml` includes Rust and maturin to build the optional compiled read-local discovery backend. No production global k-mer counter is bundled.

## Install

After activating `tandemx-dev`, install the package in editable mode:

```bash
pip install -e .
```

The editable install uses maturin and builds the Rust extension in addition to the Python package. During backend development, the equivalent explicit command is:

```bash
maturin develop --release
```

Use `--release` for performance measurements. The Python backend remains available if a compiled extension is not installed.

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
tandemx validate --help
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

`tandemx simulate toy` generates a reproducible simulated toy dataset. `tandemx discover` implements toy-scale de novo repeat discovery from reads and writes `candidate_reads.tsv`, `monomers.fa`, and `families.tsv`. `tandemx quantify` estimates toy read-based copy number from the discovered repeat catalog and writes `copy_number.tsv`. `tandemx locate` scans a toy assembly with discovered monomer k-mers and writes `repeat_density.bedgraph`, `arrays.bed`, and `assembly_vs_read_cn.tsv`. `tandemx probe` ranks toy FISH probe candidates from the discovered catalog and writes `probes.fa`, `probes.rank.tsv`, and `in_silico_fish.tsv`. `tandemx visualize` writes basic SVG/PDF static plots. `tandemx validate` checks recognized MVP outputs under a project directory.

Sequence input support is centralized in `tandemx.io.sequences`. Analysis commands can read `.fa`, `.fasta`, `.fq`, `.fastq`, and gzip-compressed `.fa.gz`, `.fasta.gz`, `.fq.gz`, and `.fastq.gz` inputs where that file type is appropriate. Readers stream records with a shared `SequenceRecord` structure and validate empty files, malformed FASTQ records, duplicate IDs, and sequence/quality length mismatches.

The `compare` command is deferred in the current MVP. Deferred commands currently:

1. parse arguments;
2. check that required input files exist;
3. create the output directory;
4. write `run_config.yaml`;
5. write `run.log`; and
6. report that the algorithm is deferred in this MVP.

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
tandemx validate --project results
```

In this workflow, `results/discover/monomers.fa` is a de novo discovery output. Passing it to `--catalog` in downstream commands reuses TandemX's discovered catalog; it does not mean TandemX needs repeat sequences before discovery.

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

The tests currently validate CLI help, missing-input errors, simulator reproducibility, streaming FASTA/FASTQ/gzip input behavior, toy-scale discover output behavior, toy-scale quantify behavior, toy-scale assembly localization behavior, toy-scale probe ranking behavior, output schema validation, and basic static visualization output.

The test suite also includes non-default toy repeat lengths, currently 421 bp and 729 bp, to check that the MVP workflow is not tied to the simulator defaults. Randomized toy workflow tests run fixed seeds `1`, `7`, `13`, `42`, and `99`; each seed generates two monomer lengths and controlled copy counts, then runs simulate, discover, quantify, locate, probe, and validate. Simulator truth files are used only for simulator checks and test assertions, not as command inputs in the analysis workflow.

These tests improve engineering reliability for the toy MVP. They are not validation for real 7-20 Gb plant genomes. Real large-genome support still requires streaming optimization, parallel execution, external benchmarking, and validation on real reads and assemblies.

## MVP Documentation

Current status and limitations:

```text
docs/mvp_status.md
docs/known_limitations.md
docs/roadmap.md
docs/benchmark_plan.md
docs/real_data_pilot_plan.md
```

## Synthetic Benchmarks

TandemX includes a synthetic benchmark harness for engineering checks before real-data pilots:

```bash
python benchmarks/scripts/run_synthetic_benchmark.py \
  --config benchmarks/configs/synthetic_scale.yaml \
  --scale tiny \
  --outdir /tmp/tandemx_benchmark_tiny
```

The runner executes `simulate -> discover -> quantify -> locate -> probe -> validate`, writes `benchmark_summary.tsv` and `accuracy_summary.tsv`, and stores per-command logs. Only the `tiny` scale is intended for pytest. Larger synthetic scales are manual tests and do not imply real 7-20 Gb production readiness.

## Discover Pilot Controls

Discover uses a repeated-k-mer spacing prefilter and bounded local period refinement; it no longer scans every possible period against every base. For real HiFi subset pilots, limit work explicitly and monitor `run.log`:

```bash
tandemx discover \
  --reads subset.fastq.gz \
  --outdir pilot_discover \
  --max-reads 10000 \
  --max-read-bases 200000000 \
  --min-read-length 1000 \
  --min-period 20 \
  --max-period 2000 \
  --kmer-backend rust \
  --progress-every 1000
```

`candidate_reads.tsv` and `run.log` are created at startup and flushed during processing. `--kmer-backend python` remains the default and fallback; `--kmer-backend rust` accelerates the same read-local algorithm when the extension is installed. See `docs/performance.md` for parity results and scaling limits.

Inspect and benchmark a local real-read subset without running downstream biological analyses:

```bash
python benchmarks/scripts/inspect_reads.py \
  --reads subset.fastq.gz \
  --output /tmp/read_stats.tsv
python benchmarks/scripts/run_real_read_pilot_benchmark.py \
  --reads subset.fastq.gz \
  --max-reads 1000,5000,10000,25000 \
  --outdir /tmp/tandemx_real_pilot
```

The real-read runner executes only `discover` and `validate`; it does not read simulator truth files.
