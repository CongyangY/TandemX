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
raw reads -> tandemx discover -> de novo repeat catalog -> quantify/locate/compare/probe/visualize
```

`tandemx discover` starts from reads only. It does not require pre-existing repeat sequences, a user library or simulator truth file. Its `monomers.fa`, `families.tsv` and `candidate_reads.tsv` outputs are the de novo repeat catalog used by downstream commands.

## Current Status

This repository currently contains a toy dataset simulator, toy-scale `discover`, `quantify`, `locate`, `compare`, `probe`, and `visualize` MVPs, and a step-level `tandemx run` orchestrator.

No production-scale tandem repeat discovery, copy-number estimation, assembly localization, assembly/read comparison, probe scoring, or visualization algorithm is available yet.

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
tandemx run --help
tandemx discover --help
tandemx quantify --help
tandemx locate --help
tandemx probe --help
tandemx compare --help
tandemx visualize --help
tandemx validate --help
tandemx annotate-repeats --help
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

`tandemx simulate toy` generates a reproducible simulated toy dataset. `tandemx discover` implements toy-scale de novo repeat discovery from reads and writes `candidate_reads.tsv`, `monomers.fa`, `families.tsv`, and `family_similarity.tsv` for pairwise catalog redundancy review. `tandemx quantify` estimates toy read-based copy number from the discovered repeat catalog and writes `copy_number.tsv`. `tandemx locate` scans a toy assembly with discovered monomer k-mers and writes `repeat_density.bedgraph`, `arrays.bed`, and a backward-compatible `assembly_vs_read_cn.tsv`. `tandemx compare` compares read-based abundance from `copy_number.tsv` with family-level assembly array abundance from `arrays.bed` and writes `assembly_vs_read_cn.tsv`. `tandemx probe` ranks toy FISH probe candidates from the discovered catalog and writes `probes.fa`, `probes.rank.tsv`, and `in_silico_fish.tsv`. `tandemx visualize` writes basic SVG/PDF static plots. `tandemx annotate-repeats` performs post hoc known-repeat annotation after discovery. `tandemx validate` checks recognized MVP outputs under a project directory.

Sequence input support is centralized in `tandemx.io.sequences`. Analysis commands can read `.fa`, `.fasta`, `.fq`, `.fastq`, and gzip-compressed `.fa.gz`, `.fasta.gz`, `.fq.gz`, and `.fastq.gz` inputs where that file type is appropriate. Readers stream records with a shared `SequenceRecord` structure and validate empty files, malformed FASTQ records, duplicate IDs, and sequence/quality length mismatches.

The `compare` MVP is an assembly/read abundance comparison for one run, not a multi-sample population comparison. It uses `copy_number.tsv` and `arrays.bed`; `repeat_density.bedgraph` is not the primary compare input because it does not contain `family_id`.

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
tandemx compare \
  --copy-number results/quantify/copy_number.tsv \
  --arrays results/locate/arrays.bed \
  --outdir results/compare
tandemx probe \
  --catalog results/discover/monomers.fa \
  --assembly results/toy/assembly.fa \
  --copy-number results/quantify/copy_number.tsv \
  --arrays results/locate/arrays.bed \
  --outdir results/probe
tandemx visualize \
  --catalog results/discover/monomers.fa \
  --copy-number results/quantify/copy_number.tsv \
  --comparison results/compare/assembly_vs_read_cn.tsv \
  --probes results/probe/probes.rank.tsv \
  --fish results/probe/in_silico_fish.tsv \
  --outdir results/visualize
tandemx validate --project results
```

In this workflow, `results/discover/monomers.fa` is a de novo discovery output. Passing it to `--catalog` in downstream commands reuses TandemX's discovered catalog; it does not mean TandemX needs repeat sequences before discovery.

The same dependency chain can be run in one command:

```bash
tandemx run \
  --reads results/toy/reads.fa \
  --assembly results/toy/assembly.fa \
  --genome-size 7744 \
  --outdir results/run \
  --steps discover,quantify,locate,compare,probe,visualize,validate \
  --kmer-backend rust
```

Without `--assembly`, locate, compare, probe, and assembly-dependent visualization steps are recorded as skipped. `--resume` is basic output-level resume, not an intra-step checkpoint. `--force` reruns selected steps. Pipeline runs write per-step logs plus `pipeline_summary.tsv`, `pipeline_summary.json`, `output_manifest.tsv`, and `run_report.md`. Read limits are passed to both discover and quantify so copy-number depth uses the same input prefix.

## Where are my outputs?

`tandemx run --outdir results/run1` uses this standard structure:

```text
results/run1/
├── discover/
├── quantify/
├── locate/
├── compare/
├── probe/
├── visualize/
├── validate/
├── logs/
├── profiles/                 # only with --profile
├── output_manifest.tsv       # file inventory, sizes, status and dependencies
├── run_report.md             # human-readable results and next commands
├── pipeline_summary.tsv      # per-step timing and status
├── pipeline_summary.json
└── pipeline.log
```

Start with `run_report.md` for a concise run overview and use `output_manifest.tsv` to locate individual files or diagnose skipped/missing outputs.

## Comparing two run directories

When two outputs were generated from the same reads but with different discovery
parameters, family counts can differ without indicating a biological conflict.
Use the run-comparison helper to make that explicit:

```bash
python benchmarks/scripts/compare_tandemx_runs.py \
  --run-a test_data/output/real_hifi_100k_discover \
  --run-b test_data/output/real_hifi_100k_reads_only \
  --outdir test_data/output/compare_real_hifi_100k
```

The script writes `compare_runs.tsv` and `compare_runs.md`. It compares the
discover `run_config.yaml`, optional `pipeline_summary.tsv`, `families.tsv`,
`candidate_reads.tsv`, and `monomers.fa`. The report states whether reads and
read limits match, lists result-affecting discover parameter differences such as
`min_support_reads`, reports candidate/family/monomer-length differences, and
marks whether the two catalogs are directly comparable.

Known repeat sequences can be compared only after de novo discovery:

```bash
tandemx annotate-repeats \
  --catalog results/run1/discover/monomers.fa \
  --known known_repeats.fa \
  --out results/run1/repeat_annotation.tsv
```

This post hoc check does not pass known repeats to `tandemx discover` and does not make them templates for candidate detection. It reports the best known-repeat match for each discovered family with Dice, Jaccard, containment and local-identity metrics.

`tandemx discover --collapse-redundant-families` is optional and off by default. When enabled, it writes `collapsed_families.tsv`, `collapsed_monomers.fa`, and `family_collapse.tsv`, but only collapses pairs classified as `likely_redundant`. Pairs labelled `possible_higher_order_or_partial` are retained and should be reviewed; TandemX does not claim they are definitely redundant or definitely higher-order repeats.

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
docs/sensitivity_validation.md
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

Discover uses a repeated-k-mer spacing prefilter and bounded local period refinement; it no longer scans every possible period against every base. For real HiFi subset pilots, limit work explicitly and monitor the live terminal progress plus `run.log`:

```bash
tandemx discover \
  --reads subset_lane1.fastq.gz subset_lane2.fastq.gz \
  --outdir pilot_discover \
  --max-reads 10000 \
  --max-read-bases 200000000 \
  --min-read-length 1000 \
  --min-period 2 \
  --max-period 2000 \
  --kmer-backend rust \
  --threads 8 \
  --count-threads 4 \
  --progress-every 1000
```

`--reads` accepts one or more FASTA/FASTQ files, including gzip-compressed files. Multiple files are streamed in the order supplied and analyzed as one merged read set. Duplicate read IDs across input files are treated as an input error, which helps catch accidental repeated file arguments. `candidate_reads.tsv` and `run.log` are created at startup and flushed during processing. The terminal progress line refreshes in place and reports the current step, processed reads and bases, elapsed time, estimated total runtime, remaining time, reads/min and MB/min. Discover pre-counts input reads and bases before scanning, using up to `--count-threads 4`, so total and remaining time are available even without `--max-reads`. Use `--no-progress` for non-interactive batch logs. `--kmer-backend python` remains the default and fallback; `--kmer-backend rust` accelerates the same read-local algorithm when the extension is installed. With the Rust backend, `--threads` parallelizes read-local scanning; the default request is 8 threads, capped at the smaller of 64 and half of available logical CPUs. See `docs/performance.md` for parity results and scaling limits.

The default minimum period is 2 bp so short tandem repeats such as di-, tri- and heptanucleotide repeats can be reported when they span enough read sequence. Short or low-complexity candidates are flagged with warnings. Set `--min-period 20` when a run should focus only on longer satellite-like monomers.

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

For step-level timing, use:

```bash
python benchmarks/scripts/run_pipeline_benchmark.py \
  --reads subset.fastq.gz \
  --genome-size 16000000000 \
  --outdir /tmp/tandemx_pipeline_pilot \
  --steps discover,quantify,validate \
  --kmer-backend rust \
  --max-reads 100000 \
  --profile
```
