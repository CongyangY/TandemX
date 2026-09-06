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

The expanded research/release programme is tracked in [docs/release_program.md](docs/release_program.md).
It includes an independent challenge benchmark and public-data provenance, with
development and held-out seeds separated. Current development benchmarks expose
important indel-boundary and multiple-array limitations; see
[benchmarks/README.md](benchmarks/README.md). An opt-in indel-aware multiple-array
method is now available for validation (`--discovery-method elastic`); the
original method remains available as `legacy` and is still the default.
Elastic mode now uses explicit circular edit-similarity monomer clusters
(default 95%), with original candidate sequences and assignment evidence retained.
These clusters do not define biological family ancestry.
The broader comparator and multi-metric contract is in
[docs/comparator_matrix.md](docs/comparator_matrix.md). These results do not establish
production readiness or a universal advantage over external tools.

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
The package declares Matplotlib for its existing SVG/PDF commands; benchmark
scripts additionally use PyYAML (`pip install -e '.[benchmark]'`). Both are already
included in `environment.yml`. The GitHub workflow tests source installation in
`tandemx-dev` on Linux and macOS; a configured workflow is not a completed CI run.

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

To exercise the experimental elastic discovery method after creating toy data:

```bash
tandemx simulate toy --outdir /tmp/tandemx-toy
tandemx discover --reads /tmp/tandemx-toy/reads.fa --outdir /tmp/tandemx-elastic \
  --discovery-method elastic --kmer-backend rust --threads 1
tandemx validate --project /tmp/tandemx-elastic
```

See [examples/toy/elastic.md](examples/toy/elastic.md) for its evidence and limits.
`--clustering-method legacy` preserves the previous clustering for ablation;
`--cluster-identity` controls sequence resolution.

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

Without `--assembly`, locate, compare, probe, and assembly-dependent visualization steps are recorded as skipped. `--resume` skips a completed step only when its outputs validate and a SHA-256 fingerprint of its inputs and effective command still matches; it is not an intra-step checkpoint. `--force` reruns selected steps. Pipeline runs write per-step logs plus `pipeline_summary.tsv`, `pipeline_summary.json`, `output_manifest.tsv`, and `run_report.md`. Read limits are passed to both discover and quantify so copy-number depth uses the same input prefix.

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

Valid reads with no families passing the configured filters produce a successful
`discover` run with an empty catalog, header-only tables, and a checksum-backed
`discovery_summary.json`. `tandemx run` records dependent commands as
`skipped_no_discovered_families`; `validate` accepts zero results only when this
completion receipt and output checksums agree. Malformed/empty sequence input
still fails. An arbitrary empty catalog is not accepted by downstream commands.

```bash
pytest
```

Formal development and test runs should use Python 3.10 or newer.

For this repository, use the `tandemx-dev` conda environment:

```bash
conda activate tandemx-dev
pytest
```

The tests validate CLI help, missing-input errors, simulator and backend reproducibility, streaming FASTA/FASTQ/gzip input behavior, local repeat boundaries, sequence-aware clustering, circular diagnostic k-mers, toy-scale quantify/localize/probe behavior, output schema validation, resume fingerprints, and static visualization output.

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

The runner executes `simulate -> discover -> quantify -> locate -> probe -> validate`, writes `benchmark_summary.tsv` and `accuracy_summary.tsv`, records per-command wall time and peak resident memory, and stores per-command logs. Truth matching uses sequence identity plus length rather than length alone. Only the `tiny` scale is intended for pytest. Larger synthetic scales are manual tests and do not imply real 7-20 Gb production readiness.

## Discover Pilot Controls

Discover uses a repeated-k-mer spacing prefilter and bounded local period refinement; it no longer scans every possible period against every base. For real HiFi subset pilots, limit work explicitly and monitor the live terminal progress plus `run.log`:

```bash
tandemx discover \
  --reads subset_lane1.fastq.gz subset_lane2.fastq.gz \
  --outdir pilot_discover \
  --genome-size 16000000000 \
  --enable-auto-discovery-budget \
  --target-discovery-coverage 10.0 \
  --min-read-length 1000 \
  --min-period 2 \
  --max-period 2000 \
  --kmer-backend rust \
  --threads 8 \
  --count-threads 4 \
  --progress-every 1000
```

`--reads` accepts one or more FASTA/FASTQ files, including gzip-compressed files. Multiple files are streamed in the order supplied and analyzed as one merged read set. Duplicate read IDs in discovery are checked exactly with a bounded in-memory tracker that spills to a temporary SQLite index. `candidate_reads.tsv` and `run.log` are created at startup and flushed during processing. The terminal progress line refreshes in place and reports the current step, processed reads and bases, elapsed time, reads/min and MB/min. TandemX avoids an otherwise wasteful full input pre-count; a synchronous count is performed only when automatic budgeting is explicitly enabled without `--genome-size`. Use `--no-progress` for non-interactive batch logs. `--kmer-backend auto` is the default and uses Rust when the compiled extension and k-mer size are supported; use `--kmer-backend python` only for fallback/debugging. With the Rust backend, `--threads` parallelizes read-local scanning; `--chunk-size` and `--chunk-bases` jointly bound each scan batch. By default discover scans the full input unless you set explicit `--max-reads/--max-read-bases` limits. If you want bounded large-input discovery, enable it explicitly with `--enable-auto-discovery-budget`; with `--genome-size`, TandemX caps discovery to approximately `--target-discovery-coverage` genome equivalents, and with multiple files the bounded mode switches to round-robin file streaming. See `docs/performance.md` for parity results and scaling limits.

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

### Independent benchmark sequence evaluation

The benchmark extra pins the small native `edlib==1.3.9.post1` dependency for an
independent global edit-distance evaluator; TandemX discovery does not call it.
It is included in the development environment and test/benchmark extras:
`pip install -e '.[test,benchmark]'`. The evaluator exhausts rotations and strands
and is tested against a separate dynamic-programming reference. Both the original
strict equal-length endpoint and the added indel-aware endpoint are reported.
Union-of-interval coverage separates duplicated/harmonic calls from incorrectly
labelled bases. See [benchmark definitions](docs/file_formats.md#challenge-measurements).

The expanded SRF comparison executes the documented KMC-to-abundance workflow,
including input counting and mapping. See [SRF workflow evidence and scope](docs/srf_workflow.md)
for pinned tools, recorded compatibility patches, independent count checks,
empty-count handling, native versus normalized outputs, and measured costs.

Elastic discovery in Rust mode now constructs the bounded seed histogram natively
as well as aligning repeat copies. The Python reference retains the same counting
and capping rules. Rebuild the extension when updating: `pip install -e .`.
Paired source-snapshot comparisons check all six discovery outputs byte-for-byte
before treating a timing improvement as an equivalent-output optimization.
The executed 96-run development comparison preserved all six outputs and reduced
median wall time in all 16 scenarios (1.20–4.45× speedup); peak RSS changed by
-13.09% to +2.03%, including increases in six scenarios. These are same-tool
observations, not superiority to external tools. Raw pairs and the four-panel
figure are in `paper/evidence/native_seed_paired`.

Discovery also uses the selected backend for the exhaustive representative-pair
audit. Native ungapped comparisons preserve the Python audit's values and tie
rules; each representative's k-mer set is built once and pair rows are streamed
to `family_similarity.tsv`. The audit is still quadratic in catalogue size.
It is a redundancy heuristic, separate from gapped circular sequence clustering.
See [algorithms](docs/algorithms.md#representative-pair-redundancy-audit).
For larger catalogues, `discover` and `run` accept `--family-audit related`.
An exact k-mer intersection bound skips pairs that cannot satisfy the existing
rules and emits every non-distinct pair. Catalogue, warning and collapse semantics
are preserved; omitted distinct pairs and scored counts are explicit in
`family_audit_summary.json`. The default remains `full` for compatibility.

The [conditional abundance experiment](docs/abundance_benchmark.md) tests copy
number, localization and engineered assembly under-representation independently
of de novo discovery. It measures empirical interval coverage without treating
the existing diagnostic-k-mer spread as a calibrated confidence interval.
Its first full baseline completed 177 commands but covered truth in only 6/81
reported intervals; this is an identified limitation, not calibrated inference.
The [real and simulated cohort/QC programme](docs/cohort_and_qc.md) specifies
species/material breadth, scale ladders, full-file validation and independent
evidence requirements. Current development/toy results do not satisfy those gates.

[Complete FASTQ QC and nested hash sampling](docs/complete_data_qc.md) verify
whole inputs and retain reproducible read IDs, achieved sizes and distributions.
The [experimental read-cluster model](docs/read_cluster_quantification.md)
replaces within-k-mer spread by a read-level sampling calculation in development
replays; it remains separate from the default CLI pending biological calibration.

Complete Mo17 CCS input QC now covers 407,670 reads / 5.625 Gb, with nested
11.68 Mb, 111.51 Mb and 1.129 Gb samples and a checked 10-chromosome reference.
Source receipts and a four-panel QC figure are in `paper/evidence/Mo17_input_qc`.
The [real-input comparator pilot](docs/real_comparator_pilot.md) uses identical
sample sequences and preserves native outputs; call counts are not accuracy.
The complete Col-0N FASTQ also passed QC (933,904 reads, 14.647 Gb); full-library
random sampling is under way. [Published Mo17 regions](docs/published_mo17_regions.md)
retain the original spreadsheet, source cells, assembly-coordinate checks and
the distinction between a mixed repeat region and base-level satellite truth.
