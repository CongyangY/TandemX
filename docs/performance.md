# TandemX Discover Performance

## Retired Exhaustive Algorithm

The original toy MVP tested every integer period between the configured minimum and maximum for every read, then compared nearly every base against the shifted sequence. Its per-read cost was approximately:

```text
O(read_length * number_of_tested_periods)
```

For 100,000 HiFi reads with an N50 near 16 kb and a 20-2,000 bp period range, this implies trillions of Python-level base comparisons. The implementation also wrote candidate output only after all reads completed, so a long run appeared to produce no output. That design was not acceptable for real-read subset pilots and has been removed from the discover execution path.

## Current Spacing-prefilter Algorithm

`tandemx discover` now uses two bounded stages per read:

1. stream one read from `tandemx.io.sequences`;
2. skip reads that are too short;
3. directly scan bounded short periods from 2-19 bp so high-copy STR-like repeats are not discarded by the k-mer seed filter;
4. extract canonical, non-low-complexity k-mers with a streaming 2-bit rolling encoder in the selected Python or Rust backend and retain only repeated within-read seeds;
5. cap stored positions and inspected pairs per k-mer;
6. build a bounded seed-spacing histogram;
7. select only the strongest spacing peaks;
8. refine peaks and plausible fundamental divisors against a linear-time local shifted-identity interval, requiring evidence across at least two repeat units;
9. release the read-local seed structures before reading the next record.

The algorithm no longer compares every base for every possible period. Its Python pilot cost is dominated by per-read k-mer extraction plus bounded spacing/refinement work.

## Streaming And Interruption

At command start, discover creates:

1. `run.log`;
2. `run_config.yaml` with running status;
3. `candidate_reads.tsv` with its header.

Each accepted candidate is appended and flushed immediately. Discover normally starts scanning without a separate full-input count, avoiding doubled input I/O and an uncancellable background task. A synchronous pre-count is used only when the user explicitly enables automatic budgeting without providing `--genome-size`. Progress logs report processed reads, processed bases, candidate reads, elapsed time, reads/s and MB/s; a percentage and remaining-time estimate are shown only when a configured limit or counted total is available. Ctrl-C leaves existing candidate output available for diagnosis.

## Pilot Controls

Use these controls for real-read subsets:

```text
--max-reads
--max-read-bases
--sample-rate
--seed
--progress-every
--count-threads
--min-read-length
--min-period
--max-period
--kmer-size
--top-periods
--min-seed-occurrences
--min-spacing-support
--max-pairs-per-kmer
--chunk-size
--chunk-bases
--threads
--no-progress
```

`--threads` defaults to a request of 8 for `discover`, capped at the smaller of 64 and half of available logical CPUs. `--kmer-backend auto` is the default and uses Rust when the compiled extension is available and the requested k-mer size is supported by Rust, so normal installs with the default k-mer size get the multi-threaded scan path by default. Multi-threaded scanning is enabled only with the Rust backend, whose PyO3 implementation releases the Python GIL during read-local scanning. The Python backend records the requested thread setting but scans reads serially because its CPU-heavy path is GIL-bound.

`--count-threads` is used only for the explicit automatic-budget pre-count when no genome size is available. A single gzip stream is counted by one worker unless an external tool such as `seqkit` is available.

`--reads` accepts multiple files. They are streamed in the supplied order and
merged for the run without preloading all reads into memory. Duplicate read IDs
across files are treated as input errors.

For large real-read inputs, discover now has an optional bounded mode. By default it scans the full input unless the user sets explicit `--max-reads` or `--max-read-bases`. If `--enable-auto-discovery-budget` is set and the user does not set explicit `--max-reads` or `--max-read-bases`, then:

1. with `--genome-size`, discover targets approximately `--target-discovery-coverage` genome equivalents, bounded by `--auto-discovery-max-bases`;
2. without `--genome-size`, discover enables the same bounded mode once the counted input exceeds `--auto-discovery-trigger-bases`.

When this bounded mode is active and multiple files are supplied, discover reads files in round-robin order so early stopping is less biased toward the first file.

The default `--min-period` is 2 bp. Use `--min-period 20` when engineering pilots
should focus on longer satellite-like monomers and ignore STR-like periods.

`--chunk-size` and `--chunk-bases` jointly bound how many selected reads and sequence bases are submitted to the Rust thread pool at once. They are not checkpoint boundaries. Step-level `tandemx run --resume` uses validated outputs plus SHA-256 input/command fingerprints, but intra-step checkpointing is not implemented.

`tandemx quantify` uses the same live progress style while scanning reads for diagnostic k-mers. When commands are launched through `tandemx run`, child command output is streamed to the terminal while still being saved under `logs/`.

## Backend Boundary

`--kmer-backend auto` is the default. It resolves to `rust` when the compiled extension is importable and the k-mer size is supported by Rust, and to `python` otherwise; the resolved backend is written to `run_config.yaml` and `run.log`. `--kmer-backend rust` implements single-read rolling canonical k-mers, repeated positions, bounded spacing histograms, top periods and local scoring behind a PyO3 interface. FASTA/FASTQ parsing, filters, clustering and file output remain in Python. Rust scan results are written in input order even when multiple worker threads finish out of order.

Production-scale global k-mer counting is a separate problem and should use optional mature backends such as KMC, meryl or Jellyfish instead of a new TandemX counter. The Rust backend does not replace those tools.

Quantify derives circular diagnostic k-mers from the discovered families before scanning reads. Both backends count only those targets, so memory scales with the catalogue rather than all distinct read k-mers. Rust receives bounded batches and counts them while the GIL is released. The implementation remains a target counter, not a production global k-mer counter.

Locate builds one family-specific k-mer index, streams assembly FASTA in bounded
chunks with exact cross-chunk rolling context, merges hits online, and calculates
window coverage from interval unions. Probe ranking likewise indexes all candidates
before one bounded-chunk assembly pass. These changes remove the
previous `families × assembly` and `probes × assembly` rescans. Result tables can
still grow with the number of reported arrays, windows and probe regions.

## Remaining Limits

The spacing prefilter and Rust read-scanning threads enable larger subsets but do not make TandemX ready for full 7-20 Gb genomes or 100 Gb read sets. Remaining work includes multiprocessing or distributed chunks, intra-step checkpoints, bounded on-disk result tables and validation on real HiFi subsets. The synthetic benchmark runner now records per-command peak resident memory with `wait4` on supported Unix platforms.

## Synthetic Reference Measurements

Current correctness-and-streaming reference on the development workstation on
2026-07-15 (`tiny`, 1,000 × 1.8 kb reads; 421/729 bp truth monomers):

| Step | Runtime | Peak RSS | Result |
|---|---:|---:|---|
| discover | 1.32 s | 31.0 MB | 869 candidates, 2 families |
| quantify | 0.12 s | 25.0 MB | 2 estimates |
| locate | 0.08 s | 22.8 MB | validated |
| probe | 0.09 s | 22.4 MB | validated |

Both truth monomers were recovered at exact length and sequence identity 1.0;
copy-number relative errors were 8.84% and 3.94% in this engineering dataset.
Two independent end-to-end runs produced byte-identical content for all 11
primary analysis outputs. Runtime and RSS naturally vary by host and run.

## Single-thread External Tool Comparison

The current truth-aware comparison was rerun on 2026-07-15 against TRF
4.10.0-rc.2 and TideHunter 1.5.5 using the five existing repository benchmark
FASTA files. Inputs were reused without rewriting and their SHA-256 digests
were recorded. TandemX and TideHunter were explicitly limited to one thread;
TRF ran as a single process. Each tool ran three times per dataset, and the
table reports medians. Direct-process peak RSS was measured with `wait4`.

| Dataset | Tool | Runtime | Peak RSS | Recall | Precision |
|---|---|---:|---:|---:|---:|
| 2.5 Mb, 171 bp, 0% error | TandemX | 1.48 s | 30.5 MB | 1.000 | 1.000 |
|  | TRF | 2.73 s | 23.7 MB | 1.000 | 1.000 |
|  | TideHunter | 2.27 s | 37.6 MB | 1.000 | 1.000 |
| 5 Mb, 171 bp, 5% error | TandemX | 2.94 s | 37.4 MB | 1.000 | 1.000 |
|  | TRF | 6.74 s | 23.8 MB | 1.000 | 1.000 |
|  | TideHunter | 4.48 s | 147.6 MB | 1.000 | 0.9987 |
| 10 Mb, 421 bp, 5% error | TandemX | 2.46 s | 49.1 MB | 1.000 | 1.000 |
|  | TRF | 10.32 s | 56.0 MB | 1.000 | 1.000 |
|  | TideHunter | 11.77 s | 277.6 MB | 1.000 | 0.9960 |

Across the 2.5-10 Mb workloads, TandemX was 1.84-4.19x faster than TRF and
1.52-4.78x faster than TideHunter. TandemX used less peak memory than
TideHunter on all five datasets, but TRF used less memory on four of the five;
TandemX used 12.2% less than TRF only on the 10 Mb dataset. TandemX and TRF
were tied on the current simple accuracy model, so these data support a scoped
throughput claim but not a general accuracy-superiority claim.

The inputs contain a single whole-read tandem family with substitution-only
errors and random-background negatives. They do not model indels, partial or
mixed arrays, low-complexity genomic backgrounds, higher-order repeats or a
real plant HiFi quality distribution. The smallest 0.4 Mb workload also showed
process cold-start variability, so the strongest runtime statement is limited
to the larger workloads.

A separate 10.13 MB assembly scaling check reduced peak RSS from 60.3 to
40.9 MB for locate and from 59.3 to 37.5 MB for probe after bounded FASTA
chunking, while retaining identical validated outputs. The remaining fixed
Python/import overhead dominates these small absolute sizes.

Measurements on the development workstation on 2026-06-18:

| Scale | Reads | Read bp | Discover runtime | Throughput | Families |
|---|---:|---:|---:|---:|---:|
| tiny | 1,000 | 1.8 Mb | 3.85 s | 266.5 reads/s; 0.480 MB/s | 2 |
| small | 10,000 | 18 Mb | 30.23 s | 349.5 reads/s; 0.629 MB/s | 2 |

These synthetic reads are 1.8 kb, much shorter than 16 kb HiFi reads. Real-read runtime must be estimated from processed bases, not read count alone.

## Real HiFi Pilot Profile

A local 100,000-read FASTQ was used only as an engineering workload on 2026-06-18. It contains 1,406,681,945 read bases with a 14,069 bp read N50. No biological interpretation was performed and the data are not part of the repository.

The 1,000-read runtime with `--min-period 50 --max-period 1000 --top-periods 3` improved from 22.63 seconds at commit `199d885` to 8.63 seconds after the read-local optimizations below. Candidate counts were unchanged in the comparison.

| Reads | Processed bp | Runtime | Throughput | Candidate rate | Validation |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 14,054,102 | 8.63 s | 1.639 MB/s | 0.2000% | passed |
| 5,000 | 70,323,743 | 43.17 s | 1.631 MB/s | 0.0800% | passed |
| 10,000 | 140,633,200 | 92.89 s | 1.515 MB/s | 0.0600% | passed |
| 25,000 | 351,499,399 | 230.53 s | 1.525 MB/s | 0.0640% | passed |
| 50,000 | 703,466,091 | 464.69 s | 1.514 MB/s | 0.0600% | passed |

`cProfile` on 10,000 reads attributed about 78% of cumulative runtime to read-local rolling k-mer extraction, 14% to bounded shifted-identity scoring, 4% to spacing histograms, and less than 1% to FASTQ parsing. File output and reverse-complement string construction were not material costs; canonical reverse complements are encoded incrementally as integers.

## Rust Backend Measurements

The release-mode Rust extension was measured on the same local FASTQ and parameters. Python and Rust returned identical candidate and family counts at 1,000, 5,000 and 10,000 reads. On a synthetic 1,000-read parity run, family monomer lengths were identical and candidate counts differed by 3.3%.

| Reads | Python runtime | Rust runtime | Rust throughput | Speedup | Rust candidates | Rust families |
|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 9.51 s | 0.89 s | 16.961 MB/s | 10.64x | 2 | 1 |
| 5,000 | 47.88 s | 4.40 s | 16.198 MB/s | 10.87x | 4 | 2 |
| 10,000 | 92.57 s | 8.44 s | 16.795 MB/s | 10.97x | 6 | 2 |
| 25,000 | NA | 20.94 s | 16.841 MB/s | NA | 16 | 3 |
| 50,000 | NA | 42.84 s | 16.443 MB/s | NA | 30 | 3 |

All Rust benchmark outputs passed schema validation. These are engineering measurements, not biological validation. The Rust backend remains single-process and per-read; full 7–20 Gb production use still requires checkpoint/resume, memory profiling, robust clustering, downstream scaling and external tool comparisons.

The 100,000-read standalone discover run processed 1,406,681,945 bp in 83.25 seconds (16.917 MB/s), found 50 candidate reads and four families, and passed output validation. Its progress log reached all 100,000 reads. This is an engineering pilot, not evidence that the recovered families are biologically correct.

## Step-level Pipeline Measurements

The reads-only 100,000-read run used cProfile for each step and the Rust backend:

| Step | Runtime | Output validation |
|---|---:|---:|
| discover | 101.08 s | passed |
| quantify | 32.68 s | passed |
| validate | 0.09 s | passed |

Profiling overhead explains why discover is slower here than the standalone measurement. Total measured step time was 133.85 seconds. Discover remains the main real-read pilot bottleneck; quantify is material but bounded-memory after target-only counting.

The profiled toy full workflow recorded discover 0.10 s, quantify 0.17 s, locate 0.12 s, probe 0.15 s, visualize 9.04 s, and validate 0.08 s. Matplotlib and process startup dominate this tiny visualization measurement, so it does not predict chromosome-scale behavior.
