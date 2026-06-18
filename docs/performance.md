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
2. skip reads that are too short or low complexity;
3. extract canonical, non-low-complexity k-mers with a streaming 2-bit rolling encoder and retain only repeated within-read seeds;
4. cap stored positions and inspected pairs per k-mer;
5. build a bounded seed-spacing histogram;
6. select only the strongest spacing peaks;
7. refine each peak in a small local neighborhood with modulo phase support from at most 128 high-occurrence seed groups and at most 1,024 evenly sampled base comparisons;
8. release the read-local seed structures before reading the next record.

The algorithm no longer compares every base for every possible period. Its Python pilot cost is dominated by per-read k-mer extraction plus bounded spacing/refinement work.

## Streaming And Interruption

At command start, discover creates:

1. `run.log`;
2. `run_config.yaml` with running status;
3. `candidate_reads.tsv` with its header.

Each accepted candidate is appended and flushed immediately. Progress logs report processed reads, processed bases, candidate reads, elapsed time, reads/s, MB/s and an estimated remaining time when `--max-reads` is set. Ctrl-C leaves existing candidate output available for diagnosis.

## Pilot Controls

Use these controls for real-read subsets:

```text
--max-reads
--max-read-bases
--sample-rate
--seed
--progress-every
--min-read-length
--min-period
--max-period
--kmer-size
--top-periods
--min-seed-occurrences
--min-spacing-support
--max-pairs-per-kmer
--chunk-size
```

`--chunk-size` currently defines a logical boundary for future parallel/checkpoint work. Resume is not implemented and the CLI does not expose a non-functional `--resume` flag.

## Backend Boundary

`--kmer-backend python` is the current toy/pilot implementation. Production-scale global k-mer counting should use optional mature backends such as KMC, meryl or Jellyfish instead of a new Python production counter. Read-local position-aware seed spacing is distinct from global counting; profiling indicates that its rolling extraction loop is the primary Rust/C++ migration target. These external or compiled backends are planned interfaces, not current dependencies.

## Remaining Limits

The spacing prefilter enables larger subsets but does not make TandemX ready for full 7-20 Gb genomes or 100 Gb read sets. Remaining work includes multiprocessing, chunk checkpoints, resumable execution, faster compiled seed extraction, portable peak-memory reporting and validation on real HiFi subsets.

## Synthetic Reference Measurements

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
