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
7. refine each peak in a small local neighborhood with modulo phase support from at most 128 high-occurrence seed groups and at most 4,096 sampled base comparisons;
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

`--kmer-backend python` is the current toy/pilot implementation. Production-scale global k-mer counting should use optional mature backends such as KMC, meryl or Jellyfish instead of a new Python production counter. These external backends are planned interfaces, not current dependencies.

## Remaining Limits

The spacing prefilter enables larger subsets but does not make TandemX ready for full 7-20 Gb genomes or 100 Gb read sets. Remaining work includes multiprocessing, chunk checkpoints, resumable execution, faster compiled seed extraction, portable peak-memory reporting and validation on real HiFi subsets.

## Synthetic Reference Measurements

Measurements on the development workstation on 2026-06-18:

| Scale | Reads | Read bp | Discover runtime | Throughput | Families |
|---|---:|---:|---:|---:|---:|
| tiny | 1,000 | 1.8 Mb | 3.85 s | 266.5 reads/s; 0.480 MB/s | 2 |
| small | 10,000 | 18 Mb | 30.23 s | 349.5 reads/s; 0.629 MB/s | 2 |

These synthetic reads are 1.8 kb, much shorter than 16 kb HiFi reads. Real-read runtime must be estimated from processed bases, not read count alone.
