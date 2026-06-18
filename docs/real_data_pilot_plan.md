# TandemX Real Data Pilot Plan

This plan defines staged real-data pilots before TandemX is considered for full large-genome production analysis. It does not authorize direct analysis of full 7-20 Gb plant genomes or 100 Gb raw read sets with the current MVP.

## Phase 1: Small Read Subsets

Goal: test whether de novo discovery produces plausible candidate families on real HiFi read subsets.

Inputs:

1. 10,000 reads sampled from a HiFi `.fastq.gz`.
2. 50,000 reads sampled from the same dataset.
3. 100,000 reads sampled from the same dataset.

Run only:

```bash
tandemx discover \
  --reads sampled_reads.fastq.gz \
  --outdir discover_subset \
  --max-reads 10000 \
  --progress-every 1000
```

Start with 10,000 reads. Increase to 50,000 only after reviewing runtime, candidate rate and logs. Treat 100,000 reads as a manual stress pilot until the Python backend has measured acceptable throughput.

Checks:

1. `candidate_reads.tsv`, `monomers.fa` and `families.tsv` are produced.
2. candidate family counts are plausible, not empty and not dominated only by low-complexity sequences.
3. logs and `run_config.yaml` are complete.
4. `tandemx validate --project discover_subset` passes for recognized outputs.

Do not make final copy-number, assembly-collapse or biological abundance conclusions from this phase.

## Phase 2: Known Repeat Sanity Check

Known satellite sequences such as pSc200, pSc250 or CentT566 may be used only for interpretation after de novo discovery.

Allowed use:

1. compare discovered `monomers.fa` against curated sequences after discovery;
2. check whether the de novo catalog contains similar sequence lengths or motifs;
3. report matches as sanity checks with caveats.

Not allowed in the default workflow:

1. using curated sequences as `tandemx discover` input;
2. presenting guided recovery as de novo discovery;
3. claiming validation from sequence similarity alone.

## Phase 3: Assembly Region Localization

Goal: test whether localization outputs are structurally reasonable on a constrained assembly subset.

Inputs:

1. one chromosome; or
2. one 50-200 Mb assembly region.

Run:

```bash
tandemx locate \
  --assembly assembly_region.fa \
  --catalog discover_subset/monomers.fa \
  --outdir locate_region
```

Checks:

1. `repeat_density.bedgraph` is non-empty and uses valid 0-based half-open intervals.
2. `arrays.bed` intervals are plausible and validate cleanly.
3. density tracks are visually and numerically inspected for obvious artifacts.

This phase still does not establish whole-genome production readiness.

## Phase 4: Whole-genome Preconditions

Before any full genome run, TandemX must have:

1. runtime and memory benchmarks on synthetic scale tests;
2. chunked read and assembly processing;
3. resumable runs with restartable intermediate outputs;
4. robust command logs and `run_config.yaml` for every step;
5. `tandemx validate` coverage for all generated outputs;
6. clear failure behavior for partial or malformed inputs;
7. documented pilot results on read subsets and assembly regions.

Do not directly run the current MVP on 100 Gb raw reads or full 7-20 Gb assemblies.
