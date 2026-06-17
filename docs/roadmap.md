# TandemX Roadmap

## Near-term MVP Completion

1. Split `assembly_vs_read_cn.tsv` comparison logic into a dedicated `tandemx compare` implementation.
2. Expand output validation toward cross-file consistency checks.
3. Add richer confidence and warning propagation across modules.
4. Add versioned toy benchmark configs.
5. Keep randomized toy workflow tests small, reproducible and independent of simulator truth files as algorithm input.

## Large-genome Optimization Path

Future real plant genome support will require:

1. streaming algorithm design beyond the parser layer;
2. chunked assembly scanning;
3. scalable k-mer counting;
4. multiprocessing with bounded memory;
5. resumable intermediate files;
6. indexed monomer and assembly search;
7. robust logging and resource reporting.

FASTA, FASTQ and gzip input parsing is now available for the MVP. This does not remove the need for streaming optimization, parallelization, bounded-memory k-mer counting, external benchmarking, and real read validation before 7-20 Gb plant genome analysis.

## Possible Rust or C++ Rewrite Targets

Python should remain the prototype layer until algorithms stabilize. Future compiled components may be useful for:

1. high-throughput k-mer counting;
2. large-read periodicity scanning;
3. assembly k-mer indexing;
4. interval merging at chromosome scale;
5. probe off-target scanning.

## Benchmark Roadmap

Simulated benchmark:

1. monomer length accuracy;
2. family recall;
3. false positive rate;
4. copy-number error;
5. assembly-collapse detection;
6. runtime and peak memory.

Real plant benchmark:

1. wheat, rye, barley, oat and maize datasets;
2. known centromeric repeats;
3. subtelomeric repeats;
4. chromosome-enriched repeats;
5. published FISH probe-associated repeats.

Comparison tools:

1. TRF;
2. TideHunter;
3. TRASH;
4. RepeatExplorer2/TAREAN.

## Real Data Integration Plan

Rye, wheat and triticale integration should start only after simulated benchmarks are stable.

Required inputs:

1. public HiFi reads;
2. matching assembly versions;
3. known satellite repeat annotations;
4. published FISH probes or cytogenetic signal descriptions;
5. clear genome size and depth metadata.

## Publication Evidence Gaps

Bioinformatics:

1. robust software benchmark;
2. reproducible workflows;
3. clear comparison with existing tools;
4. documented file formats and tests.

Genome Biology or Genome Research:

1. stronger real-genome biological case studies;
2. assembly representation analysis;
3. repeat biology interpretation;
4. validation against known annotations and literature.

Plant Communications or Molecular Plant:

1. plant-focused biological story;
2. strong wheat, rye, barley, oat or maize examples;
3. FISH probe recovery or experimental validation;
4. relevance to chromosome biology, genome evolution or breeding resources.
