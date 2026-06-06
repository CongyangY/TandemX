# TandemX Algorithm Design

This document describes current MVP algorithms and planned future algorithms. The repository currently implements the toy simulator, toy-scale `discover`, `quantify`, `locate`, and `probe` MVPs.

## Candidate Periodic k-mer Discovery

MVP goal: identify simple candidate tandem repeat monomers from toy HiFi-like FASTA reads.

Current MVP implementation:

1. parse FASTA reads only;
2. for each read, test candidate periods from `--min-monomer-len` to `--max-monomer-len`, limited to at most half the read length;
3. compute a simple shifted periodicity identity score for each candidate period:
   `matches(sequence[i], sequence[i + period]) / compared_bases`;
4. keep the best period for a read when the score is at least 0.75;
5. write each retained read-level candidate to `candidate_reads.tsv`;
6. cluster candidates by period length using a small bp tolerance;
7. keep clusters with at least `--min-support-reads` distinct reads;
8. use the highest-scoring candidate sequence as the representative monomer for each family;
9. write `monomers.fa` and `families.tsv`.

This method is intentionally simple. It is designed to recover the 566 bp and 350 bp simulated repeat families from the toy dataset, not to analyze real large plant genomes.

MVP constraints:

1. toy data only;
2. simple tandem arrays only;
3. FASTA only; FASTQ support is future work;
4. one best candidate period per read;
5. no local repeat-boundary refinement;
6. no multiple-alignment consensus;
7. no higher-order repeat inference;
8. no production-scale memory optimization claims.

Future work:

1. robust clustering for related satellite families;
2. strand-aware consensus refinement;
3. parallel chunked read processing;
4. uncertainty modeling for ambiguous monomer periods.

## Diagnostic k-mer Copy-number Calibration

MVP goal: estimate repeat family copy number from diagnostic k-mer depth on toy FASTA reads.

Current MVP implementation:

1. read monomer FASTA from `--catalog/--catalogue` or `--monomers`;
2. enumerate canonical k-mers from each consensus monomer;
3. count k-mer multiplicity within each monomer;
4. remove low-complexity k-mers and k-mers shared by multiple families;
5. count canonical k-mers in the input reads;
6. correct each diagnostic k-mer depth by its multiplicity within the monomer;
7. summarize corrected diagnostic k-mer depth with the median;
8. if `--haploid-depth` is provided, use it directly;
9. otherwise estimate haploid depth as total read bases divided by `--genome-size` and add a warning;
10. estimate copy number as `median_kmer_depth / haploid_depth`;
11. estimate repeat bp as `estimated_copy_number * monomer_length`.

MVP constraints:

1. toy FASTA reads only;
2. no complex ploidy model;
3. no genome-wide unique k-mer depth estimation;
4. no external k-mer counter;
5. missing `--haploid-depth` uses a rough total-bases/genome-size estimate and is labeled with a warning;
6. confidence labels are based on diagnostic k-mer count and whether haploid depth was provided.

Future work:

1. contamination-aware k-mer filtering;
2. depth modeling across multiple samples;
3. uncertainty intervals from bootstrap or Bayesian models;
4. scalable k-mer counting backends.

## Assembly Density Localization

MVP goal: locate simple repeat evidence on a toy assembly and summarize density in windows.

Current MVP implementation:

1. read assembly FASTA and monomer FASTA;
2. enumerate non-low-complexity canonical k-mers from each monomer;
3. scan assembly sequence for matching canonical k-mers;
4. convert matching k-mers into intervals;
5. merge nearby intervals for each family;
6. filter very short intervals;
7. write candidate arrays as 0-based half-open `arrays.bed`;
8. compute merged-interval coverage per sliding window;
9. write `repeat_density.bedgraph`;
10. if `copy_number.tsv` is provided, compare read-estimated bp and assembly-estimated bp in `assembly_vs_read_cn.tsv`.

MVP constraints:

1. toy assemblies only;
2. k-mer evidence only, no read mapping;
3. no exact per-copy placement claim;
4. simple threshold classification for `possible_collapse` and `possible_overexpansion`;
5. bigWig output is future work.

Future work:

1. alignment-backed localization;
2. chromosome-scale streaming;
3. family-specific tracks;
4. confidence labels for ambiguous or multi-family hits.

## Assembly-vs-read Comparison

MVP goal: compare toy read-based and assembly-based repeat abundance.

Planned approach:

1. convert read-based copy number into estimated repeat bp;
2. summarize assembly repeat bp from localization outputs;
3. calculate read-to-assembly ratio;
4. classify status as consistent, possible under-representation, possible over-expansion or uncertain;
5. preserve warnings from input estimates.

MVP constraints:

1. simple threshold-based classification;
2. no claim of proof of collapse;
3. no complex assembly quality model.

Future work:

1. chromosome-level discrepancy summaries;
2. support for polyploid subgenome interpretation;
3. uncertainty propagation from read and assembly estimates.

## FISH Probe Scoring

MVP goal: rank simple candidate probes for toy repeat families.

Current MVP implementation:

1. read monomer FASTA, assembly FASTA, `copy_number.tsv`, and `arrays.bed`;
2. generate fixed windows from each monomer using `--min-len` and `--max-len`;
3. exclude probes with high single-base low-complexity ratio;
4. compute probe length, GC content and a simple Wallace-rule Tm estimate;
5. scan assembly with probe k-mers to find predicted target and off-target regions;
6. treat hits overlapping same-family `arrays.bed` intervals as target-array hits;
7. estimate `arrayiness_score` as target hits divided by all predicted hits;
8. estimate `specificity_score` as `1 / (1 + off_target_hits)`;
9. combine normalized copy number, specificity, arrayiness and GC balance into `probe_score`;
10. write ranked probes and toy in silico FISH predicted signal regions.

MVP constraints:

1. score is a prioritization heuristic;
2. no guarantee of experimental FISH success;
3. no thermodynamic hybridization model;
4. no full off-target alignment;
5. no probe tiling optimization;
6. no experimental calibration.

Future work:

1. off-target search against full assemblies;
2. probe tiling and multiplex design;
3. empirical calibration against published and new FISH data.
