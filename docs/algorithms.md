# TandemX Algorithm Design

This document describes planned algorithms. The current repository skeleton does not implement these algorithms yet.

## Candidate Periodic k-mer Discovery

MVP goal: identify simple candidate tandem repeat monomers from toy HiFi-like reads.

Planned approach:

1. stream reads instead of loading all reads at once;
2. scan each read for periodic k-mer patterns within a configured monomer length range;
3. record candidate repeat spans and period estimates;
4. filter low-complexity candidates;
5. cluster similar candidate monomers into repeat families;
6. derive a simple consensus monomer for each toy-scale family.

MVP constraints:

1. toy data only;
2. simple tandem arrays only;
3. no higher-order repeat inference;
4. no production-scale memory optimization claims.

Future work:

1. robust clustering for related satellite families;
2. strand-aware consensus refinement;
3. parallel chunked read processing;
4. uncertainty modeling for ambiguous monomer periods.

## Diagnostic k-mer Copy-number Calibration

MVP goal: estimate repeat family copy number from diagnostic k-mer depth on toy reads.

Planned approach:

1. enumerate canonical k-mers from each consensus monomer;
2. exclude low-complexity and non-specific k-mers;
3. count selected diagnostic k-mers in reads;
4. summarize depth per family using robust statistics such as median depth;
5. normalize against an estimated genome-wide depth;
6. report estimated repeat bp and monomer copy number with confidence labels.

MVP constraints:

1. simple genome-size input;
2. no complex ploidy model;
3. confidence labels based on diagnostic k-mer count and depth dispersion.

Future work:

1. contamination-aware k-mer filtering;
2. depth modeling across multiple samples;
3. uncertainty intervals from bootstrap or Bayesian models;
4. scalable k-mer counting backends.

## Assembly Density Localization

MVP goal: locate simple repeat evidence on a toy assembly and summarize density in windows.

Planned approach:

1. scan assembly sequence for monomer matches or approximate matches;
2. write repeat hit intervals as 0-based half-open BED;
3. merge overlapping hits for density summaries;
4. compute covered bp per window;
5. write bedGraph density tracks for visualization.

MVP constraints:

1. toy assemblies only;
2. simple matching is acceptable;
3. bigWig output is future work.

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

Planned approach:

1. generate probe candidates from consensus monomers or concatenated monomer sequence;
2. compute probe length, GC fraction and simple melting temperature estimate;
3. evaluate target abundance from `copy_number.tsv`;
4. estimate specificity from k-mer sharing with non-target families;
5. summarize predicted signal regions from assembly localization when available;
6. combine abundance, specificity and localization evidence into a probe score.

MVP constraints:

1. score is a prioritization heuristic;
2. no guarantee of experimental FISH success;
3. no thermodynamic hybridization model.

Future work:

1. off-target search against full assemblies;
2. probe tiling and multiplex design;
3. empirical calibration against published and new FISH data.
