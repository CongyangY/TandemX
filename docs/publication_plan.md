# TandemX Publication Plan

This document lists evidence needed before targeting different journals. It does not imply that the current skeleton has produced any result.

## Bioinformatics

Likely emphasis: method, software design, reproducibility and benchmarked utility.

Required evidence:

1. documented CLI and file formats;
2. unit and integration tests;
3. simulated read benchmark;
4. simulated assembly collapse benchmark;
5. comparison with TRF, TideHunter, TRASH and RepeatExplorer2/TAREAN where applicable;
6. runtime and memory profiling;
7. at least one reproducible toy tutorial and one realistic small benchmark workflow.

FISH validation is helpful but may not be mandatory if the paper is framed as a computational method with clear limitations.

## Genome Biology

Likely emphasis: method plus biological insight across important systems.

Required evidence:

1. strong simulated and real-genome benchmarks;
2. multiple real plant genome case studies;
3. recovery of known centromeric, pericentromeric or subtelomeric repeats;
4. clear examples where read-vs-assembly comparison reveals possible under-representation;
5. external biological validation, preferably including FISH evidence;
6. robust comparison with existing repeat-analysis tools.

## Genome Research

Likely emphasis: genome biology, assembly interpretation and evidence depth.

Required evidence:

1. rigorous assembly-vs-read copy-number comparison;
2. case studies showing how satellite arrays are represented or under-represented in assemblies;
3. validation against high-quality assemblies and curated satellite annotations;
4. careful uncertainty labels;
5. biological interpretation for real plant genomes;
6. comparison to established repeat discovery workflows.

## Plant Communications

Likely emphasis: plant genomics application and practical biological utility.

Required evidence:

1. plant-focused benchmark datasets;
2. wheat, rye, barley, oat or maize case studies;
3. recovery of known plant satellite repeats;
4. FISH probe prioritization examples;
5. at least published FISH consistency, ideally new validation;
6. accessible tutorial and reproducible workflow.

## Molecular Plant

Likely emphasis: strong plant biology story and experimental validation.

Required evidence:

1. compelling plant genome case studies;
2. experimentally validated FISH probes or strong independent cytogenetic evidence;
3. biological insight into centromeric, pericentromeric or subtelomeric repeats;
4. high-confidence comparison between read-based copy number and assembly representation;
5. clear relevance to plant genome evolution, breeding resources or chromosome biology.

## Required Benchmark Evidence

Minimum computational evidence before manuscript submission:

1. simulated reads benchmark;
2. simulated assembly collapse benchmark;
3. real plant genome benchmark;
4. runtime and memory benchmark;
5. comparison with TRF, TideHunter, TRASH and RepeatExplorer2/TAREAN;
6. reproducibility benchmark with fixed seeds and versioned configs.

## Current Advantage Claims

The 2026-07-15 external-tool benchmark supports the following provisional,
scope-limited engineering claim:

> On controlled single-family synthetic long-array reads, single-threaded
> TandemX retained the same period-recovery accuracy as TRF while achieving
> higher throughput, and used less peak memory than TideHunter.

The current evidence does not support saying that TandemX is generally more
accurate than TRF or TideHunter, always uses the least memory, or outperforms
TRASH or TAREAN/RepeatExplorer2. TRASH and TAREAN require task-specific
comparisons because they target assembly hierarchy/HOR annotation and
short-read graph-based family recovery, respectively.

The stronger publication-level differentiation is the integrated read-first,
assembly-aware workflow: candidate monomer discovery, diagnostic k-mer copy
number estimation, assembly localization, possible assembly
under-representation checks, FISH probe prioritization, confidence labels and
reproducible outputs. Biological utility claims for that integration still
require independent real HiFi data, known-repeat or expert truth, and FISH or
other orthogonal validation.

## Required Biological Case Studies

Potential case-study categories:

1. centromeric repeat family in a cereal genome;
2. subtelomeric repeat family with chromosome-enriched signal;
3. species-specific or accession-enriched repeat family;
4. assembly under-representation example supported by read-vs-assembly discrepancy;
5. known FISH probe family recovered by TandemX.

## Required FISH Validation

Validation options:

1. recovery of published FISH probe sequences or repeat families;
2. agreement with published signal regions;
3. independent off-target assessment against assembly or repeat catalogue;
4. new FISH experiment for top-ranked candidate probes if resources allow.

Claims must remain limited to tested datasets and validation evidence. Probe ranking should be described as prioritization, not guaranteed experimental success.
