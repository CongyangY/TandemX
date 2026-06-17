# TandemX Benchmark Plan

## Purpose

Benchmarks should test whether TandemX can identify candidate repeat families, estimate read-based copy number, detect possible assembly under-representation and rank FISH probe candidates. Benchmarks must separate toy MVP validation from real plant genome claims.

## Compared Tools

Benchmark comparisons should include TRF, TideHunter, TRASH and RepeatExplorer2/TAREAN when applicable. Comparisons must be task-specific because these tools have different assumptions and outputs.

## Stage 1: Simulated Reads Benchmark

Goal: evaluate candidate monomer discovery and diagnostic k-mer copy-number calibration on controlled reads.

Inputs:

1. simulated HiFi-like reads;
2. truth repeat families;
3. truth monomer lengths;
4. truth copy numbers.

Metrics:

| Metric | Unit | Definition |
|---|---:|---|
| family_recall | fraction | Truth families recovered by the tool |
| false_positive_rate | fraction | Predicted families without truth support divided by predicted families |
| monomer_length_error_bp | bp | Absolute difference between predicted and truth monomer length |
| copy_number_relative_error | fraction | Absolute estimated-vs-truth copy-number difference divided by truth |
| reproducibility | boolean | Whether fixed seed gives stable output |

## Stage 2: Simulated Assembly Collapse Benchmark

Goal: evaluate assembly-vs-read comparison under known assembly representation errors.

Inputs:

1. simulated reads;
2. simulated assemblies;
3. truth assembly collapse or expansion labels;
4. truth repeat bp in reads and assemblies.

Metrics:

| Metric | Unit | Definition |
|---|---:|---|
| assembly_collapse_detection_accuracy | fraction | Fraction of simulated cases correctly classified |
| under_assembly_precision | fraction | True under-assembly calls divided by all under-assembly calls |
| under_assembly_recall | fraction | True under-assembly calls recovered |
| read_to_assembly_ratio_error | ratio | Difference between observed and truth read/assembly ratio |

## Stage 3: Real Plant Genome Benchmark

Goal: evaluate behavior on real plant genome datasets after the toy and simulated stages are stable. This is not part of the MVP.

Candidate species:

1. wheat;
2. rye;
3. barley;
4. oat;
5. maize.

Required metadata:

| Field | Definition |
|---|---|
| species | Species name |
| accession | Accession or cultivar |
| genome_size_bp | Estimated genome size |
| read_type | HiFi, CLR, ONT or other |
| assembly_version | Assembly accession or version |
| known_repeats | Literature-supported repeat families |
| fish_validation | Whether FISH validation exists |
| citation | Dataset or publication source |

Metrics:

| Metric | Unit | Definition |
|---|---:|---|
| known_family_recovery | fraction | Known repeat families recovered |
| known_location_consistency | qualitative | Agreement with published chromosomal localization |
| runtime_sec | seconds | Wall-clock runtime |
| peak_memory_mb | MB | Peak resident memory |
| disk_usage_mb | MB | Output and intermediate disk usage |

## Stage 4: FISH Validation Benchmark

Goal: test whether probe ranking is useful for FISH probe prioritization.

Validation sources:

1. published FISH probes;
2. published satellite repeat locations;
3. known centromeric or subtelomeric repeats;
4. new experimental FISH validation if available.

Metrics:

| Metric | Unit | Definition |
|---|---:|---|
| known_probe_recovery | fraction | Known probe families or sequences recovered |
| probe_specificity | fraction | Predicted target signal relative to predicted off-target signal |
| off_target_count | hits | Number of predicted off-target regions |
| predicted_signal_consistency | qualitative | Agreement with published FISH signal regions |
| validated_probe_rank | rank | Rank of experimentally validated probe |

## Tool Comparison Metrics

TRF:

1. tandem interval detection;
2. monomer length accuracy;
3. runtime and memory on sequence-level inputs.

TideHunter:

1. long-read monomer recovery;
2. family recall from read-level candidates;
3. monomer length accuracy.

TRASH:

1. satellite family recovery;
2. curated satellite recovery in plant genomes;
3. runtime and memory.

RepeatExplorer2/TAREAN:

1. known satellite family recovery;
2. biological consistency with literature;
3. comparison against graph/clustering-based repeat characterization.

TandemX-specific metrics:

1. diagnostic k-mer copy-number error;
2. assembly-vs-read discrepancy classification;
3. probe specificity and validated probe rank;
4. output field completeness;
5. reproducibility.

## Claims Policy

After toy benchmarks, claim only that the toy workflow runs and produces documented outputs.

After simulated benchmarks, claim only behavior under tested simulated conditions.

After real plant and FISH benchmarks, claim biological usefulness only for tested species, datasets and validation scenarios.
