# TandemX Benchmark Plan

## Current Synthetic Harness

TandemX now includes a synthetic benchmark runner:

```bash
python benchmarks/scripts/run_synthetic_benchmark.py \
  --config benchmarks/configs/synthetic_scale.yaml \
  --scale tiny \
  --outdir /tmp/tandemx_benchmark_tiny
```

Configured scales:

1. `tiny`: 1,000 reads for CI.
2. `small`: 10,000 reads for manual runtime measurement.
3. `pilot`: 50,000 reads for manual subset scaling.
4. `real_pilot_manual`: 100,000 reads, excluded from default tests and commands.

The runner records wall-clock runtime, peak resident memory from `wait4`, and discover throughput metrics in `benchmark_summary.tsv`: processed reads/bases, candidate reads, candidates/MB, reads/s, MB/s and `algorithm_mode=spacing_prefilter`.

`accuracy_summary.tsv` uses simulator truth files only after analysis commands finish. Recovered families are paired to truth by length-aware sequence identity rather than monomer length alone. Truth files are benchmark evaluation metadata and must not be passed to `discover`, `quantify`, `locate`, `probe` or `validate`.

Synthetic benchmark success does not validate real 7-20 Gb plant genome production analysis.

## Purpose

Benchmarks should test whether TandemX can identify candidate repeat families, estimate read-based copy number, detect possible assembly under-representation and rank FISH probe candidates. Benchmarks must separate toy MVP validation from real plant genome claims.

## Compared Tools

Benchmark comparisons should include TRF, TideHunter, TRASH and RepeatExplorer2/TAREAN when applicable. Comparisons must be task-specific because these tools have different assumptions and outputs.

The implemented external comparison harness covers the directly comparable
read-level discovery subset for TandemX, TRF, and TideHunter. It uses the same
synthetic FASTA inputs and one thread per tool. TRASH and
RepeatExplorer2/TAREAN remain separate task-specific benchmarks because their
assembly and low-coverage short-read graph workflows are not interchangeable
with per-read tandem-repeat detection.

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
| selected_source_query_recovery | fraction | Source-backed repeat queries recovered from actual consensus outputs; call this known-family recall only when independent same-material presence and a declared family denominator support it |
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

## Divergence-aware localization validation

The first IID-proxy development run on seeds 5301--5303 failed its predeclared
full-assembly mean base-recall gate (0.903608 < 0.95). Development v2 reused only
those consumed development seeds, retained the 0.90 IID proxy threshold and
bridged exact anchors across at most one monomer length. It passed the three
development gates: full-assembly mean base recall 0.975728, positive-assembly
mean precision 0.999441 and absent-family false-positive rate 0/54.

`abundance_localizer_heldout_v1.json` freezes that localizer, its development
artifact hashes, the earlier multi-k collapse rule and untouched seeds
5501--5503. The runner re-hashed and recomputed the development gate metrics
before the first IID held-out run. All 1,062 commands completed. Held-out
full-assembly recall was 0.981033, positive-assembly precision was 0.999565 and
absent-family false-positive rate was 0/54. The frozen multi-k classifier raised
sensitivity from 0.600137 to 0.918381 while raising false-positive rate from
0.012346 to 0.085391 and lowering precision from 0.986471 to 0.941632. The
complete declared matrix and adverse rows are retained without retuning on
seeds 5501--5503. This is a known-catalogue substitution model; exact-k-mer IID
identity is not alignment identity or biological satellite validation.

The next classifier development is predeclared in
`abundance_classifier_development_v1.json`. Seeds 5601--5603 are development
data; 5701--5703 are reserved and must remain untouched until a selected rule is
frozen. The candidate set is a 5-by-4 grid of log-space single/multi-k blends
and decision thresholds. Selection maximizes sensitivity only among candidates
with no worse false-positive rate or precision than single k=21, with fixed
tie-breakers and leave-one-seed-out stability gates. The alpha-zero, threshold-
0.6 candidate exactly anchors the current baseline. This transparent experiment
is not an AI model and must be retained as a failed development result if its
predeclared gates are not met.

## Endpoint and clustering audit (2026-09-06)

Retain both strict equal-length recovery and independent cyclic edit recovery;
small consensus indels are not automatically missing families. Report raw array
precision together with period-independent base unions and duplicate-base burden.
The 95% operational monomer-cluster threshold is a declared sequence resolution;
its sensitivity sweep is separate from the 90% planted-sequence recovery endpoint.
Seed 2101 was inspected during repair and cannot be presented as untouched test
performance. Seeds 3101/3102/3103 remain unused.

SRF requires its whole KMC/count-dump/circle-assembly/elongation/minimap2/filter/
abundance workflow to be measured. Preliminary ci100 and ci20 settings are an
illustrative README preset and a lower-count sensitivity setting; the synthetic
library has no known average genome coverage. Record both, including no-catalogue
outcomes. Native higher-order units must remain available and cannot be called
false monomers simply because they exceed a period-restricted score. A dedicated
HOR decomposition endpoint is still needed before interpreting that difference.
