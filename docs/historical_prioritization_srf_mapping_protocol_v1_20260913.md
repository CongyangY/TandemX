# Frozen historical-assembly SRF and mapping prioritization protocol v1

Status: **preregistered before formal historical native comparator execution**  
Freeze date: 2026-09-13  
Execution gate: main-task review plus the runner acknowledgement string
`HISTORICAL_PROTOCOL_V1_REVIEWED`  
Protocol scope: benchmark-only; no TandemX production algorithm is run or changed.

## 1. Question and interpretation boundary

This protocol asks whether two native comparators prioritize the same frozen
Ey15-2 and *Macadamia jansenii* TandemX families as a retrospective
historical-assembly reference proxy. The proxy-positive label is the already
frozen condition

`old_assembly_bp / new_assembly_bp < 0.6`

among families with `new_assembly_bp >= 15,000`. It is a newer-assembly proxy,
not absolute copy-number truth. Confusion matrices therefore measure agreement
with that proxy and must not be described as biological or physical accuracy.
Magnitude estimates remain read--assembly abundance deficits.

The methods are exactly:

1. SRF with `k=151`, `ci=20` (primary SRF cell);
2. SRF with `k=101`, `ci=20` (predeclared sensitivity cell);
3. ordinary competitive mapping to the complete frozen TandemX catalogue.

No native comparator result existed in the endpoint when these rules were
selected. Existing direction-only Macadamia mapping artifacts, including any
old PAF, are excluded. The runner has no option to import a PAF or abundance
table.

## 2. Immutable input binding

The machine-readable configuration
`benchmarks/configs/historical_prioritization_srf_mapping_v1_20260913.json`
binds every accepted file by absolute path, byte size, and SHA-256. Formal
execution must recompute all hashes before creating the result directory.

### Ey15-2

- Reads: complete `ERR8666125.fastq.gz`, 837,586 reads and
  18,636,790,429 bases.
- Frozen catalogue: all 2,133 records in the frozen donor-matched
  `discover/monomers.fa`.
- Frozen endpoint table: the primary total-bases `family_metrics.tsv`; exactly
  19 rows must be `eligible`.
- Genome-size denominator: 143,120,000 bp, inherited from the frozen source
  configuration.

### *Macadamia jansenii*

- Reads: the two complete FASTQ files in this exact order:
  `SRR13557763_subreads.fastq.gz`, then `SRR13557762_subreads.fastq.gz`.
  Together they contain 1,642,394 reads and 22,546,488,654 bases.
- Every native method consumes both files in that order. No concatenated,
  sampled, or one-file surrogate is allowed.
- Frozen catalogue: all 1,227 records in the frozen donor-matched
  `discover/monomers.fa`.
- Frozen endpoint table: the primary total-bases `family_metrics.tsv`; exactly
  43 rows must be `eligible`.
- Genome-size denominator: 780,000,000 bp, inherited from the frozen source
  configuration.

The full catalogue is the competitive search space for both methods. The
19/43 eligible subsets are applied only after native measurement and
correspondence. Mapping only eligible families would remove competitors and is
forbidden.

## 3. Fixed native workflows

All six cells run serially with one external-tool thread. Commands and their
argument order are written before launch.

### 3.1 SRF k151 and k101

For each species and k value, the workflow is:

1. Write `reads.list` with one absolute FASTQ path per line in configuration
   order.
2. `kmc -fq -k<K> -t1 -m2 -sm -ci20 -cs100000 @reads.list counts tmp`
3. `kmc_dump counts counts.txt`
4. `srf -p srf counts.txt > srf.fa`
5. `k8 srfutils.js enlong srf.fa > srf.enlong.fa`
6. `minimap2 -c -N1000000 -f1000 -r100,100 -t1 srf.enlong.fa <FASTQ...> > srf.paf`
7. `k8 srfutils.js paf2bed srf.paf > srf.bed`
8. `k8 srfutils.js bed2abun -g <TOTAL_LIBRARY_BASES> srf.bed > srf.abundance.tsv`

KMC receives gzip FASTQ explicitly through `-fq` and its native `@file-list`
interface. The second column of `srf.abundance.tsv` is the native retained read
base count. The runner verifies it against the sum of positive-flag intervals
in `srf.bed`. Columns four and five are not used as independent measurements.

An empty KMC dump is `no_eligible_kmers`; an empty SRF catalogue is
`no_catalogue`. Both are successful native zero-discovery states, not process
failures.

### 3.2 Ordinary competitive mapping

The runner builds one 10,000-bp tandemized template for every record in the
complete frozen TandemX catalogue. It then runs:

`minimap2 -x map-hifi -c -N1000000 -f1000 -r100,100 -t1 templates.fa <FASTQ...> > mapping.paf`

The frozen simulation rule is then applied unchanged:

- retain primary and secondary PAF rows;
- require alignment block length at least 100 bp;
- require `matches / alignment_block >= 0.90`;
- union intervals within each read/family;
- exclude each query base covered by two or more families from every family
  estimate and report it once as `ambiguous_native_bp`.

Discovery and family recovery are `N/A_shared_frozen_TandemX_catalogue` for
this method.

## 4. Family and HOR correspondence

SRF native motifs are compared with the complete frozen TandemX catalogue by
the frozen independent correspondence implementation in
`benchmarks/challenge/unified_correspondence.py`:

- cyclic rotations and reverse complements are allowed;
- Needleman--Wunsch edit identity must be at least 0.90;
- shorter and longer motifs may correspond by an integer repeat multiple;
- the maximum accepted multiple is 64;
- length disagreement beyond 10% of the longer sequence is rejected.

A native motif is assigned only when exactly one TandemX family passes.
Statuses `ambiguous`, `unmatched`, `unmatched_or_composite`, and `out_of_scope`
remain unassigned. Their counts and retained read bases are reported by status;
they are never redistributed to the best hit. Multiple uniquely assigned SRF
motifs may contribute to the same TandemX family because they represent native
SRF catalogue output rather than post-hoc family selection.

Ordinary mapping already uses TandemX family identifiers, so no sequence
correspondence step is applied after mapping.

## 5. Depth normalization and endpoint scoring

For either method and each assigned family:

`read_estimated_bp = native_retained_read_bp / total_library_bases * genome_size_bp`

The total-library and genome-size values are fixed per species and are not
estimated from comparator output. No HiFi quality correction, haploid-depth
fit, or calibration against the newer assembly is applied.

For every frozen eligible family with a positive estimate:

`predicted_proxy_positive = old_assembly_bp / read_estimated_bp < 0.6`

The label is recomputed from frozen old/new assembly values and checked against
the frozen table. Summary rows report TP/FN/FP/TN only over families with a
defined prediction. The eligible denominator, available denominator, and
unavailable count are always separate.

## 6. Zero, unmatched, unavailable, and failure semantics

- **zero native support**: the method completed, but no uniquely assigned
  retained read bases support the eligible family. Store numeric
  `read_estimated_bp=0`; `old/read` and the binary prediction are `N/A` because
  the ratio denominator is zero. Exclude the row from TP/FN/FP/TN. Count it in
  `zero_supported_eligible_families` and in method-specific family recovery.
- **unmatched/ambiguous/out-of-scope SRF motif**: retain the native motif and
  its raw and normalized abundance in the unassigned output. Do not convert it
  to a family-level false positive and do not give its bases to an eligible
  family.
- **N/A**: used only for a mathematically undefined quantity, an endpoint that
  does not apply, or a technical cell failure. It is never serialized as zero.
- **process failure/timeout/resource stop**: preserve prior stages and logs;
  write every eligible family with `read_estimated_bp=N/A`, prediction `N/A`,
  and the technical state. Do not reuse output from another attempt.
- **unresolved frozen reference row**: would remain `N/A`, but formal input
  validation requires every selected eligible row to be either
  `reference_collapse` or `reference_retained`.

Zero recovery is a method result. Technical failure is not.

## 7. Resource and stop rules

- cells and stages are serial; `threads=1` everywhere;
- KMC strict memory is 2 GiB (`-m2 -sm`);
- per-stage wall limits are fixed in the configuration (8 h counting/mapping,
  2 h dump/assembly/filter, 1 h elongation/abundance);
- the complete six-cell launch deadline is 96 h;
- at least 250 GB free must remain on the output volume;
- a cell may occupy at most 500 GB after any completed stage;
- no later cell launches after a failure or a breached limit;
- a native process is killed on its wall timeout; all prior evidence remains;
- no parameter, threshold, input, or method order may change after output
  inspection. A changed protocol requires a new version and all cells rerun.

The disk limits are checked before launch and after every stage. This is an
explicit conservative stop rule; it is not a claim that peak temporary disk is
measured continuously.

## 8. Required receipts and outputs

The result root must be new and empty. It contains:

- `run_config.json`: exact configuration snapshot and its SHA-256;
- `environment.json`: platform, Git state, Python, tool paths and hashes;
- `input_manifest.json`: all input sizes/hashes and validated read ordering;
- `source_snapshot/`: this protocol, runner, correspondence and PAF parser;
- `plan.tsv`: fixed six-cell order;
- per-cell `reads.list`, `*.command.json`, stdout/stderr logs,
  `*.receipt.json`, native files, catalogue correspondence, assigned and
  unassigned abundance, family estimates, family prioritization, and
  `result.json`;
- `summary.tsv`: cell status, denominators, confusion counts, zero support,
  unassigned abundance, runtime and resources;
- `completion.json` and `archive_manifest.json` with hashes of every final
  artifact.

Each stage receipt records command, start/end UTC, exit code, timeout state,
wall/CPU/RSS measurement, output byte count and SHA-256, cell disk usage and
free disk after the stage. A cell is resumable only from a finalized result
whose complete artifact hashes, input hashes, tool hashes, configuration hash,
protocol hash, and source hashes still match. The v1 runner deliberately
requires a new output root instead of partial in-place resume; failed attempts
are retained as evidence.

## 9. Execution gate

Before review, only structural planning and unit tests may run:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.run_historical_prioritization_srf_mapping \
  --config benchmarks/configs/historical_prioritization_srf_mapping_v1_20260913.json \
  --plan
```

After main-task review, the formal run requires both flags:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.run_historical_prioritization_srf_mapping \
  --config benchmarks/configs/historical_prioritization_srf_mapping_v1_20260913.json \
  --outdir /Volumes/T7/Codex/TandemX/benchmarks/results/historical_prioritization_srf_mapping_v1 \
  --execute --ack HISTORICAL_PROTOCOL_V1_REVIEWED
```

This acknowledgement is an execution latch, not evidence that the result is
valid. Completion still requires all receipts and zero failed cells.
