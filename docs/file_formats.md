# TandemX File Formats

All TSV files must use tab separators, include a header line and use stable column names. Fields with uncertain interpretation should include `confidence` or `warning` where practical.

All BED files use 0-based half-open coordinates: `start` is included and `end` is excluded.

The current repository implements the toy simulator and toy-scale `discover`, `quantify`, `locate`, `probe`, and `visualize` MVPs. The formats below define current core outputs.

The analysis file flow is de novo by default:

```text
raw reads -> tandemx discover -> monomers.fa/families.tsv -> downstream commands
```

`monomers.fa` and `families.tsv` are discovery outputs, not prerequisite inputs for `tandemx discover`. Simulator truth files are documented only for toy simulation and tests.

## Sequence Inputs

TandemX uses `tandemx.io.sequences` for streaming sequence input where possible. Supported extensions are:

1. `.fa`
2. `.fasta`
3. `.fq`
4. `.fastq`
5. `.fa.gz`
6. `.fasta.gz`
7. `.fq.gz`
8. `.fastq.gz`

Each parsed record is normalized as `SequenceRecord(id, sequence, quality=None, description="...")`. FASTQ records retain quality strings and must have matching sequence and quality lengths. The reader reports clear errors for empty files, malformed FASTA/FASTQ syntax, duplicate record IDs, unsupported extensions, and unsupported bases.

Commands that accept `--reads` may receive one or more read files. Multiple
files are streamed in the order supplied. Discovery requires globally unique
read identifiers and uses an exact duplicate tracker that spills to a temporary
SQLite index after a bounded in-memory threshold. Quantification does not use
read IDs analytically, so IDs may repeat across distinct files, but a duplicate
path argument is rejected and duplicates within each file remain invalid.

## candidate_reads.tsv

Produced by: `tandemx discover`

This is a table of read-local evidence. Downstream catalogue arguments use
`monomers.fa`; this table itself contains no monomer sequence.
`read_id` may occur on several rows in elastic mode; `candidate_id` is unique.

| Field | Type | Unit | Description |
|---|---|---:|---|
| read_id | string | NA | Input read identifier |
| candidate_id | string | NA | Candidate repeat span identifier |
| read_start | integer | bp | 0-based start in read |
| read_end | integer | bp | 0-based half-open end in read |
| strand | string | NA | `+`, `-` or `.` |
| period_bp | integer | bp | Estimated period; observed-unit consensus length in elastic mode |
| repeat_span_bp | integer | bp | Length of candidate repeat span |
| unit_count | float | copies | Approximate repeat unit count |
| score | float | unitless | Legacy: weighted seed/shifted-identity score. Elastic: matches / alignment columns including gaps |
| low_complexity_flag | boolean | NA | Whether candidate is low complexity |
| confidence | string | NA | `high`, `medium` or `low` |
| warning | string | NA | Semicolon-separated warnings or empty, such as `short_period_candidate` or `low_complexity_candidate` |

Elastic warnings include `elastic_alignment`, `uncalibrated_confidence`,
`consensus_units=N` (at most 32 sampled observed units), and
`alignment_median_offset=N` when that offset differs from consensus length.
Neither the confidence label nor the composition filter is a calibrated
probability. Coordinates are read-local estimates, not assembly placements.

## monomers.fa

Produced by: `tandemx discover`

Header format:

```text
>family_id=TXF000001;monomer_id=TXM000001;length_bp=156;confidence=high
ACGT...
```

| Header Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Repeat family identifier |
| monomer_id | string | NA | Monomer identifier |
| length_bp | integer | bp | Monomer length |
| confidence | string | NA | Confidence label |

## families.tsv

Produced by: `tandemx discover`

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Repeat family identifier |
| monomer_id | string | NA | Representative monomer identifier |
| monomer_length_bp | integer | bp | Consensus monomer length |
| consensus_md5 | string | NA | MD5 hash of monomer sequence |
| gc_fraction | float | fraction | GC fraction from 0 to 1 |
| support_read_count | integer | reads | Number of supporting reads |
| support_span_bp | integer | bp | Total supporting read span |
| mean_identity | float | fraction | Historical field name: arithmetic mean of candidate scores, not a fresh candidate-to-family alignment. Elastic scores are adjacent-copy alignment identities; legacy scores also include seed support |
| low_complexity_flag | boolean | NA | Whether family is low complexity |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty, such as `low_complexity_family` |

Families with possible sequence-level redundancy are not automatically removed in
the MVP. Instead, `warning` may include labels such as
`possible_higher_order_or_partial:TXF000001-TXF000004`; inspect
`family_similarity.tsv` before collapsing or excluding a monomer.
Elastic families retain `uncalibrated_confidence` in their warning field.

## family_similarity.tsv

Produced by: `tandemx discover`

This file compares discovered representative monomers against each other. It is
a catalog-quality check that helps identify possible redundant monomers,
higher-order units, partial duplicates, or related families. Known repeats are
not used.

In default `--family-audit full`, all `F*(F-1)/2` pairs are emitted in catalogue order. A `.partial`
file during execution is incomplete; the final filename is atomically replaced
after the audit. Python and Rust use the same values and tie rules. This remains
an ungapped redundancy heuristic, including historical equality of ambiguous
characters; it must not be substituted for the gapped circular clustering score.
`--family-audit related` emits every non-distinct pair only, in the same order,
with identical values. Missing pair rows in that mode are distinct under the
current heuristic; they are not untested biological negatives.

`family_audit_summary.json` records `schema_version`, `complete`, `mode`,
`backend`, `family_count`, `possible_pairs`, `pairs_scored` (exact ungapped
comparisons), `pairs_pruned_by_kmer_gate` (provably distinct without alignment),
`emitted_pairs`, `related_pairs`, `omitted_distinct_pairs`, and `warning`.
Possible = scored + pruned; possible = emitted + omitted. Full mode has no
pruned/omitted pairs. The discovery completion receipt hashes this audit receipt.

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_a | string | NA | First family identifier |
| family_b | string | NA | Second family identifier |
| length_a_bp | integer | bp | Monomer length for `family_a` |
| length_b_bp | integer | bp | Monomer length for `family_b` |
| kmer_jaccard | float | fraction | Jaccard similarity of canonical monomer k-mer sets |
| shared_kmer_fraction | float | fraction | Shared canonical k-mers divided by the smaller k-mer set |
| local_identity | float | fraction | Best ungapped local identity over both orientations |
| local_overlap_bp | integer | bp | Aligned overlap length for the best local identity |
| local_overlap_fraction_shorter | float | fraction | Local overlap divided by the shorter monomer length |
| length_ratio | float | ratio | Longer monomer length divided by shorter monomer length |
| orientation | string | NA | `forward` or `reverse` orientation for the best local match |
| relationship | string | NA | `distinct`, `possible_higher_order_or_partial`, or `likely_redundant` |
| redundant_candidate | boolean | NA | Whether TandemX considers the pair a likely redundant representative |
| notes | string | NA | Interpretation notes for non-distinct pairs |

## collapsed_families.tsv and collapsed_monomers.fa

Produced by: `tandemx discover --collapse-redundant-families`

These optional files use the same schema as `families.tsv` and `monomers.fa`,
but include only families retained after collapsing pairs labelled
`likely_redundant` in `family_similarity.tsv`. They are not written by default.
`possible_higher_order_or_partial` families are not collapsed.

## family_collapse.tsv

Produced by: `tandemx discover --collapse-redundant-families`

This audit table records how the optional collapse mode treated each original
family. If no likely redundant family exists, rows are retained records and no
family is removed from the collapsed catalog.

| Field | Type | Description |
|---|---|---|
| original_family_id | string | Family identifier in the original de novo catalog |
| retained_family_id | string | Family retained in the collapsed catalog |
| action | string | `retained` or `collapsed` |
| reason | string | Decision reason |
| relationship | string | Relationship used for collapse; only `likely_redundant` should collapse |
| similarity_metrics | string | Semicolon-separated metrics supporting a collapse decision |
| notes | string | Interpretation notes |

## copy_number.tsv

Produced by: `tandemx quantify`

Input catalog: `monomers.fa` generated by `tandemx discover`.

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Repeat family identifier |
| monomer_length | integer | bp | Consensus monomer length |
| diagnostic_kmer_count | integer | k-mers | Number of selected diagnostic k-mers |
| median_kmer_depth | float | counts | Median diagnostic k-mer count in reads |
| haploid_depth | float | X | Haploid sequencing depth provided by the user or estimated as total read bases divided by genome size |
| estimated_copy_number | float | copies | Estimated monomer copy number |
| estimated_bp | float | bp | Estimated total repeat bp in reads |
| depth_mad | float | counts | Median absolute deviation of corrected diagnostic k-mer depths |
| copy_number_interval_low | float | copies | Empirical 10th-percentile copy-number estimate across diagnostic k-mers |
| copy_number_interval_high | float | copies | Empirical 90th-percentile copy-number estimate across diagnostic k-mers |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings, including `genome_background_uniqueness_not_verified` when diagnostic uniqueness was checked only within the catalog |

## repeat_density.bedgraph

Produced by: `tandemx locate`

bedGraph uses 0-based half-open intervals.

| Field | Type | Unit | Description |
|---|---|---:|---|
| chrom | string | NA | Chromosome or contig name |
| start | integer | bp | 0-based window start |
| end | integer | bp | 0-based half-open window end |
| value | float | fraction | Union repeat coverage divided by window size; always from 0 to 1 even when family arrays overlap |

Family-specific bedGraph tracks should use family-specific filenames or documented track metadata.

## arrays.bed

Produced by: `tandemx locate`

BED6 plus confidence fields, 0-based half-open. This file has no header line.

| Field | Type | Unit | Description |
|---|---|---:|---|
| chrom | string | NA | Chromosome or contig name |
| start | integer | bp | 0-based interval start |
| end | integer | bp | 0-based half-open interval end |
| family_id | string | NA | Repeat family identifier |
| score | integer | unitless | Exact diagnostic-k-mer support fraction scaled from 0 to 1000 |
| strand | string | NA | `+`, `-` or `.` |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated limitations such as `exact_kmer_identity_proxy` and `shared_family_kmers_excluded` |

## assembly_vs_read_cn.tsv

Produced by: `tandemx compare`. Also produced by `tandemx locate` for backward compatibility.

The formal compare MVP input is `copy_number.tsv` plus `arrays.bed`. `repeat_density.bedgraph` does not include `family_id`, so it is not suitable as the primary family-level comparison input.

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Repeat family identifier |
| read_estimated_bp | float | bp | Repeat bp estimated from reads |
| assembly_estimated_bp | float | bp | Repeat bp observed or inferred from assembly |
| assembly_read_ratio | float | ratio | Assembly estimate divided by read estimate |
| status | string | NA | `consistent`, `possible_collapse`, `possible_overexpansion`, `assembly_only`, `reads_only` or `low_confidence` |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

## probes.fa

Produced by: `tandemx probe`

Header format:

```text
>probe_id=TXP000001;family_id=TXF000001;length_bp=120;probe_score=0.87;confidence=medium
ACGT...
```

| Header Field | Type | Unit | Description |
|---|---|---:|---|
| probe_id | string | NA | Probe identifier |
| family_id | string | NA | Target repeat family |
| length_bp | integer | bp | Probe length |
| probe_score | float | unitless | Ranking score |
| confidence | string | NA | Confidence label |

## probes.rank.tsv

Produced by: `tandemx probe`

| Field | Type | Unit | Description |
|---|---|---:|---|
| probe_id | string | NA | Probe identifier |
| family_id | string | NA | Target repeat family |
| sequence_length | integer | bp | Probe sequence length |
| gc_content | float | fraction | GC fraction from 0 to 1 |
| tm | float | degrees C | Salt/formamide-adjusted long-oligo melting-temperature approximation |
| estimated_copy_number | float | copies | Target family read-based copy number |
| arrayiness_score | float | unitless | Fraction of predicted probe hits overlapping target arrays |
| specificity_score | float | unitless | Higher values indicate fewer predicted off-target hits |
| off_target_hits | integer | hits | Predicted off-target count |
| predicted_regions | string | NA | Semicolon-separated predicted target regions |
| probe_score | float | unitless | Combined ranking score |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings; the MVP includes `heuristic_probe_score_not_experimentally_calibrated` |

## in_silico_fish.tsv

Produced by: `tandemx probe`

| Field | Type | Unit | Description |
|---|---|---:|---|
| probe_id | string | NA | Probe identifier |
| chrom | string | NA | Chromosome or contig |
| start | integer | bp | 0-based predicted signal-region start |
| end | integer | bp | 0-based half-open predicted signal-region end |
| predicted_signal | float | unitless | Predicted signal strength score |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings, normally including `heuristic_signal_prediction` |

## Output Validation

Run:

```bash
tandemx validate --project results
```

The validator scans the project directory for recognized TandemX output filenames and checks:

1. required TSV fields;
2. numeric fields that must parse as integers or floats;
3. required `confidence`, `status`, or `warning` fields where defined by the schema;
4. BED and bedGraph 0-based half-open coordinates;
5. BED scores and strand values;
6. TandemX FASTA header structure for `monomers.fa` and `probes.fa`;
7. non-empty recognized output files, except pairwise/audit tables that can legitimately have no data rows when there are no pairs or no collapse events.

Currently recognized files are `candidate_reads.tsv`, `families.tsv`, `family_similarity.tsv`, `collapsed_families.tsv`, `family_collapse.tsv`, `repeat_annotation.tsv`, `copy_number.tsv`, `repeat_density.bedgraph`, `arrays.bed`, `assembly_vs_read_cn.tsv`, `output_manifest.tsv`, `probes.rank.tsv`, `in_silico_fish.tsv`, `monomers.fa`, `collapsed_monomers.fa`, and `probes.fa`.

## Pipeline Summaries

Produced by: `tandemx run` and `benchmarks/scripts/run_pipeline_benchmark.py`.

`pipeline_summary.tsv` has one row per requested step. `pipeline_summary.json` contains the same rows as a JSON array. Skipped and failed steps remain explicit rows.

| Field | Type | Description |
|---|---|---|
| run_id | string | Unique pipeline invocation identifier |
| input_reads | path list | Reads supplied to the pipeline, separated by semicolons when multiple files were provided |
| input_assembly | path or empty | Optional assembly supplied to the pipeline |
| max_reads | integer or empty | Configured read limit |
| max_read_bases | integer or empty | Configured cumulative read-base limit |
| kmer_backend | string | Selected `python` or `rust` backend |
| step | string | Pipeline step name |
| command | string | Shell-escaped executed command; empty for skipped steps |
| start_time | ISO-8601 string | UTC step start time |
| end_time | ISO-8601 string | UTC step end time |
| runtime_seconds | float | Measured wall-clock seconds |
| exit_status | integer | Process exit status; zero for successful or intentionally skipped steps |
| output_dir | path | Step output directory |
| output_validated | boolean | Whether expected outputs passed current validation |
| notes | string | Skip reason, failure reason, profiling state, or recorded thread setting |

Each successfully completed non-validation step also stores a hidden
`.resume_fingerprint.json` beside its outputs. `tandemx run --resume` compares
the effective command and SHA-256 hashes of all direct inputs before skipping
that step; changed or missing fingerprints cause the step to rerun.

## output_manifest.tsv

Produced by: every `tandemx run` invocation after pipeline setup succeeds.

The manifest inventories requested step outputs and pipeline-level reports. Missing outputs remain explicit rows. If an assembly-dependent step was requested without an assembly, each expected output row records `skipped_missing_assembly`.

| Field | Type | Description |
|---|---|---|
| step | string | Producing step or `pipeline` |
| output_type | string | Stable semantic output label |
| file_path | path | Output path under the configured run directory |
| exists | boolean | Whether the file existed when the manifest was finalized |
| file_size_bytes | integer | File size, or zero when absent |
| description | string | Human-readable output purpose |
| required_for_next_step | string | Downstream consumer or review purpose |
| notes | string | Step status, skip reason, or missing-output warning |

## run_report.md

Produced by: every `tandemx run` invocation after pipeline setup succeeds.

The Markdown report summarizes inputs, requested and completed steps, timings, output row counts, validation status, warnings, skips, primary paths, and suggested validation or post hoc matching commands. It is a run overview, not a biological claim report.

## known_repeat_matches.tsv

Produced by: `benchmarks/scripts/check_known_repeats_against_catalog.py`.

Known sequences are read only after discovery and are never inputs to `tandemx discover`.

| Field | Type | Description |
|---|---|---|
| known_repeat_id | string | Known repeat FASTA identifier |
| known_repeat_length | integer | Known repeat sequence length in bp |
| best_family_id | string | Best discovered catalog family |
| best_monomer_length | integer | Best monomer length in bp |
| similarity_score | float | Orientation-aware k-mer Dice similarity from 0 to 1 |
| shared_kmer_fraction | float | Fraction of known-repeat k-mers shared with the best monomer |
| orientation | string | `forward` or `reverse` orientation of the best match |
| interpretation | string | Calibrated MVP label: strong, possible, weak, or no k-mer match |

## repeat_annotation.tsv

Produced by: `tandemx annotate-repeats`.

This is a post hoc interpretation file. Known-repeat libraries are compared only
after `tandemx discover` has produced `monomers.fa`; they are not discovery
templates. The table reports the best known-repeat match for each discovered
family.

| Field | Type | Description |
|---|---|---|
| family_id | string | Discovered family identifier from `monomers.fa` |
| monomer_length | integer | Discovered monomer length in bp |
| best_known_id | string | Best matching known-repeat identifier |
| best_known_length | integer | Known-repeat length in bp |
| best_orientation | string | `forward` or `reverse` orientation of the best known-repeat comparison |
| shared_kmer_fraction | float | Shared k-mers divided by the smaller k-mer set |
| jaccard | float | K-mer Jaccard similarity |
| dice | float | K-mer Dice similarity |
| containment_discovered_in_known | float | Fraction of discovered monomer k-mers found in the known repeat |
| containment_known_in_discovered | float | Fraction of known-repeat k-mers found in the discovered monomer |
| local_identity | float | Best ungapped local identity |
| local_overlap_bp | integer | Local overlap length supporting `local_identity` |
| annotation_status | string | `strong_known_match`, `weak_known_match`, `no_known_match`, `possible_partial_match`, or `possible_higher_order_match` |
| notes | string | Interpretation note; does not convert annotation into discovery evidence |

## compare_runs.tsv

Produced by: `benchmarks/scripts/compare_tandemx_runs.py`.

This is a post hoc run-consistency report. It compares two TandemX discover or
pipeline output directories so users can distinguish parameter-driven catalog
differences from apparent result conflicts. It does not rerun discovery and does
not use external truth sequences.

The companion `compare_runs.md` is a human-readable summary of the same checks.

| Field | Type | Description |
|---|---|---|
| category | string | Check group: `input`, `metadata`, `discover_parameter`, `result`, or `interpretation` |
| item | string | Specific checked value, such as `reads`, `max_reads`, `threads`, `min_support_reads`, or `family_count` |
| run_a_value | string | Value observed in the first run |
| run_b_value | string | Value observed in the second run |
| same | boolean | Whether the normalized values match |
| direct_comparison_impact | string | `blocking_if_different`, `result_difference`, `informational`, or `summary` |
| notes | string | Reason a difference matters, including why direct comparison is not valid |

Discovery thread count is treated as an execution parameter. It may affect runtime
but should not change deterministic candidate or family output for the same
backend and biological parameters.

## Toy Simulation Outputs

Produced by: `tandemx simulate toy`

The toy simulation command uses only simulated sequences. It does not use real biological sequences and does not implement repeat discovery.

Files with the `truth_` prefix are simulator metadata for tests and benchmark assertions. They are not inputs to the real analysis workflow and are not required by `tandemx discover`, `quantify`, `locate`, `probe`, or `visualize`.

### reads.fa

Simulated FASTA reads sampled from the toy source sequence.

Header format:

```text
>toy_read_0001;source_start=1000;strand=+;error_rate=0.01
ACGT...
```

| Header Field | Type | Unit | Description |
|---|---|---:|---|
| toy_read_id | string | NA | Read identifier before the first semicolon |
| source_start | integer | bp | 0-based source-sequence start used for simulation |
| strand | string | NA | `+` or `-`; `-` reads are reverse-complemented before errors |
| error_rate | float | fraction | Per-base substitution error rate used for reads |

### assembly.fa

Simulated toy assembly FASTA.

| Header | Type | Unit | Description |
|---|---|---:|---|
| toy_chr1 | string | NA | Single simulated assembly contig |

### truth_monomers.fa

Truth monomer FASTA used by the simulator.

Header format:

```text
>family_id=TXF000001;monomer_id=TXM000001;length_bp=566;source=simulated
ACGT...
```

| Header Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Simulated repeat family identifier |
| monomer_id | string | NA | Simulated monomer identifier |
| length_bp | integer | bp | Monomer length |
| source | string | NA | Always `simulated` for toy data |

### truth_arrays.bed

Truth assembly arrays. Coordinates are BED-style 0-based half-open. This file has no header line.

| Field | Type | Unit | Description |
|---|---|---:|---|
| chrom | string | NA | Assembly contig name |
| start | integer | bp | 0-based array start |
| end | integer | bp | 0-based half-open array end |
| name | string | NA | Repeat family identifier |
| score | integer | unitless | BED score, 0 to 1000 |
| strand | string | NA | `+`, `-` or `.` |
| monomer_length_bp | integer | bp | Truth monomer length |
| assembly_copies | integer | copies | Number of monomer copies represented in the assembly |
| array_status | string | NA | `compressed` or `truth_like` |

### truth_copy_number.tsv

Truth read and assembly copy-number table.

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Simulated repeat family identifier |
| monomer_id | string | NA | Simulated monomer identifier |
| monomer_length_bp | integer | bp | Truth monomer length |
| read_copies | integer | copies | Number of monomer copies in the read-generating source sequence |
| assembly_copies | integer | copies | Number of monomer copies represented in the assembly |
| read_repeat_bp | integer | bp | `monomer_length_bp * read_copies` |
| assembly_repeat_bp | integer | bp | `monomer_length_bp * assembly_copies` |
| assembly_status | string | NA | `simulated_under_assembly` or `truth_like` |

### simulation_config.yaml

Stable YAML record of the simulation parameters and simulated family metadata. It intentionally excludes timestamps so that fixed-seed outputs are byte-for-byte reproducible.

Top-level fields:

| Field | Type | Unit | Description |
|---|---|---:|---|
| command | string | NA | `tandemx simulate toy` |
| seed | integer | NA | Random seed |
| num_reads | integer | reads | Number of simulated reads |
| read_length | integer | bp | Simulated read length |
| background_length | integer | bp | Random background length |
| error_rate | float | fraction | Per-base substitution error rate |
| monomer_lengths | list[integer] | bp | Simulated monomer lengths |
| copies | list[integer] | copies | Read-truth copy counts |
| families | list[object] | NA | Family IDs, monomer IDs, monomer lengths and assembly copies |

## Static Visualization Outputs

Produced by: `tandemx visualize`

The MVP writes SVG and PDF files. SVG text should remain editable where matplotlib supports it.

| File | Format | Description |
|---|---|---|
| catalogue_summary.svg | SVG | Family-level read bp, assembly bp and probe-score summary |
| catalogue_summary.pdf | PDF | PDF version of catalogue summary |
| assembly_vs_read.svg | SVG | Scatter plot of read-estimated bp versus assembly-estimated bp |
| assembly_vs_read.pdf | PDF | PDF version of assembly-vs-read plot |
| in_silico_fish.svg | SVG | Toy ideogram-like predicted signal plot |
| in_silico_fish.pdf | PDF | PDF version of in silico FISH plot |
# Research/release additions (2026-09-06)

All new benchmark TSVs have a header. `NA` denotes an unavailable/undefined
measurement (including zero denominators); JSON uses `null`. Coordinates below
are 0-based and half-open, including normalized TRF/TideHunter output.

## Discovery completion

`discovery_summary.json`: `schema_version`, `status` (`completed` or
`no_families`), `processed_reads`, `processed_bases`, `candidate_count`,
`family_count`, `warning`, and `output_sha256` keyed by output filename.
It is written only after successful analysis. Zero-result `families.tsv` and
`candidate_reads.tsv` retain their original headers; `monomers.fa` is empty.
They validate only when receipt counts and hashes permit the corresponding
empty output. `run_config.yaml` status is `discover_no_families` for this case.
Pipeline summary `notes=skipped_no_discovered_families` means dependent analysis
was not performed, with `output_validated=false`.

## Challenge datasets and predictions

- `truth_reads.tsv`: `read_id` (neutral identifier), `length_bp` (observed
  post-error length), `truth_positive` (0/1 planted array), `strand` (+/-),
  `negative_kind` (`NA`, `random`, `at_rich`, `dispersed`, `low_complexity`).
- `truth_arrays.tsv` and normalized `predictions.tsv`: `read_id`, `start`,
  `end`, `period` (bp), `sequence` (consensus, empty if not in source format),
  `family_id` (truth family or empty for unassigned predictions).
- `manifest.json`: scenario parameters, seed, actual reads/bases/array counts,
  generator description, warning, per-file sizes and SHA-256 hashes.
- `truth_monomers.fa`: `>truth_fN` identifiers; strictly evaluation input only.
- `matches.tsv`: `prediction_index` (0-based normalized prediction index),
  `read_id`, `start`, `end`, `period`, `matched_truth_index` (0-based or NA),
  `truth_family_id`, `iou`, `status` (`matched`/`unmatched`). An empty table may
  contain only its `prediction_index,status` header.
- `family_recovery.tsv`: `family_id`, `best_identity` (best circular identity),
  `recovered` (0/1 one-to-one assignment), `assigned_sequence_index` (0-based in
  sorted distinct canonical consensuses or NA), `criterion`. With no truth
  families its empty header is `family_id,recovered,criterion`.

## Challenge measurements

`raw_runs.tsv` has one row per scenario x seed x tool x repetition:

- Identification: `scenario`, `dataset_id`, `seed`, `split`, `tool`,
  `repetition`, `total_bases`.
- Execution: `exit_code`, `runtime_seconds`, `peak_rss_mib` (MiB, direct child
  only), `timed_out`, `status` (`ok`, `failed`, `invalid_output`), `error`.
- Arrays: `truth_array_count`, `predicted_array_count`, `matched_array_count`,
  `array_recall` = matched/truth, `array_precision` = matched/predicted,
  `array_f1` = 2*matched/(truth+predicted).
- Reads: `read_detection_recall` = called positive reads/positive reads,
  `read_detection_precision` = called positive reads/called reads,
  `negative_read_count`, `false_positive_read_count`,
  `negative_read_call_rate` = called negative reads/negative reads. These are
  planted-truth endpoints; incidental natural patterns in negative controls
  are not independently curated biological false positives.
- Errors: `matched_period_mae_bp` and `matched_boundary_mae_bp` (mean of start
  and end errors), both conditional on successfully matched arrays.
- Descriptive Wilson 95% intervals: `array_recall_ci_low/high`,
  `array_precision_ci_low/high`, `negative_read_call_rate_ci_low/high`.
  Arrays within a read and reads within a family are correlated; these are
  not confidence intervals across independent genomes or simulation seeds.
- Families: `truth_family_count`, `recovered_family_count`,
  `sequence_family_recall` = recovered/truth, `family_criterion`.
- Audit: `warning`, `prediction_sha256`, `catalog_sha256`.

`summary.tsv` groups scenario/dataset/tool with `seed`, `split`,
`successful_runs`, `attempted_runs`, `deterministic`, `median_runtime_seconds`,
`median_peak_rss_mib` and the above unqualified accuracy/error fields. Accuracy
is retained only when all attempted runs succeed and their predictions agree.
A single run has `deterministic=not_tested_single_run`. Failure rows never enter
successful-run accuracy aggregates.

`run_config.yaml` records the full effective suite and selected split/scenarios;
`environment.json` records source and executable hashes, Git HEAD, Python,
platform and measurement method. Each execution has `command.json`, logs and
`receipt.json`; `validation.json` records completion and success/failure counts.

New runs also save `source_snapshot/` with the actual Python/Rust source, native
extensions and build metadata, verifying each copy against its source hash.
TandemX subprocesses load this snapshot through PYTHONPATH and use it as their
working directory, preventing a checkout in the current directory from taking
import precedence. Early development
v1/v2 and validation runs recorded hashes only; their intermediate dirty source
is not a complete archived release. Final publication runs require a committed
revision and the full snapshot. External tools retain executable hashes and
separate pinned-source/build provenance.

The suite optionally records `discovery_method` (`legacy` or `elastic`) and
`ultra_options` (`window_size`, `windows`, `tune`, `tune_indel`). ULTRA rows use
`tool=ultra`; buffer/tuning configurations are separate suites with separate
receipts. Its normalized `predictions.tsv` uses the same half-open ArrayRecord
schema. Unresolved consensus symbols are N; missing consensus is an empty string.

## Figure source data and ENA subsets

`heatmap_source.tsv`: `panel`, `scenario`, `tool`, `seed`, `metric`, `value`,
`source_dataset`. `interval_source.tsv`: `track`, `read_id`, `start`, `end`,
`period`, `source` (source-table path; newer multi-run figures use absolute paths). `figure_provenance.json` ties
the figure to the input table hashes and documents evidence limits.

`ena_metadata.tsv` preserves ENA's requested columns: `run_accession`,
`sample_accession`, `scientific_name`, `instrument_model`, `read_count`,
`base_count`, `fastq_ftp`, `fastq_md5`, `fastq_bytes`. These describe the full
remote archive, not the extracted subset. `subset_receipt.json` separately
records the accession, retrieval UTC, URLs, sampling method, requested/observed
read counts and bases, downloaded compressed bytes and cap, output size/SHA-256,
first headers, and `source_full_file_md5_verified=false` for prefix extraction.
No complete-library checksum claim is made. `ena_run.xml` and `ena_sample.xml`
preserve retrieved public metadata. Partial extraction failures retain a
`.partial` file and never emit a successful subset receipt.

## Candidate sequence and monomer assignment evidence

`discover/candidate_monomers.fa` is emitted by both discovery methods. Headers:
`>candidate_id=TXC000001;read_id=<percent-escaped-read-ID>;length_bp=<length>`.
Sequences are the original read-local consensuses, prior to between-read
clustering. Each ID joins exactly one `candidate_reads.tsv` row. Read IDs are
percent-encoded to preserve delimiters. Empty output requires a verified
zero-candidate completion receipt.

`discover/monomer_membership.tsv` is emitted only by sequence clustering:

| Field | Definition |
| --- | --- |
| read_id, candidate_id | Input read and read-local candidate identifiers |
| cluster_id | Internal operational cluster, TXG-prefixed; NA for unresolved sequence |
| family_id | Retained output TXF cluster, or NA when not retained |
| representative_sha256 | SHA-256 of canonical representative sequence bytes; NA if unresolved |
| edit_distance_upper_bound | Integer edit count of a valid witnessed comparison; not necessarily minimum |
| similarity_lower_bound | 1 minus witnessed edits divided by the larger length, in [0,1] |
| minimum_cluster_identity | Configured operational similarity threshold, in (0,1] |
| compatible_cluster_count | Number of eligible existing clusters at assignment, or 1 for a new cluster, 0 if unresolved |
| alternative_cluster_ids | Semicolon-separated additional eligible clusters, or empty |
| status | assigned, below_minimum_support, or unresolved_sequence |
| warning | multiple_compatible_clusters, ambiguous_bases_not_clustered, or empty |

For unresolved candidates, edit upper bound is candidate length and similarity
lower bound is zero, with no comparison claimed. Internal cluster IDs and output
family IDs have separate rankings. `families.tsv` warnings identify operational
clusters, the similarity threshold, observed representatives and uncalibrated
confidence. `run_config.yaml` records requested and resolved clustering methods.
The completion receipt hashes both new files when present. A header-only
membership file is valid only with a verified zero-candidate receipt.

## Added independent challenge endpoints

All rates below are fractions (0–1); undefined denominators remain NaN in TSV,
null in JSON. Failures remain unavailable and are never replaced with zeros.

- `base_union_recall`, `base_union_precision`, `base_union_f1`: shared base pairs
  divided by the union of truth, union of predictions, or their harmonic mean.
  Intervals are merged per read first; this endpoint ignores repeat period.
- `predicted_union_bp`, `truth_union_bp`: per-read union lengths summed across
  reads. `duplicated_prediction_bp`: sum of raw predicted interval lengths minus
  predicted union. `duplicate_bp_fraction`: duplicated/raw predicted bases.
- `cyclic_monomer_recall`: maximum one-to-one matching at >=0.90 cyclic global
  unit-cost edit similarity to planted monomers. Every rotation and both strands
  are searched by independent edlib; N never matches N. A single consensus
  cannot recover two distinct planted monomers. This is a sequence recovery
  endpoint, not independently inferred family taxonomy.
- `distinct_consensus_count`: exact canonical distinct prediction sequences.
  `homologous_consensus_fraction`: proportion matching any planted monomer;
  multiple variants can match one truth, so this is not family precision.
  `unmatched_distinct_consensus_count` counts the remainder.
- `mean_best_cyclic_edit_similarity`: mean best similarity per planted monomer,
  irrespective of one-to-one assignment.
- `cyclic_monomer_recovery.tsv`: truth_id, best_cyclic_edit_similarity, recovered
  (0/1), assigned_sequence_index (canonical sorted 0-based index or NA), criterion.
  With no truth the empty header is truth_id,recovered,criterion.

`benchmarks/scripts/rescore_challenge.py` emits `rescored_metrics.tsv` plus
`provenance.json` and a source snapshot. It re-evaluates first-repetition archived
outputs, retaining source_run, scenario, seed, tool and original_status; it does
not rerun tools or measure new runtime. Original failures remain NA. Input files,
original configuration, scoring script and source snapshots are hash recorded.
Challenge configuration also accepts clustering_method and cluster_identity for
TandemX-only ablations. Other tool commands do not receive these settings.

## SRF workflow comparison outputs

See [srf_workflow.md](srf_workflow.md#measurements-and-files) for native field
semantics and all workflow receipt fields. The normalized predictions.tsv uses
the challenge ArrayRecord header. SRF workflow status is ok, no_catalogue,
no_eligible_kmers, or failed. The latter has unavailable accuracy; no_eligible_kmers
means the count stage completed with an observed empty dump, native graph assembly
was explicitly skipped, and no native catalogue is fabricated. Normalized zero
predictions in that state are conditional on the configured k/count filter.
Raw native BED/FASTA/abundance files retain upstream formats, including headerless
native tables; TandemX-normalized TSV files always have a header.

## Paired discovery performance experiments

`benchmarks/scripts/compare_discovery_snapshots.py` saves complete baseline/new
source snapshots, input hashes, commands, per-file output hashes, raw_runs.tsv,
summary.tsv and validation.json. Each dataset has randomized execution order and
at least two repetitions. It compares six discovery data files, excluding
run-specific timestamps/logs. Failure or any output difference fails the paired
gate. `speedup_baseline_over_native` is the ratio of median wall times;
`rss_ratio_native_over_baseline` is the ratio of median RSS. Ratios remain
measured diagnostics even if output equivalence fails, with the failed gate
clearly recorded. These are same-tool optimizations, not external-tool rankings.

The process measurement helper now also emits cpu_user_seconds and
cpu_system_seconds from wait4. They are absent in earlier records and unavailable
on non-wait4 platforms. raw_runs.tsv includes these values; paired summaries
report baseline/native medians. Sampling, compiler builds and profilers must not
run concurrently with the timing experiment.

## Conditional abundance experiment files

Generated by `python -m benchmarks.abundance.run`; see
[abundance_benchmark.md](abundance_benchmark.md) for the sampling and scoring
contract. Every normalized TSV has a header. Missing denominators are `NA`.

| File | Fields and interpretation |
| --- | --- |
| `genomes/*/truth_copy_number.tsv`, `assembly_*.truth.tsv` | `chrom`, `family_id`; `start/end` are 0-based half-open planted coordinates (equal for absent arrays); `period` in bp; integer `copies`; `repeat_bp=end-start` |
| `reads/*/*/sampling.tsv` | Neutral `read_id`, 0-based circular `genome_start`, pre-error `source_length` in bp, `strand`, count of observed `substitutions` |
| `copy_number_metrics.tsv` | Context `seed`, requested `coverage`, `substitution_rate`, `read_length`; `family_id`, `truth_copies`, `estimate`; `signed_relative_error=(estimate-truth)/truth`, its `absolute_relative_error`; `interval_low/high`, Boolean `interval_contains_truth`, `(high-low)/truth` as `interval_relative_width`; `sampling_oracle_copy_estimate` from sampled planted bases/period/actual depth, `oracle_relative_error`, `(estimate-oracle)/truth` as `estimator_minus_sampling_oracle`; explicit `interval_kind` |
| `localization_metrics.tsv` | `seed`, requested `assembly_fraction`, `family_id`; `true_assembly_bp`, union `predicted_assembly_bp`, `overlap_bp`; overlap/true as `base_recall`, overlap/prediction as `base_precision`, and merged interval count `fragments` |
| `comparison_metrics.tsv` | Read context above plus requested `assembly_fraction`, `family_id`; actual integer-copy `truth_assembly_read_ratio`, Boolean `truth_underrepresented`, native `predicted_assembly_read_ratio` and `native_status`, Boolean `predicted_underrepresented`, `outcome` TP/FN/FP/TN |
| `copy_number_summary.tsv` | Grouped by `coverage/substitution_rate`; `family_observations`; mean signed, mean/median absolute relative error; `interval_empirical_coverage` (fraction of truth-containing intervals, not nominal CI coverage), `mean_interval_relative_width`, `mean_oracle_relative_error`, `mean_estimator_minus_oracle` |
| `comparison_summary.tsv` | Grouped by `coverage/substitution_rate/assembly_fraction`; `family_observations`, TP/FN/FP/TN counts, `sensitivity=TP/(TP+FN)`, `false_positive_rate=FP/(FP+TN)`, `precision=TP/(TP+FP)` |

`environment.json` records detector/benchmark source hashes, Git revision,
config hash and platform. `source_snapshot/` retains executable source/native
extension; `run_config.json` is the exact experiment configuration. Dataset
manifests record file hashes and observed sampling/error totals. Each
`receipt.json` includes the complete command, exit code, timeout, wall/CPU/RSS.
Root `validation.json` counts completed executions and emitted score rows; its
execution success does not establish scientific acceptance. Original failed
receipts remain present and are not represented as zero accuracy.

Performance figure `panel_source.tsv`: `panel`, `dataset`, `metric`, `variant`,
`value`; panel C transforms the RSS ratio to percent change. All other values
are the corresponding per-dataset median or ratio in the archived summary.
`figure_provenance.json` records input/figure/script hashes and derived ranges.

## Complete-data QC and cohort screening

See [complete_data_qc.md](complete_data_qc.md) for commands and JSON receipts.
`length_histogram.tsv`: integer `length_bp`, `read_count`.
`base_quality_histogram.tsv`: integer `phred` (ASCII minus 33), `base_count`.
`joint_distribution.tsv`: lower bounds `length_bin_kb` (1-kb width),
`gc_bin_percent` (1% width, denominator includes N), `mean_quality_bin_phred`
(5-Phred width, read quality from mean reported error probabilities), and
`read_count`. Histograms include all validated reads/bases.

Cohort `screen_summary.tsv`: queried `accession`, `metadata_rows`,
`pacbio_wgs_run_rows`, `distinct_nonempty_biosample_ids`,
`missing_biosample_run_rows`, sum of `reported_fastq_bytes` over those runs,
archived `metadata_file`, `metadata_sha256`, and `warning`. PacBio+WGS alone
does not establish HiFi. BioSample/run counts do not equal independent
specimens. Original `screening.json` used the exploratory label `samples`,
counting empty strings as well; explicit summary fields supersede that label
without altering raw history.

Whole-file sampling produces `sample_NNN.fastq.gz` plus `.ids.tsv` columns
`read_id_hex` (lossless hexadecimal primary ID bytes), `hash128_hex` (32 hex
digits), `length_bp`. Numeric suffixes follow ascending inclusion probability.
Subset `.length_histogram.tsv` and `.joint_distribution.tsv` use the QC fields
above. The plan stores exact integer thresholds; the receipt defines probability,
expected/observed read and base totals, hashes, and nominal total-base coverage.
An empty random sample has zero counts and JSON null distribution statistics.

## Experimental read-cluster replay

See [read_cluster_quantification.md](read_cluster_quantification.md) for formulas
and boundaries. This does not modify the public `copy_number.tsv` schema.

- `estimates.tsv`: `seed`, `coverage`, `substitution_rate`, `family_id`,
  `diagnostic_kmer_count`, `estimated_copy_number`, `sampling_standard_error`,
  `sampling_interval_low/high`, `positive_reads`, `effective_positive_reads`,
  `interval_status`, `warning`. Missing numerical estimates/intervals are `NA`.
- `metrics.tsv`: context/family, `truth_copies`, `estimate`, signed and absolute
  relative error, `sampling_oracle_copy_estimate`, `estimator_minus_sampling_oracle`,
  Boolean `interval_available`, `interval_contains_truth` (NA if unavailable),
  `interval_relative_width`, `interval_status`.
- `inputs.tsv`: context, read/catalogue/truth/sampling SHA256, observed
  `read_count`, `total_bases`, and k-mer opportunity `total_exposure`.
- `summary.tsv`: coverage/error groups, `total_family_conditions`,
  `estimates_available`, `intervals_available/unavailable`, mean signed/absolute
  error, `mean_estimator_minus_oracle`, `conditional_interval_coverage` and
  `mean_available_interval_relative_width`. Conditional coverage never includes
  unavailable intervals; their explicit count must accompany the result.

Real discovery pilot files are defined in
[real_comparator_pilot.md](real_comparator_pilot.md). Its observed call/base fields
are not recall/precision. Complete-input figure `panel_source.tsv` contains
`panel`, `dataset`, `metric`, `x`, `value`, `read_count`: A-C are length/GC/quality
histogram or cumulative probabilities; D is observed/expected Gb with dataset
index x. Figure provenance records all input/output/source hashes. Reference
receipts store contig names, exact lengths and per-base/ambiguity counts; the
original NCBI assembly report defines chromosome naming equivalence.

Published region extraction fields, source cells and coordinate-inference rules
are defined in [published_mo17_regions.md](published_mo17_regions.md). These
records are mixed regions, not complete satellite-base truth.

## Streamed factorial simulation

See [factorial_scale_simulation.md](factorial_scale_simulation.md).
All genomic and read coordinates below are 0-based, half-open.

- `genome_index.json`: contig, length, sequence byte offset, line bases/bytes;
  valid only for the hash-verified single-contig fixed-width `genome.fa`.
- `truth_copy_number.tsv`: `chrom`, `family_id`, `start`, `end`, `period`,
  `copies`, `repeat_bp`, `requested_gc`, `founder_gc`,
  `unit_substitution_rate`, `observed_unit_substitutions`, `founder_sha256`.
  Copy count is planted ancestral-unit count despite same-length substitutions.
- `sampling.tsv`: `read_id`, `genome_start`, `source_length`, `observed_length`,
  `strand`. A read may wrap the circular source; lengths and starts precede errors.
- `truth_read_segments.tsv`: `read_id`, observed `start/end`, `period`,
  `family_id`, `sampled_source_repeat_bp`, `at_least_two_source_units`.
  Partial segments are retained; the Boolean is not a guarantee of detectability.
- Genome/read manifests record parameters, input/script/output hashes, realized
  biological/sequencing substitutions, insertions/deletions, sampled repeat
  bases and separate source/observed coverage. Error counts obey observed bases
  = source bases + insertions - deletions. Sequence only, no synthetic quality
  scores. The controller records completed condition IDs and refuses reserved
  seeds or requests above the declared worst-case sequence-base budget.

The conditional scorer `benchmarks.abundance.run_stream_quantify` writes:

- `copy_number_metrics.tsv`: `seed`, `condition_id`, requested source `coverage`,
  `error_model`, sequencing `substitution_rate`, `insertion_rate`, `deletion_rate`,
  biological `unit_substitution_rate`, `period`, `requested_gc`, `founder_gc`,
  `array_scope` (`factorial` or `megabase`), actual `source_coverage`,
  `observed_coverage`, `read_count`, `total_bases`, and the previously documented
  conditional copy-number score fields (`family_id`, planted copies, estimate,
  signed/absolute errors, k-mer spread bounds, oracle and spread interpretation).
- `copy_number_summary.tsv`: group keys `coverage`, `error_model`,
  `unit_substitution_rate`, `array_scope`, with `family_observations` and the
  conditional-abundance summary statistics above. Counts are correlated family
  conditions, not independent plants.
- `runs/sSEED/condition_ID/execution.json`: exact command, input hashes, exit,
  timeout and resource measurements; `output/` contains unmodified public CLI
  products. `environment.json` holds source/dataset provenance and interpretation
  limits; `validation.json` records completed/successful executions, independent
  simulated genomes and scored family conditions. Scientific acceptance is not
  inferred from exit status. `run.log` records each attempted execution.

## Experimental multi-k replay

See [multik_quantification.md](multik_quantification.md). `NA` represents missing
model estimates, never zero error. This is separate from public quantify output.

- `estimates.tsv`: `seed`, `condition_id`, `family_id`, `monomer_length`,
  `extrapolated_copy_number`, `uncorrected_mean_k21`, `log_slope_per_base`,
  `effective_word_loss_probability` (only for nonpositive slopes),
  `max_absolute_log_residual`, `leave_one_k_out_log_range`, `status`, `warning`.
  Loss combines divergence/error/model effects; it is not measured read accuracy.
- `per_k.tsv`: seed/condition/family plus `k`, `diagnostic_kmer_count`,
  `total_exposure`, `valid_windows`, `mean_corrected_count`,
  `uncorrected_copy_number`. Counts are corrected for founder word multiplicity;
  ambiguous-base opportunities remain in exposure.
- `metrics.tsv`: seed/condition/family, coverage/error/divergence/period/GC/scope
  strata, `method` (`native_median_k21`, `mean_k21_exposure`, `multik_loglinear`),
  `truth_copies`, `sampling_oracle_copy_estimate`, `estimate`,
  `signed_relative_error`, `absolute_relative_error`,
  `estimator_minus_sampling_oracle`, `status`, `warning`. No sampling CI is supplied.
- `summary.tsv`: method/coverage/error/divergence/scope groups, total
  `family_conditions`, `estimates_available`, mean signed/absolute error, median
  absolute error and mean estimator-minus-oracle; means use available estimates.
- `inputs.tsv`: seed/condition, reads/catalogue SHA-256, observed read/base counts.
  Environment/execution/validation receipts retain source and baseline hashes,
  fixed k values, exit/timeout/resources and complete pairing checks.

## Experimental joint-read multi-k calibration

See [joint_multik_uncertainty.md](joint_multik_uncertainty.md). All missing
numbers are `NA`; Boolean coverage is also `NA` if no interval exists.

- `metrics.tsv`: `seed`, `condition_id`, `family_id`, `estimated_copy_number`,
  `log_sampling_variance`, `sampling_interval_low`, `sampling_interval_high`,
  `minimum_effective_reads` (minimum observed (sum Y)²/sum Y² over k),
  `positive_reads_by_k` (comma-separated counts in fixed increasing k order),
  `fit_status`, `interval_status`, `warning`; strata `coverage`, `error_model`,
  `unit_substitution_rate`, `array_scope`, `period`; `truth_copies`,
  `truth_covered` (inclusive bounds), `relative_interval_width=(high-low)/truth`,
  and `signed_relative_error=(estimate-truth)/truth`.
- `per_k.tsv` and `inputs.tsv`: same point/exposure and hash/count fields as
  the multi-k replay above; inputs also retain collector ambiguity/model warnings.
- `summary.tsv`: groups by seed/coverage/error/divergence/scope; total
  `family_conditions`, `intervals_available`, `truth_covered`,
  `conditional_coverage_fraction=covered/available` (NA when none available),
  `available_and_covering_fraction=covered/all`, `mean_relative_interval_width`
  over available intervals, and a dependent-family `warning`.
- JSON environment/execution/validation receipts retain source and previous
  result hashes, exit/resources, fixed k/support guard, pairing completion and
  held-out exclusion. Execution success alone is not scientific acceptance.

## Reference mapping QC

The manuscript's `paper/tables/input_cohort.tsv` uses `reported_species`,
`reported_material`, `run_accession`, `read_count`, `total_bases`,
`median_read_length`, `read_n50`, `gc_fraction`, `n_fraction`,
`exact_duplicate_read_ids`, `raw_sha256`, `qc_receipt`, `qc_receipt_sha256`
and `warning`. Lengths/counts are exact file-derived values; species/material
labels are source metadata. This table is not a biological-validation outcome.

See [reference_mapping_qc.md](reference_mapping_qc.md). `alignments.paf` is
unchanged native0-based, half-open minimap2 output; `mapping.sqlite` is scratch.

- `read_mapping.tsv`: `read_id`, `length_bp`, `gc_bin_5pct`,
  `max_purine_pyrimidine_bin_5pct` (integer bin indices0–19),
  `any_alignment_span_bp`, `primary_span_bp`, `primary_mapq20_span_bp`
  (per-read query unions), `primary_rows`, `secondary_rows`,
  `organelle_primary_rows`; `primary_organelle_span_bp`,
  `primary_other_reference_span_bp` and `primary_compartment_overlap_bp`.
  The last is organelle span plus other span minus total primary span.
  Unmapped reads have zero spans/counts.
- `composition_mapping.tsv`: `dimension`, `bin_start_percent`,
  `bin_end_percent`; sums of `read_count`, `total_bases`, `any_mapped_reads`,
  `primary_mapped_reads`, `primary_mapq20_reads`, the three span fields above,
  `organelle_primary_reads` and the three compartment span fields above.
  Counts are per dimension, not independent data.
- `reference_mapping.tsv`: `contig`, `reference_length_bp`, `compartment`
  (`organelle` explicitly declared or `other_reference`); primary and MAPQ20
  `aligned_target_bases`, `reference_union_bp`, `aligned_base_depth` fields.
  Prefixes are `primary_` and `mapq20_`. Depth is the sum of CIGAR-aligned
  target M/= /X blocks divided by contig length; it is not unique-copy depth.
- `environment.json`: source/helper/input hashes, exact command, native version,
  declared organelle list and limits. `reference_qc.json` retains full reference
  sequence totals. `execution.json` retains exit/time/RSS/CPU/timeout.
- `validation.json`: `complete`, exact query totals/hash, native alignment
  rows/types, all-read summary counts, mapped fractions and output hashes.
  Failed validation records `error`; partial PAF is not an accuracy result.

## Factorial discovery evaluation

See [factorial_discovery_benchmark.md](factorial_discovery_benchmark.md).
`inputs.tsv` has seed/condition/coverage/error labels, generation/condition/read/
truth/catalogue SHA-256 fields, observed read/base counts, all/eligible segment
counts, excluded partial/out-of-scope read count, and all/observed-eligible family
counts. Coordinates are observed read coordinates, 0-based half-open.

`summary.tsv` retains tool, seed/condition, process exit/timeout/resources, scoring
`status`, `evaluation_seconds` and warning. Metric prefixes define populations:
`eligible_read_*` applies the truth-only read exclusion to the documented challenge
array/read/base endpoints; `all_input_*` keeps every planted observed base for
union/duplicate endpoints; `all_planted_*` and `observed_eligible_*` contain family
count, recovered count, cyclic monomer recall and distinct-consensus support.
`all_input_prediction_count` and `predictions_on_excluded_reads` expose filtering.
Missing failed-stage metrics are `NA`, not zero sensitivity.

Per-tool `normalized_arrays.tsv` uses the standard read_id/start/end/period/
sequence/family_id schema. `eligible_array_details.tsv` indexes the filtered
prediction/truth lists, with matched coordinates/family and IoU. Family-recovery
TSVs hold `truth_id`, `recovered`, `assigned_sequence_index`, `threshold` and
`criterion`; indices refer to sorted exact canonical distinct consensuses.
Below-threshold exact sequence identities are not calculated or reported.
`metrics.json` uses null for undefined values. Standard source/execution/validation
receipts and logs accompany every run.

### Isolated clustering replay

`replay_clustering_isolated.py` writes `environment.json` with exact input,
core/helper source and baseline hashes, frozen snapshots and the precision/resource
scope. `baseline.json` / `packed.json` contain complete sorted-key family and
membership payloads. Their `.receipt.json` files record candidate/family counts,
clustering seconds, output and clustering-source hashes. `validation.json` keeps
each command, child exit/timeout, wall/CPU/peak RSS, worker metrics and exact
payload-hash parity, plus `complete`/`error`. Stage seconds exclude input loading
and serialization; child wall/RSS includes them. Failed comparisons are retained.

### Disk-backed real benchmark evaluation

`run_real_comparators.py` defaults to `--evaluation-backend disk` and records the
backend, evaluator source hash and SQLite cache policy in `environment.json`.
`evaluation.sqlite` retains a `reads` table (unique `read_id`, `length` in bp),
`metadata(name,value)` with a JSON `ready` input receipt, and `arrays` with tool,
1-based native in-scope ordinal, read ID and 0-based half-open start/end. It is a
scratch/evaluation index, not a genome annotation or accuracy truth source.
Native call order and `normalized_arrays.tsv` schema are unchanged. Failed
conversion/normalization leaves `.partial` artifacts, never a completed TSV or
FASTA. Unknown read IDs, out-of-bounds calls and late malformed rows fail.

`replay_real_normalization.py` writes `replay_receipt.json`: input/source hashes,
preparation/total seconds, optional previous-run path/environment hash, and one
entry per tool with descriptive metrics, `byte_identical`, `metrics_identical`,
native/normalized hashes and normalization seconds. `complete=false` plus
`error` preserves failure; preparation-only runs have no normalized tool rows.
It is evaluator validation, not a new discovery or accuracy measurement.

`fetch_morex_reference.py` retains `source_metadata.html` and a `reference_plan.json`
with original DOI/license, exact file/landing/download URLs, official SHA-256,
byte budget and source hashes. `reference_receipt.json` records `complete`,
transfer bytes/hash/state, full FASTA QC, plan/script/parser hashes and a material
matching warning. Transfer/QC errors keep `complete=false` and `error`; an
incomplete `.partial` file is not a reference release.

### Curated repeat-query source records

`source_units.tsv` has one row per accession: `known_id`, exact `accession`,
`source_length`, 1-based inclusive `start1`/`end1` (NA if excluded), official
`sequence_md5`, `species`, deposited `material`, `role`, source `evidence` URL,
`boundary`, full-sequence `source_sha256`, original `embl_sha256`/`fasta_sha256`,
`included_monomer`, derived `monomer_length`/`monomer_sha256` (NA if excluded),
and `warning`. `known_monomers.fa` contains only included queries with accession
and extraction coordinates in each header. `curation_receipt.json` identifies
source/query counts, curator/retrieval/output hashes, completion and scope.
These coordinates concern the small deposited record, not the test assembly.

Native-index ablations retain the isolated replay schema and add `comparison`
to the environment and `index_backend` (`python`/`native`) to each worker receipt.
Both variants use the same source/native binary and native alignment. Historical
source comparisons use `source_default`; earlier storage replays predate this
field and labelled the current variant `packed`.

### TRASH comparator execution and coordinate controls

Native catalogue sensitivity tables: `catalogue_metrics.tsv` contains
`consensus_policy` (`primary`, `primary_and_secondary` for TRASH1),
`native_region_count`, `secondary_missing_regions`,
`excluded_sequence_length_count`, `in_scope_native_sequences`, cyclic family
recovery/count/fraction fields and an explicit scope warning. Sequence length
30–1000 bp determines eligibility independently of predicted period.
`catalogue_POLICY_recovery.tsv` retains truth ID, recovered Boolean, assigned
native sequence index, threshold and matching criterion for every true family.

The container runner retains `image_inspect.json`, `input_qc.json`, its exact
source, stdout/stderr, original `native/` products and `container_provenance/`
package/version/hash records. `execution.json` has the immutable image ID,
command, input/helper hashes, host-client resource record, container exit/OOM
state, Linux command resource record, cgroup memory peak in bytes, native product
hashes and explicit completion/error/cleanup state. Startup failure has no
inference success; GNU time resource parsing failures remain `resource_error`.
`native/linux_time.tsv` reports elapsed/user/system seconds, maximum RSS KiB,
and exit code. GNU time may also retain a nonzero-exit status line. Cgroup peak
includes cache and is not interchangeable with maximum RSS.
New executions also retain `cgroup_cpu_before.txt`/`cgroup_cpu_after.txt`.
`execution.json.cgroup_cpu` has all monotonic counter differences, aggregate
`usage_seconds`, `user_seconds`, `system_seconds` and measurement scope.
Counter units ending in `_usec` remain microseconds in `delta_counters`; other
counters retain their native count units. Missing/invalid counters are explicit
resource errors, never invented totals. Historical executions lack this field.

`coordinate_offsets.tsv` reports native1-based offsets-2..2, total monomer rows,
width discrepancies, forward/strand-adjusted exact sequence matches, and
out-of-reference extractions. The control-only audit preserves original tables
and writes source/input/output hashes in `audit_receipt.json`; it does not apply
automatic correction or calculate biological recall/precision.

`replay_discovery.py` writes a frozen environment plus `validation.json` with
command, child execution, profile flag and all seven expected/observed product
hashes and equality flags. Profile mode additionally creates `discovery.prof`.
Any changed input/baseline product, tool failure or output mismatch leaves
`complete=false` and a concrete error.

### Factorial assembly comparator evaluation

`score_trash_factorial.py` consumes a generated `genome/` directory and a complete
container run. It verifies all source/native hashes before and after evaluation.
`evaluation.json` records parameters, source/input hashes, seed, completion/errors
and the following TSV hashes. This development evaluator retains at most100,000
interval records and10 MB of truth/catalogue input; it never loads the genome.

- `region_metrics.tsv`: `coordinate_policy` is TRASH1 `window_grid` (native
  zero-based inclusive window converted to half-open) or `r_extraction` (native
  R one-based inclusive extraction, with its start=0 special case); TRASH2 uses
  `one_based`. `period_source` is either `native_peak` or `consensus_length`.
  Native periodicities can be fractional (TRASH1 emitted171.5 in the factorial
  run); retain these decimals for tolerance and error scoring, without rounding.
  Both alternatives use the previously documented array and cyclic-recovery
  metrics. `raw_region_count` and `excluded_scope_count` show the effect of the
  fixed30–1000-bp period and100-bp span scope. One chromosome is not a read-level
  negative-control experiment, so read-level metrics are omitted.
  TRASH2 additionally reports `unit_extent`, using its native repeat-to-array IDs
  and the outermost observed unit boundaries; it does not infer membership from
  overlap with truth. Native approximate regions and actual unit unions remain
  separate, and arrays without units cannot be silently dropped.
- Each policy directory contains `normalized_arrays.tsv`, `array_details.tsv`
  and `family_recovery.tsv`, with the existing independent evaluator schemas.
- `unit_coordinate_audit.tsv`: offsets−2..2 from native1-based coordinates,
  `monomer_rows`, `width_discrepancies`, `strand_adjusted_matches` and
  `outside_reference`; exact sequence checks use bounded fixed-FASTA seeks.
- `unit_coverage_metrics.tsv`: explicit `offset_from_native_1based` (TRASH1 0
  and−1, TRASH2 0), `status`, `warning`, `unit_count` and base-union metrics.
  Unit width is not an inferred array period; units never enter array/family
  recall. An invalid sensitivity policy is unavailable, never silently clipped
  or scored as zero accuracy. Native CSVs remain unchanged.

The generic NCBI reference downloader writes `reference_plan.json` with exact
GCA version/BioProject, official directory, archived metadata hashes, byte/MD5
expectations and explicit transfer budget. `reference_receipt.json` retains
complete/error state, transfer and streaming FASTA QC, plan/source hashes.
Optional FCS metadata absent from the official manifest is explicitly unavailable.
