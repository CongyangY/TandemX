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

## candidate_reads.tsv

Produced by: `tandemx discover`

This is a de novo discovery output used as the repeat catalog for downstream commands.

| Field | Type | Unit | Description |
|---|---|---:|---|
| read_id | string | NA | Input read identifier |
| candidate_id | string | NA | Candidate repeat span identifier |
| read_start | integer | bp | 0-based start in read |
| read_end | integer | bp | 0-based half-open end in read |
| strand | string | NA | `+`, `-` or `.` |
| period_bp | integer | bp | Estimated repeat period |
| repeat_span_bp | integer | bp | Length of candidate repeat span |
| unit_count | float | copies | Approximate repeat unit count |
| score | float | unitless | Internal candidate score |
| low_complexity_flag | boolean | NA | Whether candidate is low complexity |
| confidence | string | NA | `high`, `medium` or `low` |
| warning | string | NA | Semicolon-separated warnings or empty |

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
| mean_identity | float | fraction | Mean candidate-to-consensus identity |
| low_complexity_flag | boolean | NA | Whether family is low complexity |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

Families with possible sequence-level redundancy are not automatically removed in
the MVP. Instead, `warning` may include labels such as
`possible_higher_order_or_partial:TXF000001-TXF000004`; inspect
`family_similarity.tsv` before collapsing or excluding a monomer.

## family_similarity.tsv

Produced by: `tandemx discover`

This file compares discovered representative monomers against each other. It is
a catalog-quality check that helps identify possible redundant monomers,
higher-order units, partial duplicates, or related families. Known repeats are
not used.

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
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

## repeat_density.bedgraph

Produced by: `tandemx locate`

bedGraph uses 0-based half-open intervals.

| Field | Type | Unit | Description |
|---|---|---:|---|
| chrom | string | NA | Chromosome or contig name |
| start | integer | bp | 0-based window start |
| end | integer | bp | 0-based half-open window end |
| value | float | fraction | Repeat density, typically covered bp divided by window size |

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
| score | integer | unitless | Scaled score from 0 to 1000 |
| strand | string | NA | `+`, `-` or `.` |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

## assembly_vs_read_cn.tsv

Produced by: `tandemx locate` in the current MVP. A dedicated `tandemx compare` command is planned.

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Repeat family identifier |
| read_estimated_bp | float | bp | Repeat bp estimated from reads |
| assembly_estimated_bp | float | bp | Repeat bp observed or inferred from assembly |
| assembly_read_ratio | float | ratio | Assembly estimate divided by read estimate |
| status | string | NA | `consistent`, `possible_collapse`, `possible_overexpansion` or `low_confidence` |
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
| tm | float | degrees C | Simple melting temperature estimate |
| estimated_copy_number | float | copies | Target family read-based copy number |
| arrayiness_score | float | unitless | Fraction of predicted probe hits overlapping target arrays |
| specificity_score | float | unitless | Higher values indicate fewer predicted off-target hits |
| off_target_hits | integer | hits | Predicted off-target count |
| predicted_regions | string | NA | Semicolon-separated predicted target regions |
| probe_score | float | unitless | Combined ranking score |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

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
| warning | string | NA | Semicolon-separated warnings or empty |

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

Currently recognized files are `candidate_reads.tsv`, `families.tsv`, `family_similarity.tsv`, `collapsed_families.tsv`, `family_collapse.tsv`, `repeat_annotation.tsv`, `copy_number.tsv`, `repeat_density.bedgraph`, `arrays.bed`, `assembly_vs_read_cn.tsv`, `probes.rank.tsv`, `in_silico_fish.tsv`, `monomers.fa`, `collapsed_monomers.fa`, and `probes.fa`.

## Pipeline Summaries

Produced by: `tandemx run` and `benchmarks/scripts/run_pipeline_benchmark.py`.

`pipeline_summary.tsv` has one row per requested step. `pipeline_summary.json` contains the same rows as a JSON array. Skipped and failed steps remain explicit rows.

| Field | Type | Description |
|---|---|---|
| run_id | string | Unique pipeline invocation identifier |
| input_reads | path | Reads supplied to the pipeline |
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
| item | string | Specific checked value, such as `reads`, `max_reads`, `min_support_reads`, or `family_count` |
| run_a_value | string | Value observed in the first run |
| run_b_value | string | Value observed in the second run |
| same | boolean | Whether the normalized values match |
| direct_comparison_impact | string | `blocking_if_different`, `result_difference`, `informational`, or `summary` |
| notes | string | Reason a difference matters, including why direct comparison is not valid |

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
