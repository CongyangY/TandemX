# TandemX File Formats

All TSV files must use tab separators, include a header line and use stable column names. Fields with uncertain interpretation should include `confidence` or `warning` where practical.

All BED files use 0-based half-open coordinates: `start` is included and `end` is excluded.

The current repository implements the toy simulator and the toy-scale `discover` MVP. The formats below define current and planned core outputs.

## candidate_reads.tsv

Produced by: `tandemx discover`

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

## copy_number.tsv

Produced by: `tandemx quantify`

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

BED6, 0-based half-open.

| Field | Type | Unit | Description |
|---|---|---:|---|
| chrom | string | NA | Chromosome or contig name |
| start | integer | bp | 0-based interval start |
| end | integer | bp | 0-based half-open interval end |
| name | string | NA | Repeat family or array identifier |
| score | integer | unitless | Scaled score from 0 to 1000 |
| strand | string | NA | `+`, `-` or `.` |

## assembly_vs_read_cn.tsv

Produced by: `tandemx compare`

| Field | Type | Unit | Description |
|---|---|---:|---|
| family_id | string | NA | Repeat family identifier |
| read_estimated_repeat_bp | float | bp | Repeat bp estimated from reads |
| assembly_estimated_repeat_bp | float | bp | Repeat bp observed or inferred from assembly |
| read_to_assembly_ratio | float | ratio | Read estimate divided by assembly estimate |
| status | string | NA | `consistent`, `possible_under_assembly`, `possible_over_expansion` or `uncertain` |
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
| rank | integer | rank | Rank, with 1 as the highest ranked probe |
| probe_length_bp | integer | bp | Probe length |
| gc_fraction | float | fraction | GC fraction from 0 to 1 |
| tm_estimate_c | float | degrees C | Simple melting temperature estimate |
| estimated_copy_number | float | copies | Target family read-based copy number |
| specificity_score | float | unitless | Higher values indicate stronger predicted specificity |
| off_target_count | integer | hits | Predicted off-target count |
| predicted_signal_region | string | NA | Region summary, or empty if unavailable |
| probe_score | float | unitless | Combined ranking score |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

## in_silico_fish.tsv

Produced by: `tandemx probe`

| Field | Type | Unit | Description |
|---|---|---:|---|
| probe_id | string | NA | Probe identifier |
| family_id | string | NA | Target repeat family |
| chrom | string | NA | Chromosome or contig |
| region_start | integer | bp | 0-based predicted signal-region start |
| region_end | integer | bp | 0-based half-open predicted signal-region end |
| signal_score | float | unitless | Predicted signal strength score |
| off_target_signal_score | float | unitless | Predicted off-target signal score |
| specificity_class | string | NA | `high`, `medium`, `low` or `uncertain` |
| confidence | string | NA | Confidence label |
| warning | string | NA | Semicolon-separated warnings or empty |

## Toy Simulation Outputs

Produced by: `tandemx simulate toy`

The toy simulation command uses only simulated sequences. It does not use real biological sequences and does not implement repeat discovery.

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
