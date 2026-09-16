# Methods benchmark input and split gate (2026-09-16)

This is a preregistration scaffold, not a new scientific benchmark or approval to run one. The machine-readable [input inventory](../benchmarks/inputs/methods_benchmark_manifest_v1.json), [protocol](../benchmarks/inputs/methods_benchmark_protocol_v1.json), and [validator](../benchmarks/scripts/validate_input_manifest.py) are source-side artifacts. They do not alter the frozen TandemX estimator. A valid manifest means its listed bytes and structure passed validation; each dataset's `status` still controls whether it can enter an analysis.

## Field contract

Every dataset records `dataset`, `sample`, `species`, `truth_type`, `family`, `array_coordinates`, numeric or null `coverage`, `read_source`, `assembly_source`, `split` (`development`, `validation`, `final-heldout`), boolean `allowed_tuning`, `metrics`, `competitors`, `software_versions`, `donor_id`, `accessions`, `assembly_lineage`, `family_homology_group`, `status`, and `files`. Each file has absolute `path`, `kind` (`fasta` or `fastq`), exact `bytes` and SHA-256. Empty family/coordinate lists mean no asserted repeat truth, not a negative biological result. `coverage=null` means unavailable; the YSD56 sample coverage is nominal against the recorded source size, not an independent read depth. Each run must add source/build/executable hashes, command/config, host/OS/compiler, threads, limits, output hashes, resource traces and explicit failure rows before formal scoring. File manifest rows never imply matched read/assembly provenance by themselves.

The validator streams file hashes and FASTA/FASTQ structure, including gzip decompression through trailer/CRC. It checks nonempty sequences, FASTQ sequence/quality length, required fields, allowed values, tuning only in development, and split leakage across species, donor, sample, accession, assembly lineage and family homology group. An `invalid` record with no files preserves a known failure without reopening damaged data. A `valid: true` report with `source-unresolved` rows is **not** scientific input eligibility.

## Partition and scoring boundary

Assign a whole species/donor/accession/assembly lineage and every derived subset, edited assembly and homologous family variant to one split. Use prior Ey15-2, Macadamia, YSD56 and public Col-CEN material only for development or positive controls. Validation is a new frozen-candidate check; if it informs a change, it becomes development history. Final held-out requires at least two previously unused plant species plus independent donor/family groups and one planned reveal after the endpoint, denominator, identifiability, primary baseline, effect threshold, power, resource budget and multiple-comparison rules are fixed. No independent result has been opened by this inventory. The protocol JSON specifies truth types, comparator tracks, endpoint units and run/failure states. A 2–4 donor feasibility study cannot substitute family counts for independent statistical replicates.

## Read-only source check

The branch at this check was `codex/srf-submission-comparison-20260910`, HEAD `00290a1`; worktree contained unrelated untracked manuscript/figure material, left untouched. `mount` showed `/dev/disk5s2 on /Volumes/T7 (exfat, ... fskit)`. `diskutil info /Volumes/T7` could not use the DiskManagement framework in this process, so current volume health was **not** established. The 2026-09-14 K30076 block-2 chunk-3 changed SHA despite an unchanged 370,177,772-byte length, failed gzip, and a later T7 write returned EIO and unmounted the volume. Its aggregate stays `invalid_no_acceptable_aggregate`; block 3 and the comparator were not started. We did no large T7 write, re-download, repair or retry.

The validator re-read only three bounded existing files from T7:

| Input | SHA-256 | Bytes | Structural readback |
| --- | --- | ---: | --- |
| ERR6210723 prefix 1000 `reads.fa` | `4adca195ae3587c52bf503751fd82fe68ae755fe15c7d10f3c09d2d1a1a2305f` | 15,710,688 | FASTA 1,000 records / 15,666,956 bp |
| Col-CEN v1.2 `Col-CEN_v1.2.fasta.gz` | `b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8` | 38,265,592 | gzip/FASTA 7 records / 132,081,078 bp |
| YSD56 SRR28726931 `sample_001.fastq.gz` | `b03acf91c941b72fc5b5b2e40a3646596944cc90ed272b1b9d9806970268ae82` | 6,635,734 | gzip/FASTQ 673 records / 11,669,565 bp |

The ERR prefix is the first 1,000 archive records, not an unbiased abundance sample; its source full-file checksum and gzip trailer were not verified by prefix extraction. Col-CEN is a separately sourced assembly with unresolved donor/extraction/haplotype matching to those reads. YSD56's small sample is a nested technical subset of SRR28726931, not a new plant. These files are `source-unresolved`, not an eligible read/assembly matched benchmark. Col-CEN includes organelles; select nuclear contigs explicitly in any future task. K30076 remains a failure record without re-reading the bad chunk.

Runtime inventory: `tandemx-dev` exists with Python 3.11.15 and editable TandemX 0.1.0; Apple clang 15.0.0 is available. Local `.conda-benchmark/bin/trf -v` returned 4.10.0-rc.2 and `.conda-benchmark/bin/TideHunter -v` returned 1.5.5. `rustc`/`cargo` and current SRF/CENdetectHOR/TideCluster executables were not located on this process's PATH; historical versions must not be treated as current installed versions. Freeze exact hashes and builds before any future same-host comparison.

Run the bounded gate with:

```bash
conda run -n tandemx-dev python benchmarks/scripts/validate_input_manifest.py benchmarks/inputs/methods_benchmark_manifest_v1.json
conda run -n tandemx-dev python -m pytest -q tests/unit/test_validate_input_manifest.py
```

On this check the manifest validator returned `valid: true`, three structurally valid small files and one `invalid` K30076 record; 9 focused pytest cases passed. Neither result establishes stable T7 write reliability, matched biological truth or benchmark performance. Before any large acquisition or compute, inspect the connection and perform a separately authorized bounded write/readback gate in a new T7 directory, then revalidate every intended input file and freeze the formal protocol.
