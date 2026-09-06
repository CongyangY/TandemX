# SRF workflow: scope, build and evidence

SRF is a required read-first satellite comparator. Its native output is a motif
catalogue that can include higher-order units. The comparison runner executes
seven native stages: KMC, kmc_dump, SRF, srfutils enlong, minimap2, srfutils paf2bed,
and srfutils bed2abun. No truth sequence or label enters these commands.

## Pinned tools for the macOS pilot

| Component | Source / version |
| --- | --- |
| SRF | lh3/srf, e54ca8c8eccf6b1f19428b0f862f2c90575290a0 |
| KMC | refresh-bio/KMC, 751ef36a3c1ccc6dda664f529ad218dc51d76f55, reports 3.2.4 |
| minimap2 | lh3/minimap2, 3c28777e7e2dcc90f825de1b9f17a89cca7d4452, reports 2.31-r1302 |
| k8 | attractivechaos/k8 release v1.2, native arm64-Darwin reports 1.2-r137 |

Build/log/hash evidence is in `paper/evidence/comparator_builds`. KMC initially
failed under Apple Clang because GNU ext/algorithm is unavailable in libc++.
`prepare_kmc_libcxx.py` replaces that header and GNU copy_n with the standard C++
equivalents; its complete diff is retained. Independent canonical counts agree
exactly at k=17/151 on reverse-complement, N-containing and multiline FASTA inputs.
The first incremental build includes a CLI object compiled as C++17 and other
objects as C++14. This is explicitly a development build; a clean consistent
rebuild is required before a final publication/release benchmark.

The actual successful KMC link uses `tandemx-dev/lib/libz.1.3.2.dylib` and a local
rpath. The unversioned libz.dylib path was absent. A failed nested zlib build and
its logs are retained. KMC help advertises a 1-GB minimum but its parser rejected
-m1 and required -m2. Memory values are measured RSS, not that allocation limit.
No system compiler or package was installed. SRF used Apple Clang -O3; minimap2
used its official aarch64=1 configuration. k8 came from its official v1.2 asset;
the observed SHA-256 is recorded, and no upstream checksum verification is claimed.

## Run and normalization

Activate `tandemx-dev`, place the pinned tools under a chosen tools/src directory,
then run a new output directory, for example:

```bash
python -m benchmarks.scripts.run_srf_pilot \
  --datasets /path/to/challenge/datasets/clean_171_s1101 \
  --tools-root /path/to/tools/src --minimum-counts 100 20 \
  --outdir /path/to/new-srf-run
```

This pilot tool layout currently expects macOS ARM64 k8. General platform/tool
path configuration remains required for the final release workflow. Use official
SRF/KMC/minimap2/k8 installation guidance for other platforms; do not relabel a
failed or unavailable dependency as zero accuracy.

KMC uses canonical k=151 counts, ci100 or ci20, counter saturation 100000, one
thread parameter, disk mode and a 2-GB strict allocation limit. ci100 is the
README's illustrative preset; ci20 is a sensitivity setting. Synthetic average
genome coverage is not known and no optimal tuning is claimed. Mapping and native
filtering use upstream defaults/recommended options. Abundance is divided by the
observed total FASTA bases, not an invented haploid genome size.

Native BED coordinates are zero-based half-open. The final interval set uses
positive upstream keep flags (1 or 2). Native outputs remain intact; common
challenge scoring retains periods 30–1000 bp and spans >=100 bp. Motifs/HORs beyond
that range are counted separately. No decomposition of native HORs is performed;
minimal-monomer recovery is not a complete assessment of SRF HOR discovery.

An empty successful KMC dump triggers SRF's native ca_kmer_read assertion. The
first 32-run attempt retains all 11 such crashes as failures. The guarded runner
records `no_eligible_kmers`, skips graph assembly and downstream commands, and
creates no fake native FASTA. Its empty normalized prediction set is explicitly
conditional on the observed count filter. An actual successful SRF invocation
with no emitted motifs has `no_catalogue`. Distinguish these two states from a
failed process. The wrapper guard is disclosed, not attributed to native SRF.

## Measurements and files

`workflow_receipt.json` includes input SHA-256, k, minimum_count, total_input_bases,
threads_parameter, status, individual stage exit codes/timings/RSS, skipped stages,
workflow_wall_seconds, external_stage_wall_seconds and temporary/output bytes
remaining after execution. Remaining bytes are not peak scratch-disk usage.
Maximum external stage RSS is max direct-child wait4 RSS across sequential native
stages; benchmark-controller memory is excluded. Tool input counting and command
orchestration are included in workflow wall time; independent truth scoring is
excluded. These pilot runs have no timing repetitions or final superiority claim.

Each stage saves its exact command and stderr. Native files are counts.kmc_pre,
counts.kmc_suf, counts.txt, srf.fa, srf.enlong.fa, srf.paf, srf.bed and
srf.abundance.tsv when the stage ran. Native files are intentionally unmodified
and retain upstream header conventions. Normalized predictions.tsv has the shared
ArrayRecord header. summary.tsv adds dataset, minimum_count, status, native motif
lengths/catalogue size, out-of-scope counts, scoring_scope, truth hashes and common
accuracy fields; failed rows are NA. environment.json snapshots source plus the
runner and records native executable hashes. validation.json records attempt and
failure counts. The four-run initial pilot's exact script was reconstructed and
verified against its original stored hash; it is archived with that evidence.

Native BED columns: read ID, start, end, motif ID, upstream mapping score, read
length, native motif length, keep flag. Native abundance columns: motif ID,
retained mapped bp, length-weighted native score, retained mapped bp/denominator,
all mapped bp/denominator. The native score is not a calibrated probability or
an independently measured sequence-identity fraction.

Primary sources: [SRF](https://github.com/lh3/srf),
[KMC](https://github.com/refresh-bio/KMC),
[minimap2](https://github.com/lh3/minimap2),
[k8](https://github.com/attractivechaos/k8).
