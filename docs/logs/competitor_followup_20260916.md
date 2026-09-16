# Local comparator follow-up, 2026-09-16

This is a bounded **availability/interface smoke** on an internal-disk synthetic input, not a common-endpoint accuracy or resource benchmark. The native files and failed attempts are archived under [`benchmarks/evidence/competitor_followup_20260916`](../../benchmarks/evidence/competitor_followup_20260916/) with `SHA256SUMS.txt`. The shared 5,700 bp input and deterministic generator are in the [earlier evidence](../../benchmarks/evidence/competitor_local_smoke_20260916/inputs/cendetecthor_pipeline.fa) and [`generate_competitor_smoke_inputs.py`](../../benchmarks/scripts/generate_competitor_smoke_inputs.py). Input SHA-256 is `0acf542594ff333c69ffa5bcb4850c8b3a75fc4c0a828a4d6b6765a7215616c5`. No T7 write, large genome run, base conda install or global shell change occurred.

| Tool/task applicability | Fixed source and environment | Observed local outcome |
| --- | --- | --- |
| TideCluster: assembly array discovery and clustering; full TAREAN/KITE HOR stages separate | [official source](https://github.com/kavonrtep/TideCluster) `8131547e75f73a1e61298a89ff4e861fbd613b7e`, reports 1.21.3; project `.venv` and `.conda-clustering`; source-built TideHunter 1.4.3 | `tidehunter` and `clustering` exited 0, produced native GFF3; 1 toy array/`TRC_1`, 60 bp consensus. Designed A/B 30 bp monomers and AB HOR **not recovered**. Full pipeline untested. |
| SRF: k-mer graph core motif discovery; full read-first KMC/mapping workflow separate | [official source](https://github.com/lh3/srf) `e54ca8c8eccf6b1f19428b0f862f2c90575290a0`; local Apple Clang `make` binary SHA-256 `1b78ea8669eb8a03c9be8a1848d354b156876a9ce3875997b617f07426d2c2fb` | Core `srf` exited 0 on deterministic canonical k=151 counts, emitted one 60 bp motif; KMC, read mapping and family scoring were not run. |
| TRASH1/2: assembly tandem repeat annotation; historical comparator exists | [TRASH1](https://github.com/vlothec/TRASH) `f7d53a013aa4664ec275f5f5544d3a3906a43c00`; [TRASH2](https://github.com/vlothec/TRASH_2) `f290a5efe2341b5525711dcaabfbfc9087edf3cb`; historical isolated Docker image `sha256:50ff2f2af62681892295d76ddd79b8c87a06135455ffe0e79f1c81dc07939243` | Prior 10 Mb synthetic run exited 0 with native outputs and input SHA-256 `15958d097a6e36fd49d7e8fc397dceb36297845b1e4b7fd4f06f0690ce15c57c` in [`TRASH_factorial_s6301`](../../paper/evidence/TRASH_factorial_s6301/). No new run: `docker info` exit 1, daemon socket unavailable. |
| HiCAT-human: human alpha-satellite read/assembly HOR annotation, **not a plant comparator** | [official source](https://github.com/xjtu-omics/HiCAT-human) `86558becd73e07a96845017b0e82fe68844ff8b7`, `VERSION` 1.0.0; pinned requirements SHA-256 `e8f86125d6994e3a0ae8963bbf2f11b13eb0579d14d902429e32780b5f51a1d6` | CLI `--help` exit 0 in `tandemx-dev`. Exact requirements `mamba create --dry-run --platform osx-arm64` exit 1 because `seqtk=1.2` is unavailable. No eligible matched human HiFi/assembly sample was enrolled, so no inference was attempted. |

## TideCluster run and retained failures

The source `TideCluster.py` SHA-256 is `6199d6d1e1066fd56c62e61b6897ab32b372dee129065a66eccb2a523ff44f29`; `tc_utils.py` is `0b796918acfc2b39c4d495c79b16d410ecd6ffcb79e68a4e7af03ca87183e1b8`. The upstream 1.21.3 conda package could not solve natively on osx-arm64 because it requires TideHunter 1.4.3. The small osx-64 conda TideHunter 1.4.3 package did install into a project prefix and reported 1.4.3, but execution under Rosetta died with `SIGILL` (archived `tidecluster_stderr.txt`). An ARM source build of [TideHunter v1.4.3](https://github.com/yangao07/TideHunter/tree/v1.4.3), commit `d4a7feb2b9edfefb267a30db5d6f0685fb79ea32` with pinned `abPOA` submodule `4dc66a346a0f98dd4858e066a1128215100631e6`, failed on x86 `immintrin.h` (archived `tandemx_tidehunter_143_make_stderr.txt`). The succeeding temporary source-tree build was:

```bash
make clean_all
make -j2 sse2=1 CC='/usr/bin/clang -arch x86_64' CXX='/usr/bin/clang++ -arch x86_64'
```

The built binary was copied into project-local `benchmarks/competitor_envs/tidecluster/.local-bin/TideHunter` (SHA-256 `67a8010d948013cd1f3e151f7f918cda3a5209edc4e5486f28d8ac671360d6a6`); it reported 1.4.3 and ran under Rosetta. Project `.venv` used Python 3.11.15, NetworkX 3.4.2. Project `.conda-clustering` used MMseqs2 18.8cc5c, BLASTN 2.17.0+ and dustmasker 1.0.0. The dependency prefixes and output workdir are gitignored; exact native outputs are archived.

From an empty temporary workdir with `PATH=<project>/.conda-clustering/bin:<project>/.local-bin:$PATH`, the commands were:

```bash
<project>/.venv/bin/python <pinned-tidecluster>/TideCluster.py tidehunter \
  -c 1 -f <evidence>/inputs/cendetecthor_pipeline.fa -pr tcsmoke2 \
  -T '-p 30 -P 100 -c 2 -e 0.25' --max_memory 1
<project>/.venv/bin/python <pinned-tidecluster>/TideCluster.py clustering \
  -c 1 -f <evidence>/inputs/cendetecthor_pipeline.fa -pr tcsmoke2 -m 100
```

`<project>` resolves to the absolute `benchmarks/competitor_envs/tidecluster` directory, `<pinned-tidecluster>` to `/private/tmp/tandemx_tidecluster_20260916`, and `<evidence>` to `benchmarks/evidence/competitor_local_smoke_20260916`. Both exits were 0. The first GFF3 is SHA-256 `300a3b69db72c02040e0668e40b25ce99b0e8166d7f0b53c0e824c10a4440eb4`; clustered GFF3 `06c315bf5285e821e33c211bdebab466b0a90223146be67f9991fd3deb94cf72`; `TRC_1.fasta` `a778761f44ec5070274f497c4295cb5ff9522847f028f3998a9d3c0710ad3311`. The native line spans 1–5700 and names `TRC_1`; TideHunter's consensus length is 60. This confirms command compatibility and an assembly array call only. `-m 100` admits the tiny input and is not a pre-registered benchmark threshold.

## SRF core run and unresolved whole-workflow gate

The 5,700 bp fixture generates canonical k=151 count text using `min(kmer, reverse_complement(kmer))`, sorted by k-mer; this is a KMC-style **synthetic count fixture**, not output of KMC. Its SHA-256 is `c4130d557ae375dd47f035a2fb3dab2310b81795aeedfece0004687d8298f7fa` (210 rows, total 5,550 k-mer observations). With the pinned `srf` source built by `make` in a temporary directory:

```bash
srf -p srftoy -c 20 <evidence>/inputs/srf_k151_counts.txt > tandemx_srf_smoke_output.fa 2> tandemx_srf_smoke_stderr.txt
```

Exit was 0, stderr reported 60 distinct retained k-mers, and native FASTA SHA-256 `44797e2592b3ff2498e4fc78b49fe848b2e2448ec435333783befc0bf2360ed0` contains one 60 bp circular unit. `-c 20` is a smoke filter, not a tuned plant setting. Historical seven-stage SRF runs and source hashes remain in [`srf_workflow.md`](../srf_workflow.md); this new check does not replace them.

## Technical and eligibility boundaries

HiCAT-human's exact upstream requirements are Python-era pins (including `seqtk=1.2`, lastz, older scikit-learn); the archived `hicat_mamba_stdout.txt` preserves the native solver failure. A generic 30 bp plant-like toy cannot serve as an eligible human alpha-satellite sample. TRASH native results are historical and should not be re-labeled as a current rerun. For all four tools, output presence and process exit are recorded separately from correctness; none of these new smokes enrolls a donor-matched validation or held-out biological truth set.
