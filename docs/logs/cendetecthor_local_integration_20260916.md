# CENdetectHOR local integration and comparator status (2026-09-16)

This is a small synthetic **interface and availability smoke**, not a formal accuracy or equal-endpoint benchmark. All work used internal disk. The T7 volume and Docker daemon were not repaired or used. Native inputs, outputs and failed-run logs are in [`benchmarks/evidence/competitor_local_smoke_20260916`](../../benchmarks/evidence/competitor_local_smoke_20260916/); `SHA256SUMS.txt` covers each archived evidence file. The fixture generator is [`generate_competitor_smoke_inputs.py`](../../benchmarks/scripts/generate_competitor_smoke_inputs.py).

## Exact upstream and local software

| Component | Official source and pinned revision | Local state |
| --- | --- | --- |
| [CENdetectHOR pipeline](https://github.com/CENdetectHOR/CENdetectHOR_pipeline) | `969f4ef515fadd8b8a035357294971ff57d4618f` | Snakemake source run locally; `Snakefile` SHA-256 `7cd4a1ba66224751b37f316fb7f3a33e7a33c18c087407cb30ce126f524ebfe1` |
| [CENdetectHOR library](https://github.com/CENdetectHOR/CENdetectHOR_lib) | `bdc3a9df5b60adafdd17eaded1b994edd22aae03`, source version `0.0.6` | Installed from that checkout into project `.venv`; `pyproject.toml` SHA-256 `74e08449c2f7205d77cc290d27052ba1e4e7562f72346788dc5e0478f11a0952` |
| [StringDecomposer](https://github.com/ablab/stringdecomposer) | `0a967f6bf131face88397445ca49b65fee49489c` | Source compiled with `/opt/homebrew/opt/llvm/bin/clang++`; binary `dp` SHA-256 `f7e38092bfe7c8f10c36340f457bc5b301ef1b65bfa76a2d42342c3ef921c0a6` |

Project-isolated prefixes: `benchmarks/competitor_envs/cendetecthor/.venv`, `.conda-tools`, `.conda-r`. They are gitignored; nothing was installed into conda base or the global shell. Python 3.11.15, Snakemake 8.14.0, `cen_detect_hor` 0.0.6 source checkout, Biopython 1.88, pandas 3.0.5 and edlib 1.3.9.post1 were observed. The complete Python package list is [`cendetecthor_python_freeze.txt`](../../benchmarks/evidence/competitor_local_smoke_20260916/cendetecthor_python_freeze.txt). The tool prefixes also contain SeqKit 2.13.0 and R 4.4.3 with GenomicRanges 1.58.0, dplyr 1.2.1 and data.table 1.18.6.1. The upstream pipeline YAML pins a Linux-oriented environment and older CENdetectHOR library 0.0.2; this local macOS integration uses current library source and resolved local dependencies, so it is not a byte-identical replay of the YAML environment.

## Inputs and observed runs

`python benchmarks/scripts/generate_competitor_smoke_inputs.py benchmarks/evidence/competitor_local_smoke_20260916/inputs` reproduces all fixtures. Monomer A is `ACGTGCAATGTCAGTACCGTACGATCGTTA`; monomer B is `TTGACCGATGCTAGACCTGAGTCATCGTAC`; both are 30 bp. The direct HOR stage uses `cendetecthor_ab.fa` (12 AB pairs, 720 bp; SHA-256 `eaf4fa49c96325d93fddbfec7a1c317a3821a4ebdf62eee6f80c84e4ab504e1a`) and known `cendetecthor_ab.bed` (24 monomer intervals; SHA-256 `798a0604462f7916aa8925875c94fee31f86a4e5f9faeb4151ca2438d252c303`). The full pipeline input is `cendetecthor_pipeline.fa`, `((A+B)×9+A)×10`, 5700 bp, header `toy:0-5700` (SHA-256 `0acf542594ff333c69ffa5bcb4850c8b3a75fc4c0a828a4d6b6765a7215616c5`). The single-A insertion per block gives the upstream window filter more than one k-mer spacing value. All are toy sequences and lack biological truth.

The direct stage command was, from the official pipeline checkout, using the isolated Python environment:

```bash
<project>/.venv/bin/python scripts/findHORsFromMonomers_new.py \
  --fasta <evidence>/inputs/cendetecthor_ab.fa \
  --bed <evidence>/inputs/cendetecthor_ab.bed \
  --out-distMatr <evidence>/cendetecthor_stage_dist \
  --out-tree <evidence>/cendetecthor_stage_tree.xml \
  --log <evidence>/cendetecthor_stage.log \
  --t 1 --max-gap 10 --max-HOR 100 --min-loops 3 --dis-lev True
```

It exited 0. The archived `cendetecthor_stage_tree.xml` (SHA-256 `1db77f4fbe79c825b0dda0a0ad2e991a8467c7a8dfa893538361fb9613c837c0`) is parseable and has HOR clades `H` and `H1` with `F1#1,F1#2`. This test **forced known monomer intervals into the HOR stage**. It does not establish that the full pipeline discovers A/B or recovers the designed AB order.

The full pipeline invocation used a clean project work directory containing `fasta/toy.chr1.fasta`, copied byte for byte from the 5700 bp archived fixture, plus upstream `config/config.yaml`. From that directory:

```bash
env PATH=<project>/.venv/bin:<project>/.conda-r/bin:<project>/.conda-tools/bin:<stringdecomposer>/bin:$PATH \
  MPLCONFIGDIR=<project>/work/mpl TMPDIR=/private/tmp PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 \
  <project>/.venv/bin/snakemake \
  --snakefile <pinned-pipeline>/Snakefile --cores 1 \
  --shared-fs-usage input-output persistence software-deployment software-deployment-cache sources storage-local-copies \
  --config workdir=<project>/work/pipeline_run_complete/ consFile=no CHRcons=toy.chr1.fasta HORdet_threads=1 windowSize=1000
```

`<project>` is the absolute `benchmarks/competitor_envs/cendetecthor` path, `<evidence>` is the absolute evidence directory linked above, and `<pinned-pipeline>` / `<stringdecomposer>` refer to the exact pinned source checkouts in `/private/tmp/tandemx_cendetecthor_20260916/`. The explicit shared-filesystem list avoids Snakemake's default source cache, which attempted a nonwritable user cache path. The full run exited 0 with 11/11 Snakemake jobs complete. The native files are archived: `cendetecthor_pipeline_tree.xml` SHA-256 `27f9eaa0c54e361511e891124eb9a7c5314f4833bd323329e3b351f3b183a2e0`; raw decomposition TSV `3b8881e39e71fe8f924670325c12b8dbbf298184567f2acc4c53ee66374e5d24`; BED `8ee38b8ca58998d7edaa21cc6bbe7abee034beb22a54bc176a73aebda6d74712`; monomer FASTA `9db67dacb8efd4712a27c1d994659edd4643bce9e3a9a9688ec6af02bdd6b47a`. The phyloXML parsed despite Biopython `monomer_clade_seq` warnings.

**Structural truth outcome: failed.** The complete pipeline emitted one `toy_mon1` 60 bp consensus equal to A+B, and decomposed the input mainly into that 60 bp unit. It did not recover the designed separate 30 bp A and B monomers or their AB HOR as such. The 11/11 status means interface execution only. The direct stage's `H1` result cannot be counted as whole-pipeline accuracy.

## Retained failures and installation limits

| Attempt | Archived stderr/stdout | Observation |
| --- | --- | --- |
| Perfect 720 bp AB, `windowSize=120` | `cendetecthor_pipeline_stderr.txt`, `_stdout.txt` | 0 selected windows; `windowsFiltering.py` `IndexError` at `unik[0]` |
| Perfect 6000 bp AB, `windowSize=1000` | `cendetecthor_pipeline_long_stderr.txt`, `_long_stdout.txt` | Same empty-window `IndexError`; length alone did not fix constant spacing |
| 5700 bp variant, first run | `cendetecthor_pipeline_variant_stderr.txt`, `_variant_stdout.txt` | Reached 8/11, then `library(GenomicRanges)` missing; installing it in project `.conda-r` allowed the clean 11/11 rerun |

The failed perfect inputs are archived as `cendetecthor_pipeline_perfect_720.fa` and `_6000.fa`. No failed row is upgraded to a successful accuracy result.

[TideCluster](https://github.com/PetrNovak/TideCluster) current source HEAD was `8131547e75f73a1e61298a89ff4e861fbd613b7e` and README indicated 1.21.3. Historical local container setup is 1.21.2; Docker reported `Cannot connect to the Docker daemon at unix:///Users/ycy/.docker/run/docker.sock`. An isolated `mamba create --dry-run` for TideCluster 1.21.3 on `osx-arm64` found the noarch package but failed because required TideHunter 1.4.3 was unavailable for that platform. Thus there is no current TideCluster execution result. [HiCAT-human](https://github.com/865699871/HiCAT-human) source HEAD was `9e1f1697b5bdc2287e21375566153555261b6d8e`; its documented Linux/Python 3.9 and `lastz`/seqtk/StringDecomposer dependencies were inspected, but no executable was installed or run, and no matched human read/assembly dataset was enrolled. These states remain `not-run`, not negative biological results.

This log postdates the Phase-1 audit in [`tandemx_GR_gap_analysis.md`](../tandemx_GR_gap_analysis.md), whose “not installed/run” statements describe the earlier audit snapshot.
