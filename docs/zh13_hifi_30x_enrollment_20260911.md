# ZH13 HiFi 30x enrollment

Status: metadata and assembly input are closed; raw-HiFi acquisition is not
started. This is a benchmark-control enrollment, not a result or an extension
of the frozen estimator.

## Corrected raw-read identity

The GWH record `GWHBWDJ00000000.1` identifies ZH13-T2T as a complete 20
chromosome assembly in `PRJCA015269`, with GWH BioSample `SAMC1127443`.
The local official FASTA has already passed its publisher-MD5, gzip/FASTA and
independent sequence-length checks: 1,007,237,669 A/C/G/T bases in 20 records.

`CRR705247` is **not** the HiFi run. It is the ONT WGS run
(`soybean_ONT_genome`) from `SAMC1127443`. The correct HiFi experiment is
`CRX631709` (PacBio Sequel II, WGS, `SAMC1127444`,
`soybean_leaf_DNA_HiFi`) and has two public BAM runs:

| run | public object | bytes verified by direct HTTPS HEAD | repository checksum |
|---|---|---:|---|
| CRR705248 | `CRR705248.bam` | 24,194,765,935 | not shown |
| CRR705249 | `CRR705249.bam` | 22,580,796,570 | not shown |

The official project and cultivar/study linkage are adequate for a controlled
positive cohort, but the assembly and HiFi records use different BioSample
accessions. A valid leading BGZF block from CRR705248 identifies the CCS
library and BAM read-group sample as `22D00017`, not `ZH13`. The public records
do not close the individual DNA-extraction chain or map this internal label to
the cultivar. This cohort must therefore not be presented as independent,
donor-matched HiFi validation.

## Acquisition and partition gate

The two BAM objects advertise `Accept-Ranges: bytes`; a direct no-proxy,
1-MiB request of CRR705248 returned exactly HTTP 206 with the requested range.
No standard `.bam.bai` or `.bai` object is public at the corresponding paths.
Thus an arbitrary byte range cannot be treated as a valid HiFi subset.

No full source transfer has started. The 46,775,562,505-byte total must first
be acquired directly to T7 and pass length, BAM EOF/record parsing and a
self-computed SHA-256 audit. A 1,000-read FASTQ conversion gate is then
required. Only after that gate, and only when it does not materially slow the
active V14167 or YSD56 work, may three disjoint source-order partitions be
created. Each targets 10,072,376,690 bases (10x against the verified ZH13
assembly); the final receipt will report actual bases rather than target depth.

The precise URLs, ETags, transfer rule and deterministic partition rule are
frozen in `benchmarks/configs/zh13_hifi_30x_enrollment_v1_20260911.json`.
