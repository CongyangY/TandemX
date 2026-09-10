# V14167 peanut execution-readiness note

**Date:** 2026-09-10
**Scope:** metadata, official source URLs, HTTP `HEAD`, and checksum-record audit only. No FASTQ/FASTA payload was downloaded and no TandemX analysis was started.

This note audits the frozen enrollment manifest [`benchmarks/configs/peanut_v14167_enrollment_manifest.json`](../benchmarks/configs/peanut_v14167_enrollment_manifest.json). The source paths are currently reachable, but the first pilot is not executable as a no-full-file operation with the repository's current sampler. This is a tool-contract limitation, rather than an access failure.

## Resolved public paths

The ENA read-run API returned one row for each requested V14167 run on 2026-09-10. The API query is archived as a reproducible source link:

<https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22SRR33996330%22%20OR%20run_accession%3D%22SRR33996240%22%20OR%20run_accession%3D%22SRR33996449%22&fields=run_accession%2Csample_accession%2Cstudy_accession%2Cinstrument_platform%2Cbase_count%2Cfastq_ftp%2Cfastq_md5%2Cfastq_bytes%2Cfastq_aspera%2Clibrary_strategy%2Cread_count&format=tsv>

All four official ENA HTTPS FASTQ objects returned `HTTP 200`, exposed `Accept-Ranges: bytes`, and had a non-zero `Content-Length`. The exact API values (bytes are compressed FASTQ bytes) were:

| Material/library | Run; BioSample | Exact bases | FASTQ object(s) and exact bytes | MD5 from ENA | Format/readiness |
|---|---|---:|---|---|---|
| V14167 HiFi | `SRR33996330`; `SAMN49012632` | 211,519,944,561 | `SRR33996330_subreads.fastq.gz`; 161,392,348,415 | `637ea54192dd85728c13037d8df97de3` | one gzip FASTQ; `HTTP 200`; range-capable |
| V14167 ONT | `SRR33996240`; `SAMN49012638` | 108,019,414,637 | `SRR33996240_1.fastq.gz`; 99,089,152,664 | `54379601091f290337376e9ab342de6e` | one gzip FASTQ; `HTTP 200`; range-capable |
| V14167 DNBSEQ | `SRR33996449`; `SAMN49012650` | 131,113,606,800 | `_1.fastq.gz`; 47,199,787,717; `_2.fastq.gz`; 50,166,933,327 | `7a0e8bf6886f5149f47db72c9b5777ca`; `f6022778cad0e27ab2007577c6e7bab5` | two gzip FASTQs; both `HTTP 200`; paired cardinality and MD5 cardinality agree |

The direct HTTPS paths are deterministic ENA paths derived from the API `fastq_ftp` field:

- [HiFi FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/030/SRR33996330/SRR33996330_subreads.fastq.gz)
- [ONT FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/040/SRR33996240/SRR33996240_1.fastq.gz)
- [DNBSEQ mate 1](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/049/SRR33996449/SRR33996449_1.fastq.gz)
- [DNBSEQ mate 2](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/049/SRR33996449/SRR33996449_2.fastq.gz)

The runs are public and carry raw study `PRJNA1259747`. The NCBI assembly report records assembly BioProject `PRJNA1259200`; this is the known project alias retained in the frozen manifest, not a run-to-sample mismatch. The assembly report also identifies isolate `V14167` and BioSample `SAMN48356820`.

The assembly directory is public:

- [NCBI assembly record](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824555.1/)
- [NCBI FTP directory](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/824/555/GCA_054824555.1_ASM5482455v1/)
- [Assembly FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/824/555/GCA_054824555.1_ASM5482455v1/GCA_054824555.1_ASM5482455v1_genomic.fna.gz)
- [NCBI MD5 record](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/824/555/GCA_054824555.1_ASM5482455v1/md5checksums.txt)
- [Assembly report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/824/555/GCA_054824555.1_ASM5482455v1/GCA_054824555.1_ASM5482455v1_assembly_report.txt)

The FASTA returned `HTTP 200`, `Accept-Ranges: bytes`, and `Content-Length: 339,638,203` bytes (about 339.638 MB decimal, listed as 324M in the directory view). The official NCBI MD5 is `546093beba44a0517d3cbfc4bafd819d`. The report identifies `ASM5482455v1`, *Arachis duranensis*, isolate V14167, and assembly BioSample `SAMN48356820`. Its `Assembly type: haploid` is the representation of the AA reference and is not evidence for a different donor.

## Payload and storage budget

The three raw libraries sum to **450,652,965,998 reported bases** and **357,848,222,123 compressed bytes (357.848 GB decimal)**. Adding the compressed assembly gives **358.187860326 GB decimal** of source payload. The frozen manifest stores rounded three-decimal summaries; the integer values above are the current ENA API/HTTP audit values.

The proposed pilot asks for 10.61 Gb decimal bases in total: 5.89 Gb HiFi, 2.36 Gb ONT, and 2.36 Gb DNBSEQ. The equivalent first sampling probabilities, if the existing fraction-based sampler is used, are approximately 0.027846, 0.021848, and 0.018000, respectively. These are Bernoulli inclusion probabilities, so they do not guarantee an exact base target; the receipt must report observed bases.

Compressed bytes for those small base targets **cannot be predicted safely** from metadata alone. A gzip prefix is not a read-boundary-preserving sample, and compression ratio varies with read length and sequence composition. Do not reserve a guessed “10.61-Gb FASTQ” byte amount as an execution guarantee.

If complete source files are staged for the current repository workflow, retain at least two copies during transfer/QC (partial plus final): about **715.696 GB** for all raw libraries, plus **0.679 GB** for two assembly copies, before QC databases, sampled outputs, and filesystem headroom. A single-file staging budget is approximately 322.785 GB for HiFi, 198.178 GB for ONT, or 194.400 GB for the paired DNBSEQ run under the same two-copy rule.

## What can and cannot run without full files

The repository contract in [`docs/complete_data_qc.md`](complete_data_qc.md) requires `fetch_ena_complete` to verify the complete ENA file size/MD5 before publication, and `qc_complete_fastq` to parse every FASTQ record through EOF. `sample_complete_fastq` then requires a successful full-file QC receipt and again scans the source through EOF while recomputing totals and the source SHA-256. Therefore the frozen pilot's deterministic, source-order-preserving sampling **cannot run against only remote metadata, a prefix, or an arbitrary HTTP range**.

ENA's `Accept-Ranges` is useful for a future resumable transfer, but these are ordinary `.fastq.gz` objects without a read index. A byte range does not map to complete FASTQ records and cannot provide the manifest's validated primary-read ID hash semantics. Implementing a remote indexed/streaming sampler would be a separate code change and is outside this read-only audit.

The source is consequently **metadata-ready and acquisition-ready for a later approved bounded staging step**, while the no-full-download pilot is `blocked_by_current_tool_contract`. No download is authorized or implied by this note. If the pilot is approved later, archive the ENA API row first, run the repository fetch/QC workflow one complete run at a time, and only then use the frozen seeds/fractions; do not replace the MD5/size gate with a HEAD, ETag, or guessed byte prefix.

## Lightweight checks and future commands

The following checks are sufficient for the metadata phase and do not transfer the payload:

```bash
curl -L --fail 'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22SRR33996330%22%20OR%20run_accession%3D%22SRR33996240%22%20OR%20run_accession%3D%22SRR33996449%22&fields=run_accession%2Csample_accession%2Cstudy_accession%2Cinstrument_platform%2Cbase_count%2Cfastq_ftp%2Cfastq_md5%2Cfastq_bytes&format=tsv'
curl -sS -I -L 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/030/SRR33996330/SRR33996330_subreads.fastq.gz'
curl -sS -I -L 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/824/555/GCA_054824555.1_ASM5482455v1/GCA_054824555.1_ASM5482455v1_genomic.fna.gz'
```

After a separate approval for complete-file staging, the repository workflow is:

```bash
python -m benchmarks.scripts.fetch_ena_complete \
  --metadata ena_metadata.tsv --accession SRR33996330 \
  --outdir /path/to/V14167/SRR33996330_complete \
  --max-download-gb 162
python -m benchmarks.scripts.qc_complete_fastq \
  --fastq /path/to/V14167/SRR33996330_complete/SRR33996330_subreads.fastq.gz \
  --outdir /path/to/V14167/SRR33996330_complete/qc \
  --expected-bases 211519944561
```

The example is deliberately a future command: it was not run here. For DNBSEQ, the archived metadata row must retain both semicolon-separated URLs, MD5 values, and byte values. The complete-run download cap must exceed the largest individual compressed object and the selected run's total, while storage reservation must follow the two-copy figures above.

## Stop rules for the next execution stage

Stop before any payload transfer if an ENA query returns zero or duplicate rows, a run/BioSample/study/material or platform field changes, a FASTQ URL is not an official ENA HTTPS path, file/MD5/byte-list cardinalities disagree, a HEAD is non-200 or zero length, or the source is restricted/withdrawn. Stop after staging if the complete-file MD5 or byte size differs, gzip/FASTQ parsing fails, DNBSEQ mates are incomplete, assembly MD5 is not `546093beba44a0517d3cbfc4bafd819d`, or the approved disk budget is exceeded.

Treat a difference between rounded manifest values and the integer API values as expected rounding, not a source change. Treat the assembly report's `PRJNA1259200` alias and `haploid` AA-reference representation as known metadata context, not donor-match proof. A successful ingestion/QC decision also does not make the assembly absolute copy-number truth.
