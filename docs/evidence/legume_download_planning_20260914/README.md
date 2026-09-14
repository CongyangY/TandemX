# Legume cohort source-integrity and download plan — 2026-09-14

This is a metadata-only checkpoint under the frozen real-30x protocol. All
repository queries used direct NCBI, ENA or NGDC routes without a proxy. No
sequencing or assembly payload was started. K30076 block-2 recovery was active
on T7 during the audit, so only small metadata pages, checksum manifests,
directory listings and run-stat XML files were fetched to temporary storage.

## Closed source records

The six peanut HiFi FASTQs still expose exact ENA bytes, base counts, read
counts, MD5 values and FTP paths. Their combined compressed FASTQ payload is
1,064,376,301,492 bytes. The six current NCBI assembly FASTA files add
3,859,022,794 bytes and retain the MD5 values in the official
`md5checksums.txt` files. This is an exact advertised payload inventory, not a
request to download all six materials.

The existing material-level linkage boundary remains unchanged: each peanut
assembly and raw run names the same accession, but the assembly and HiFi use
different BioSamples and the public record does not prove a single-plant DNA
extraction chain. V14167 already has hash-valid formal input and closed
TideHunter cells. K30076 block-2 attempt001 was actively recovering when this
receipt was written; block 3 had not started. S245, HN873, HN51 and S83 remain
planned and unstarted.

YSD56 remains the only complete local source among the small-genome controls.
Its ENA FASTQ and NCBI assembly metadata remain checksum-bearing and its
assembly and read use the same BioSample. ZH13's two HiFi BAMs expose exact
current bytes and NGDC run-stat XML now closes their combined 67,710,780,290
bases and 4,264,476 reads. They still lack a repository checksum and published
BAM index; the assembly and HiFi BioSamples differ, so ZH13 remains a
project/study control rather than an independently donor-matched validation.

A17 and R108 now have exact direct NGDC download plans. Their assembly FASTAs
are 129,989,665 and 116,010,804 bytes with official MD5 values. Their HiFi BAMs
are 20,993,038,389 and 25,092,297,745 bytes; tiny NGDC run-stat XML records
report 39,872,065,358 bases in 2,074,971 reads for A17 and 47,021,704,397 bases
in 2,334,551 reads for R108. The BAM directories publish no checksum or index,
so complete acquisition, self-computed SHA-256, EOF/record parsing and a BAM
converter gate remain mandatory.

YP4 is the lowest-friction new FASTQ control identified in this audit:
`SRR28983565` is a 24,590,641,333-byte ENA FASTQ with a repository MD5 and
31,747,954,482 bases. Its NCBI assembly FASTA is 167,792,293 bytes with an
official MD5. Assembly and read BioSamples differ but the study describes the
same YP4 plant material. Its evidence label remains
`T2T_claimed_partial_T2T_evidence` because only 9 of 11 pseudomolecules have
complete T2T support.

HJD's assembly file and MD5 are closed, but the direct NGDC BioProject page did
not yield a run accession in this session; raw acquisition remains
`metadata_blocked`. Williams 82 has an exact checksum-bearing ENA HiFi FASTQ,
but the Wm82.a5 Figshare file inventory returned HTTP 403 through the direct
API, so the assembly payload and donor chain remain unresolved. Lee has an
exact ENA FASTQ and NCBI chromosome-level assembly, but the two records use
different BioSamples and the assembly is highly fragmented; it is not an
eligible high-completeness control.

## Next safe transfers

No transfer should start while K30076 block 2/3 recovery holds the network and
T7. After that downloader stops and its aggregate receipts reach a safe
boundary, the smallest technically useful choices are:

1. YP4: 24,758,433,626 exact source-plus-assembly bytes, checksum-bearing FASTQ,
   with the partial-T2T label retained.
2. A17: 21,123,028,054 exact source-plus-assembly bytes, stronger complete-T2T
   evidence, but a checksum-less BAM requiring full acquisition and BAM QC.
3. R108: 25,208,308,549 exact source-plus-assembly bytes under the same BAM
   gate, with the existing taxon-identity reconciliation warning retained.
4. ZH13: 46,775,562,505 new raw bytes because its assembly is already local;
   checksum/index and donor-linkage gates remain unresolved.

The input staging estimate is deliberately reported separately from exact
download bytes. The current implementation can retain both a normalized 30x
FASTA and three partition FASTAs, so a conservative pre-output reserve is
`exact source + assembly bytes + 2 * target_30x_bases`, before tool outputs and
temporary files. Headers and line wrapping make this a planning estimate rather
than an exact disk prediction. The T7 snapshot during this audit showed
2,271,772,016,640 available bytes; it is volatile while K30076 recovery runs.

No family correspondence, accuracy, exclusivity or novelty result follows from
this inventory. Missing accessions and checksums are `N/A` or `unresolved`;
transport/client failures are retained as technical failures rather than zero
biological recovery.
