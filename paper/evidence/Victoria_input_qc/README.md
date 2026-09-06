# Victoria oat complete input QC

Public HiFi run ERR10422581, sample SAMEA111508775; source/raw-to-assembly
curation is in `paper/evidence/PanOat_source_audit`. The original article's
supplement confirms Victoria's HiFi and Hi-C, with assembly accession
GCA_947311595. This included run is not the full oat pangenome cohort.

The complete6,980,726,507-byte FASTQ passed official MD5 verification
(`783eaba6ed6861ebc62abdde71fd1c85`), full gzip/record validation and exact
duplicate-read-ID checks. Observed398,850 reads /7,346,159,178 bp;
median18,137 bp, N5018,295 bp, length range46–28,424 bp, GC44.1718%, N0,
duplicate IDs0. Mean reported base-error probability0.00240421 is a Phred
summary, not empirical sequencing accuracy. Full-file SHA-256:
`116bc62711ff6de969b7ed34b30a8e21b5e6bc4fcd5c8c5d194525e9bff4022e`.

Source receipts and QC histograms are byte-preserving archived with hashes.
The download receipt's `fastq_qc_status=not_run` records its earlier transfer
stage; the separate completed `qc_v1/qc.json` is authoritative for later QC.
Read selection uses whole-library seed6101 hashing. Completed nested samples
contain615/5,972/60,215 reads and11,361,028/110,202,741/1,109,305,135 bp.
Their receipts and distributions are archived under `sampling`. Data remain under
`/Volumes/T7/Codex/TandemX/data/raw/ERR10422581_complete`.

These checks establish file integrity and observed distributions. Reference
concordance, contamination, empirical bias and family-level accuracy require
additional biological analyses. Technical subsets are not independent plants.
