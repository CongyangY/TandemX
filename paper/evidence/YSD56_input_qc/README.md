# Complete YSD56 HiFi file QC

ENA run SRR28726931 is reported as PacBio Revio genomic HiFi WGS from the
YSD56 leaf sample of wild soybean (*Glycine soja*). The complete compressed
FASTQ is 21,866,640,310 bytes; its ENA MD5
`bcec44c711a45dd7dd22784e0da7cccd` and local SHA-256
`2f3565db063e26638b8cf0518ede051300f9af09b504f2eca0c52c0bd7b23975`
were verified before QC.

Full-file streaming QC passed for 2,617,227 unique read IDs and
44,193,089,411 bases. FASTQ records and the complete gzip stream are valid;
the observed counts equal ENA metadata. Median length is 16,224 bp, read N50
17,474 bp, GC fraction 0.3474941751 and N fraction zero. The Phred-derived mean
reported error probability is not empirical read accuracy. File QC does not
establish species purity, donor matching, ploidy, nuclear coverage or repeat
truth.

The source study reports assembly GCA_040083835.1 for YSD56. Its exact
1,008,523,555-bp assembly length is used only as the declared denominator for
the deterministic seed-6101 sampling ladder. Sampling is a technical scale
experiment, not an additional biological replicate. `archive_manifest.json`
hashes the acquisition plan/receipt and the four compact QC files; the raw
FASTQ and SQLite duplicate-ID audit remain at the T7 data root.
