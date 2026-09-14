# Legume real-30x cohort re-audit — 2026-09-13

This audit resumed only after the user reported that the T7 volume had passed a
read-only volume check. It followed the frozen real-30x protocol and did not
launch BLAST, add a sixth manuscript Result, alter the TandemX production
algorithm, or start any uncontrolled source download.

Direct public-database checks were repeated against ENA, NCBI SRA RunInfo,
NCBI Datasets, NCBI assembly `md5checksums.txt`, and NGDC HTTPS headers. The six
peanut HiFi runs remain public PacBio Sequel II WGS records in
`PRJNA1259747`. Their SRA sample names bind the runs to V14167, K30076, S245,
HN873, HN51 and S83. Assembly and HiFi records use distinct BioSample
accessions, so the evidence is material-level rather than a proved single-plant
DNA-extraction chain. The four tetraploid spans in the earlier recovery
manifest were approximate. This audit replaced them with exact current NCBI
assembly spans and recomputed the 30x targets.

The ZH13 assembly remains locally complete and re-passed its publisher MD5
`86fe5e3cdd12510afa61c6e3345a4a97` and local SHA-256
`db628d16c550e4f587bf834eea36b2c2325afcbe633923161e04312fa8f32ca9`.
Direct NGDC HEAD requests for `CRR705248.bam` and `CRR705249.bam` again returned
24,194,765,935 and 22,580,796,570 bytes, respectively, with the frozen ETags.
No repository checksum or standard BAM index is published. Complete BAM
acquisition and QC therefore remain required; no byte-range substitute was
used.

The complete YSD56 source FASTQ re-passed ENA MD5
`bcec44c711a45dd7dd22784e0da7cccd`. Its common 30x FASTA and three formal
partitions re-passed every SHA-256 recorded in the frozen run manifest. Existing
TandemX and TideHunter executions completed on all three partitions, while TRF
timed out at the frozen two-hour limit on all three. The real three-tool endpoint
remains unresolved.

V14167's 29.940948785x common FASTA and all three formal partitions also
re-passed their recorded SHA-256 values. Its block-02 and replacement block-03
aggregate receipts remain valid. Block-01 has no aggregate receipt and is kept
as `recoverable_invalid` at that source-block layer, although the selected
per-chunk receipts and their hashes are preserved in the completed input
manifest. K30076 block-01 remains `recoverable_invalid`: the payload is nonzero
but its aggregate receipt is zero bytes. S245, HN873, HN51 and S83 have no local
30x input and remain unstarted.

The frozen V14167 TideHunter 1.5.5 recovery is now complete for all three
hash-verified disjoint partitions. Runs were strictly serial, used four threads
and the frozen 30--1000-bp period scope, and wrote to new independent directories.
All three exited 0 without timeout and passed independent validation of the
eleven-field rows, partition read IDs, coordinate bounds, nonempty consensus and
input/output hashes. Together they contain 5,652,836 native calls, of which
4,395,253 satisfy the frozen period and span scope. Total launcher runtime was
7,763.35 seconds; the maximum per-partition peak RSS was 2.844 GB. These are
operational records, not biological family counts. The durable aggregate receipt
is `/Volumes/T7/Codex/TandemX/results/peanut_v14167_30x_comparison_v1_20260912/formal_step_resume_v1/tidehunter/partition_summary_receipt.json`.
The interrupted earlier `normalization.sqlite` is malformed and was not reused.
TideHunter is technically closed for V14167, but the frozen three-tool benchmark
remains incomplete. No family correspondence, accuracy, exclusivity or novelty
conclusion is made.
