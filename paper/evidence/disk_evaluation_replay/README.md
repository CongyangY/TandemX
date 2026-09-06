# Completed disk evaluator checks

Source 8ff7ab6, 2026-09-06. These are evaluator checks, not new biological
accuracy results or complete discovery runs on the largest input.

- Mo17 111,505,681 bp / 8,084 reads: selected FASTA matches previous input.
  Normalized TSV bytes and all four descriptive metrics match retained TandemX,
  TRF and TideHunter outputs exactly. Complete replay: 4.117186 s / 43.328125 MiB.
- Col-0N 5,280,163,009 bp / 336,613 reads: full selected FASTQ parsing, checksum,
  duplicate-ID/totals verification and FASTA conversion succeeded above the old
  100,000-read cap. 161.180459 s / 54.59375 MiB. No tool execution or array
  normalization is asserted for this preparation-only run.

Resource receipts use direct-child wait4, including interpreter/input conversion;
other acquisition/benchmark jobs overlapped. This is a bounded-memory engineering
check, not an isolated performance ranking. SQLite indexes, FASTA and normalized
TSVs remain at the T7 paths in receipts. Source snapshots, exact hashes, logs and
reproduction command arguments are retained there; compact receipts are archived
here with byte-verified archive manifests.
