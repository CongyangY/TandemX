# Exact-output discovery optimization replay

This archive records a real-data engineering replay of two changes to the
native elastic self-alignment path:

1. traceback directions use two bits per dynamic-programming cell rather than
   one byte;
2. candidate periods for one read are passed through one native call, so the
   uppercase read buffer is allocated once and reused across those periods.

Alignment scores, band widths, tie order, traceback, hit filtering, consensus
and downstream clustering are unchanged. The replay helper reran the complete
`tandemx discover` command from each historical baseline and required matching
input receipts and byte-identical discovery products.

| Mo17 input | Reads | Bases | Baseline time (s) | Optimized time (s) | Change | Baseline RSS (MiB) | Optimized RSS (MiB) | Change | Parity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 11 Mb | 840 | 11,680,888 | 17.348 | 12.361 | -28.75% | 81.703 | 62.484 | -23.52% | 6/6 core products byte-identical |
| 111 Mb | 8,084 | 111,505,681 | 134.684 | 100.039 | -25.72% | 178.438 | 154.766 | -13.27% | 7/7 products byte-identical |

The 11-Mb baseline predates `family_audit_summary.json`; that file was generated
by the replay but is listed as an uncompared auxiliary product. All six stable
core outputs still match. The 111-Mb baseline and replay both contain the audit
summary and all seven products match.

The original result locations on the originating workstation are:

- baseline: `/Volumes/T7/Codex/TandemX/results/Mo17_complete_random_11Mb_native_audit_v1_20260906`
- replay: `/Volumes/T7/Codex/TandemX/results/Mo17_11Mb_packed_trace_batch_replay_v1_20260907`
- baseline: `/Volumes/T7/Codex/TandemX/results/Mo17_complete_random_111Mb_indexed_clustering_v1_20260906`
- replay: `/Volumes/T7/Codex/TandemX/results/Mo17_111Mb_packed_trace_batch_replay_v1_20260907`

`comparison.tsv` is the compact derived table. Each scale directory retains the
baseline and replay environments, execution records, discovery summaries,
replay validation and stdout/stderr. `archive_manifest.json` records the source,
size and SHA-256 of every archived file.

The replay was run from an uncommitted optimization working tree whose parent
was `f734004`. Both replay environments therefore carry
`worktree_differs_from_git_head_use_source_digest`. They record source digest
`19600369d3a41ff36975a774d3a426ad2053ebe7d89947ecda6bdd39255b9907`
and individual source-file hashes, including the exact Rust and Python alignment
implementations used. This warning is retained rather than rewritten after the
fact. The four executed alignment source files match their versions in published
commit `a73398d`, whose work-branch/main hosted runs `34048998182` and
`34049010825` passed Ubuntu/macOS Python, Rust and wheel jobs.

Each row compares one historical baseline with one replay. Other jobs could
affect timings, and no repeated isolated timing distribution was collected.
These values support exact-output engineering improvement only. They are not a
publication-grade performance estimate, biological validation or evidence that
TandemX is faster or more memory-efficient than every external tool.
