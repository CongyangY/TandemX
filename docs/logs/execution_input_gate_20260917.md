# Execution input gate, 2026-09-17

Read-only live checks at 11:35 UTC, before new real-read experiments.

- Repository branch `codex/srf-submission-comparison-20260910`, HEAD
  `34d531f`; tracked tree clean. Existing `.codex/`, manuscript render trees
  and `tmp/` are user-owned untracked material and were not altered.
- This exact HEAD passed the prior local `tandemx-dev` Python suite (1,028
  tests) and GitHub Source validation run `35174663548` on macOS and Ubuntu.
  The next full suite is due after this execution round's code changes.
- `tandemx-dev` provides Python 3.11.15 and edlib; `rustc` and `cargo` are
  1.96.0 **inside** the conda environment. They are absent from the default
  shell PATH, so compiler checks must run through `conda run -n tandemx-dev`.
- Internal Data volume has about 21 GiB available (96% capacity used).
  Large public read downloads therefore require a bounded size budget.
- T7 is `/dev/disk5s2`, mounted as ExFAT with about 2.1 TB free. Its PSSD T7
  USB path currently negotiates **480 Mb/s** through a VIA USB 2.0 hub.
  Mount visibility and free capacity do not clear prior EIO or establish
  whole-volume integrity. No large T7 write or full-file verification was
  started by this gate.
- The committed B2 synthetic hold-out input and truth JSONL match their
  frozen receipt SHA-256 values: input
  `3cd49dded8921d3cfe15227ebce162ab509bb7243eb1105d995a1851984252bc`,
  truth
  `8076ed6b2563dabb34351dc56b65f4c05ef186ac41baf18f8c89e0a4008ab256`.
  These are integrity checks on synthetic benchmark artifacts, not native-read
  validation.
- The internal-disk Col-CEN v1.2 gzip copy independently matched SHA-256
  `b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8`,
  passed `gzip -t`, and streamed as seven nonempty FASTA records totaling
  132,081,078 bases (five nuclear chromosomes plus ChrM and ChrC). This
  validates assembly bytes/structure only; original-read pairing is a
  separate provenance gate.
- The existing methods benchmark protocol v1 already declares development,
  validation and final-heldout splits, task-specific input roles, donor and
  file-hash leakage guards, and a truth hierarchy. Its current four-item
  manifest remains an input-integrity inventory; it does not by itself enroll
  a paired original-read/assembly validation set.
- The v1 manifest validator was rerun read-only on its three enrolled files.
  Size, SHA-256 and full sequence structure passed for 1,000 ERR6210723
  prefix FASTA reads (15,666,956 bp), the T7 Col-CEN assembly (seven records,
  132,081,078 bp), and 673 YSD56 FASTQ reads (11,669,565 bp). Its known-bad
  K30076 attempt remains `invalid` with no file eligible for analysis.
  Machine output is archived at
  `docs/evidence/methods_execution_20260917/input_manifest_v1_recheck.json`.
  This bounded readback is evidence for those exact files, not T7-wide
  stability or source donor pairing.

Open gate: inspect real read/assembly source correspondence and archived
per-file receipts, then extract only bounded verified inputs for controlled
edits. K30076/V14167 and large 30x work remain separate T7-path blockers.
