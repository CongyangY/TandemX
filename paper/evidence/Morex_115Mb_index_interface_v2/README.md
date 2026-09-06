# Morex 115-Mb sequence-native interface diagnostic

This compact archive records a fixed-order ablation on the 7,094 previously
detected Morex candidates from the 115,272,341-bp read sample. It reuses the
candidate table and monomer FASTA and does not rerun read scanning. Both fresh
children used the same current clustering source and native extension; only the
Python word-list bridge versus complete-sequence Rust interface differed.

Both paths produced 4,380 families and byte-identical complete payload hashes.
The word bridge used 13.610636 clustering seconds, 13.895589 child seconds and
61.203125 MiB peak child RSS. The sequence-native path used 12.198283,
12.473152 and 58.046875, respectively: 10.38% less clustering time, 10.24% less
child wall time and 5.16% less peak RSS in this run.

The snapshot records Git commit `a10c309` and exact hashes for the two local
native binaries. Its generic revision warning is caused only by those ignored,
untracked build artifacts; no tracked source file differed from the commit.
This is one fixed-order engineering run, so it is evidence of exact output
parity and a favourable diagnostic, not publication-grade repeated timing.
The earlier v1 measurement is retained separately and differed in magnitude,
which reinforces the need for randomized repeated isolated runs.

`archive_manifest.json` covers 14 compact files. The two multi-megabyte duplicate
JSON payloads remain at the T7 data root and are represented by their common
SHA-256 in the validation and per-child receipts.
