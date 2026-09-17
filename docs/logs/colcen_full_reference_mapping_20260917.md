# Col-CEN original-read full-reference mapping audit, 2026-09-17

The six original ERR6210723 HiFi records were selected during the earlier
three-context audit, **before** this complete-reference map. The mapping
protocol and command were committed as `b788469` before execution. It pins
the rechecked Col-CEN v1.2 132,081,078-bp, seven-record reference, the
six-read FASTQ, minimap2 2.31-r1302 executable SHA-256, `-x map-hifi -c
--secondary=yes -N 20 -t 2`, target coordinates and a fixed primary-spanner
eligibility rule. The reference gzip was read from the internal disk and
decompressed to a 134-MB file under `/private/tmp`; no large T7 write occurred.
Exact paths, hashes and command arguments are in
`benchmarks/controlled_collapse/native_genomewide_audit_20260917/protocol.json`.

The complete-reference PAF has one reported primary alignment for each of
the six selected reads and **no reported secondary alignment**. All six
primaries map to the chromosome/locus nominated in the local audit, span its
array plus >=1 kb target sequence on both sides, have identity >=0.99837111,
and pass MAPQ >=20. The lowest whole-reference MAPQ is 37. Counts are C1
2/2, C2 1/1 and C3 3/3. `all_alignments.paf`, `score/per_read.tsv`, the
machine-readable summary and SHA receipt preserve the exact evidence. The
scorer rejects changed source/reference/executable hashes, undeclared reads,
ambiguous primary states, malformed PAF fields and stale output paths.
A clean score replay under `/private/tmp/colcen_full_reference_score_replay_20260917`
matched all archived score files byte for byte (`diff -qr` exit 0); five focused
PAF validation tests passed.

This improves the locus-competition check for **these six locally selected
reads** from three extracted contexts to the supplied full Col-CEN assembly.
Minimap2's reported alignments and MAPQ are mapper-dependent; absence of a
reported secondary is not a mathematical proof of genome uniqueness. This
audit does not undo local read-selection bias or establish identical plant,
extraction, phased haplotype, physical copy number, or biological collapse.
The three C3 molecules remain one pooled read library and one assembly
lineage. Their use in the native deletion/equal-length pilots remains
development evidence only.

Minimap2's own stderr recorded 1.795 s real time, 2.418 s CPU and 0.856 GB
peak RSS while indexing this 132-Mb reference and mapping only six reads with
two threads. This is descriptive provenance for the anchoring check, not a
TandemX throughput/memory benchmark or a cross-tool speed comparison.

Recompute the score on this host after recreating the decompressed FASTA at
the protocol path:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.score_native_genomewide_mapping \
  --protocol benchmarks/controlled_collapse/native_genomewide_audit_20260917/protocol.json \
  --paf benchmarks/controlled_collapse/native_genomewide_audit_20260917/all_alignments.paf \
  --outdir /tmp/colcen_full_reference_score_replay
```
