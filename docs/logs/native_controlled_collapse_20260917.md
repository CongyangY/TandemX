# Native-read controlled-collapse development input, 2026-09-17

The source qualification is archived in
`benchmarks/controlled_collapse/native_pair_audit_20260917/`. The original
Col-CEN v1.2 assembly gzip SHA-256 is
`b059bf9b589a7a6cd13c67179b80293b91809c3c61cb8b9393a518619d8b5fa8`.
Three extracted 13,560-bp reference contexts match the previously recorded
C1/C2/C3 3,560-bp array sequences byte for byte; each nominated array occupies
context `[5000,8560)`. The combined context FASTA SHA-256 is
`638b0b0503b82f28b0c166a2b5351f2a2ced06bf4d8bed38de3657ffaaea4c67`.

The unchanged six-molecule FASTQ bundle SHA-256 is
`ff5387b874021554cebc53c4b4781f1f0f4739b8260fad03a9bf074b79e527a4`.
It contains two locally anchored C1 reads, one C2 read, and three C3 reads
from the original ERR6210723 HiFi sample. The source audit records the exact
read IDs, raw-record hashes, PAF, per-context alignment identity, flank spans,
and other-context hits. These are distinct ZMW read IDs from one movie; the
plants were pooled. The alignments were against three extracted contexts, so
their MAPQ and cross-context discrimination do not establish genome-wide
unique placement, identical assembly/read donor, or independent replication.
Even C3's three spanning molecules meet only a local technical count, not a
biological pairing gate.

The frozen input for exact edits is
`benchmarks/controlled_collapse/native_edit_development_v1/source_config.json`.
The editor `benchmarks/scripts/build_native_read_collapse.py` generated nine
cases per context: 0/25/50/75/100% terminal deletion, 25/50/75% internal
deletion, and 25% left-boundary deletion. For 20 operational 178-bp units,
these correspond to exact injected deletion sizes of 0, 890, 1,780, 2,670,
and 3,560 bp as applicable. The 100% case removes the nominated array while
retaining both natural 5-kb flanks. The 178-bp phase is an operational tile,
not independently established native monomer-copy truth. All 27 edited
contexts, source/edited coordinates, per-case FASTA hashes, and the original
read-bundle hash are in `generated/receipt.json`. No original read was edited.

Reproduce and independently check the edit archive with:

```bash
conda run -n tandemx-dev python -m benchmarks.scripts.build_native_read_collapse \
  --config benchmarks/controlled_collapse/native_edit_development_v1/source_config.json \
  --outdir /tmp/colcen_native_edit_replay
conda run -n tandemx-dev python -m benchmarks.scripts.verify_native_read_collapse \
  --config benchmarks/controlled_collapse/native_edit_development_v1/source_config.json \
  --generated /tmp/colcen_native_edit_replay \
  --out /tmp/colcen_native_edit_replay_verification.json
```

The independent verifier reconstructed every edited FASTA from the original
context and ledger coordinates and rechecked source/read hashes. The archived
verification reports **27/27 exact injected edits**, three contexts, and
read-bundle SHA-256 match. This validates the engineered assembly delta only.
An independent clean regeneration under
`/private/tmp/colcen_native_edit_replay_20260917` matched every generated
artifact byte for byte (`diff -qr` exit 0); both receipt SHA-256 values were
`c9839b7b7cbf690ead5fb5c93a5832807029cba6540393fcf3aa8786e13bae90`.
It does not establish physical missing-copy truth, a biological collapse call,
or M2 detection accuracy. The three contexts are one assembly lineage and
cannot be counted as three independent materials or final held-out samples.
