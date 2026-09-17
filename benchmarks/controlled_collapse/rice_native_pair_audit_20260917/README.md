# Nipponbare rice bounded native source screen: negative (2026-09-17)

The selection protocol was committed in `65cb940`, then bounded to the top
1,000 ranked windows and unit-tested in `16538ca`, **before mapping or observing
any native read at a candidate interval**. It fixes the AGIS1.0 assembly,
SRR25241090 Bernoulli `sample_003`, a 155-bp operational shift, 4,000-bp
windows at 500-bp stride, 3,000-bp natural flanks, and three quality thresholds.
The period is an operational screen, not an asserted biological monomer.

The current full 120,079,535-byte assembly gzip SHA-256 matches the earlier
receipt (`cbe51d6615380bb441b7cf531ebc0ae58f498bc7982d2db07c549376a7638270`).
Its 12 chromosomes total 385,710,679 bases. The current complete readback of
the bounded 1,071,470,526-byte original FASTQ subset matches SHA-256
`88e05ea73487dffa014ca3f8bf87881eb58071635fcde65a75349264d20a821f`;
gzip trailer, FASTQ framing, sequence/quality lengths and 62,345 unique read
IDs passed. The 62,345 records contain 1,150,127,844 bases.

The full assembly generated **771,185** windows. The frozen ranking reviewed
the first **1,000**. Among those, 518 reached shift identity 0.85 and all
1,000 reached entropy 1.6, but none reached the required 0.8 within-context
flank unique-31-mer fraction (maximum 0.638047). The highest shift identity
was 0.971131 at CP132245.1 `[13744000,13748000)`, whose flank uniqueness was
0.034848. The frozen selection returned `null`. This is a negative result for
**this declared screen**, not proof that no rice tandem interval is usable.
No native reads were mapped, no candidate was reselected, and no controlled
edit was generated. The full denominator and all 1,000 reviewed rows remain
in `selection_receipt.json` and `ranked_screen.tsv`.

The assembly BioSample SAMN36344332 and read BioSample SAMN36368305 share a
Nipponbare/AGIS-1.0 label but have different archive identifiers and
collection dates. Exact individual and DNA-extraction pairing is unresolved.
Even a successful technical edit would not establish natural physical copy
number. This source is development only and is not enrolled as final-heldout.

Reproduce the bounded check inside `tandemx-dev`:

```bash
python -m benchmarks.scripts.screen_rice_native_interval \
  --protocol benchmarks/controlled_collapse/rice_native_pair_audit_20260917/protocol.json \
  --outdir benchmarks/controlled_collapse/rice_native_pair_audit_20260917
python -m benchmarks.scripts.audit_rice_native_subset \
  --protocol benchmarks/controlled_collapse/rice_native_pair_audit_20260917/protocol.json \
  --out benchmarks/controlled_collapse/rice_native_pair_audit_20260917/subset_integrity_receipt.json
pytest -q tests/unit/test_screen_rice_native_interval.py tests/unit/test_audit_rice_native_subset.py
```
