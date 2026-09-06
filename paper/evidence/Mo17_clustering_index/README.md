# Exact clustering-index verification

New source e57ae542f1d3ac50ada057a8d6155821f358fb5c. The canonical multiset
overlap bound prunes only pairs rejected by the existing oriented q-gram test.

`serialized_111Mb` replays 8,401 exported candidates through original-snapshot
and indexed clustering. Both produce 4,219 families and identical complete family
and membership JSON (SHA-256 dd299e270a4316c5350083761cdec378d17954cfcc26bd5f62868aa51e43adf2).
Stage times are 13.131848 and 5.293893 s, a single-process diagnostic pair.
Exported scores are rounded; this alone is not full live-pipeline parity.

`full_pipeline_111Mb` independently reruns the actual public command from reads.
All seven deterministic data products match the previous related-audit run byte
for byte: candidates, candidate sequences, representative sequences, families,
monomer membership, similarity table and audit summary. Checks are recorded in
`full_pipeline_111Mb_parity.json`. The underlying input is the same 8,084 reads /
111,505,681 bp whole-library hash sample.

Full rerun: TandemX 134.684 s / 178.44 MiB, TRF 88.068 s / 258.03 MiB,
TideHunter 33.352 s / 328.53 MiB. Previous same-input TandemX observation was
142.919 s / 155.45 MiB. Thus observed pipeline time fell but observed peak memory
increased; the index has a memory trade-off to investigate. All timings and memory
values are concurrent-job diagnostics, not final isolated performance estimates.
TandemX is still slower than both comparators on this input. All counts are
descriptive because independent real repeat ground truth is not supplied.

The older 1.129-Gb run is a separate ongoing baseline; do not infer its completion
or memory from this smaller experiment. Source/input/run/output provenance is
retained; native outputs remain at their T7 paths.
