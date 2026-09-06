# Divergence-aware localizer held-out baseline

This archive records the first and only use of predeclared seeds 5501--5503.
Configuration `abundance_localizer_heldout_v1.json`, the development evidence
hashes and a guard that recomputes all selection gates were committed as
`4e662db` and passed hosted Ubuntu/macOS CI on both the working branch and main
before execution. The matrix retained the six divergence/structure scenarios,
nine coverage/read-error conditions and five assembly retention levels.

All 1,062 commands completed: 90 localization, 162 quantification and 810
comparison commands. The archive contains 486 copy-number, 270 localization and
2,430 comparison rows. Full-assembly mean base recall was 0.981033, positive-
assembly mean precision was 0.999565 and none of 54 absent-family rows had
predicted bases. The weakest full-assembly stratum was 5% divergence with three
segments (mean recall 0.949274). Positive-assembly mean absolute repeat-bp error
was 0.024861; one short-array condition reached 0.596721 absolute relative error.

The single-k assembly/read rule yielded TP/FN/FP/TN=875/583/12/960:
sensitivity 0.600137, false-positive rate 0.012346 and precision 0.986471.
Localization therefore passed its aggregate gates while the downstream baseline
remained insensitive under divergent read-copy estimation.

`archive_manifest.json` verifies 28 compact files and every receipt. Summed
direct-child wall time was 133.751 s; the maximum direct-child RSS was
34.188 MiB. These values describe sequential toy commands, not production-scale
throughput. The compact archive is about 1.6 MiB; the approximately 1.45-GiB T7
run retains generated reads, assemblies and complete per-command outputs.
