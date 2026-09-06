# Divergence-aware localizer development v2

Development v2 reused only the already consumed seeds 5301--5303 after v1 failed.
It kept the 0.90 IID base-identity proxy and allowed exact k-mer anchors to bridge
across at most one monomer length. The planted 500-bp interruptions remain beyond
that bound for the tested 61/171/421-bp monomers.

All 90 localization commands completed. Full-assembly mean base recall was
0.975728, positive-assembly mean precision was 0.999441 and none of 54 absent-
family rows had predicted repeat bases, passing all three development gates.
Maximum predicted fragments decreased from 119 in v1 to 4. The weakest
full-assembly stratum remained 5% divergence with three segments (mean recall
0.941017), which is below the aggregate 0.95 gate and is retained explicitly.

`archive_manifest.json` verifies 24 compact files, including all localization
rows, 90 receipt hashes/resource measurements and the frozen source snapshot.
The hashes of validation, metrics, environment and configuration were committed
in `abundance_localizer_heldout_v1.json` before fresh seeds were used.
