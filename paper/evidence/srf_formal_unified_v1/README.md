# Frozen SRF unified comparator evidence

72/72 cells completed, zero process failures, 174.336 s controller elapsed.
Protocol/source commit: 7abafc3. Production remains baseline 81827c3.
Native SRF k=151: 15 ok + 3 no_catalogue; k=101: 18 ok.
TandemX and competitive mapping each completed 18 cells. No A3 was run.

`condition_summary.tsv` contains arithmetic means and min/max over three
independent simulated datasets for each condition/method. They are not three
technical repeats. Mapping-only cost and shared-discovery-plus-mapping cost
are separate fields. Peak RSS excludes the controller; stored bytes are not
peak temporary disk. No large-input speedup inference is made.

MARE scores positive-truth families; false/unmatched catalogue abundance is
retained in `unassigned_abundance_bp` and must be read alongside MARE. In the
shared-fragment challenge, near-perfect planted-family MARE for competitive
mapping does not imply high catalogue specificity. Global precision includes
unmatched and ambiguous native intervals. Background bp measures interval
attribution, not where diagnostic-kmer abundance was allocated.

Raw inputs, native PAF/BED/count databases, exact commands, stderr and complete
source snapshots remain under:
`/Volumes/T7/Codex/TandemX/results/srf_formal_unified_v1_20260910`.
This compact archive includes all method receipts, family amounts, correspondence,
interval summaries, truth manifests, input truth and environment/build identities.
The original receipts bind larger raw products through SHA-256.

The full clean KMC build/verification archive is:
`/Volumes/T7/Codex/TandemX/tools/build_receipts/kmc_clean_20260910_archive`.

Reproduce aggregation into a new directory with:
`python -m benchmarks.scripts.summarize_srf_formal --results FULL_RUN --archive NEW_ARCHIVE`.
Historical real-assembly SRF/mapping comparisons remain not_run.
