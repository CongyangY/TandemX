# Formal historical comparator terminal evidence

This package reconstructs the terminal bookkeeping state of the accepted
historical SRF/ordinary-mapping run without rerunning a native stage or changing
any scientific rule.

The preserved formal directory is:

`/Volumes/T7/Codex/TandemX/results/historical_prioritization_srf_mapping_v1_20260913`

Its frozen `run_config.json` and source snapshot remain unchanged. KMC k151
counting and `kmc_dump` completed. The first Ey15-2 SRF assembly stage then hit
the preregistered 7,200-second timeout: exit `-9`, `timed_out=true`, measured
runtime 7,201.486 seconds, peak RSS 6,814.375 MiB, and a zero-byte `srf.fa`.
This is a technical failure, not a zero-family result.

After catching that native failure, the frozen controller attempted to emit 19
eligible-family terminal rows. Each row included
`native_retained_read_bp=N/A`, while the explicit TSV field list omitted that
column. Python `csv.DictWriter(extrasaction="raise")` therefore raised the exact
error retained in `controller_terminal.stderr.txt`. The formal directory has a
header-only 147-byte family TSV and no root `summary.tsv` or `completion.json`.

The reconstructed terminal state is:

- `ey15_2/srf_k151`: `technical_failure`, all 19 eligible family estimates and
  proxy calls are `N/A`;
- the remaining five planned cells: `not_run_prior_failure`;
- overall completion: `false`; no proxy-agreement metrics are available.

`formal_artifact_inventory.json` hashes every preserved formal file.
`reconstruction_receipt.json` binds the frozen config, source snapshot, native
stage receipts, captured controller traceback, and header-only TSV. The
reconstruction script refuses a formal root that already has terminal summary
files and never invokes KMC, SRF or minimap2.

The workspace serializer now declares `native_retained_read_bp` in its stable
terminal TSV schema. A regression test exercises the exact technical-failure
row. This future-run code fix does not alter or retroactively replace the
formal source snapshot and is not authorization to rerun the comparator.
