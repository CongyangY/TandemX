# Local comparator availability and bounded smoke (2026-09-16)

This is an executable availability check on a **synthetic 360-bp FASTA**, not a formal accuracy, speed or same-endpoint benchmark. No T7 file was written or used. The exact input was `/private/tmp/tandemx_competitor_smoke_20260916.fa`, SHA-256 `217640e4a58e1fb8ebc4890d715b58ca8afa659a619aa7383e1f5b38a31e0e76`, one 30-bp motif repeated 12 times.

| Local executable | Version command/result | Executable SHA-256 | Toy command/result |
| --- | --- | --- | --- |
| `.conda-benchmark/bin/TideHunter` | `-v`: 1.5.5 | `b45940e18bd3855f0012305a27ce97fc863927f5479a9e379440da1b3037275c` | `-t 1 -p 30 -P 100 INPUT`: exit 0; one 30-bp consensus, nominal 12 copies |
| `.conda-benchmark/bin/trf` | `-v`: 4.10.0-rc.2 | `a1a23ca7fdbb4d6e5e61ca9913a80fb492486fb0ff7e6b567a7a2bb7c77ff982` | `INPUT 2 7 7 80 10 20 100 -ngs`: exit 0; one 1–360 interval, 30-bp period, nominal 12 copies |

The two commands were run separately with local preexisting executables. Their tool-specific output fields are not interchangeable and no cross-tool family or biological inference follows. `srf-load` exists in the environment but was not established as the SRF discovery executable. Current CENdetectHOR and TideCluster implementations were not located or installed in this check. Their historical/documented versions remain separate from a frozen current comparator run. The formal protocol still requires exact source/build hashes, shared input and task, resource budget, output normalization, and failure rows.
