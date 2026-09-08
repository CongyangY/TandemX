# TideCluster 1-Gb MorexV3 resource gate v1

This archive preserves the frozen preflight decision for the nested 1-Gb
MorexV3 reference sample. The decision rule was defined before execution: a
candidate run is refused when the maximum successful 100-Mb stage already uses
at least 85% of the same Docker runtime memory limit.

The completed 100-Mb TideCluster clustering stage used 7,957,438,464 bytes
(7,770,936 kB), or 96.836% of the 8,217,432,064-byte Docker limit. The 1-Gb
run was therefore not started and its retained fate is
`resource_infeasible_preflight`, not a missing or zero-valued measurement.

This establishes infeasibility only for the measured host and pinned container
runtime. It does not establish that TideCluster cannot process 1 Gb on a host
with more memory.
