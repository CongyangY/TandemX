# Two-flank array-length baseline (v1)

Frozen before the realistic B1 v3 development predictions. The input is the
same schema-v2 input-only JSONL bundle used by M2 routes. Each assembly and
read must contain one uniquely alignable left and right flank. edlib HW
alignment permits at most `max(2, ceil(0.10 * flank_length))` edits per flank;
multiple equally best placements are rejected. At least three read sequences
must span both flanks. The median read array span and assembly array span are
called discordant only when their absolute difference exceeds
`max(2, ceil(0.05 * max(read_span, assembly_span)))` bp. There is no tuning on
the v3 truth. One case represents one edited assembly, not three independent
observations merely because three read IDs exist. Event scores 0/1 are binary
decisions, not calibrated probabilities.

The optional signed bp delta uses `assembly span - median read span`, matching
the controlled-edit ledger convention. It is a span difference, not a validated
estimate of missing physical array sequence.

This baseline detects length differences, including copy loss and duplication.
It cannot detect equal-length order, inversion, or sequence replacement events.
It rejects partial reads and any case without uniquely placed paired flanks.
Its synthetic-flank performance does not establish real-read applicability.

Run: `python -m benchmarks.m2_routes.baselines.run_flank_length INPUTS.jsonl PREDICTIONS.jsonl`
