# Phase-gate fixture interface

`benchmarks.challenge.phase_gate_fixture.generate_case(seed_label, condition, coverage)`
creates one deterministic, development-only toy accuracy-gate case. `seed_label`
must begin with `dev`; labels containing `holdout` are rejected. `condition` is
one of `clean`, `substitutions`, `indels_small`, `indels_long`, `partial`,
`interruptions`, `background_homology`, `close_families`, `heterogeneous` or
`abundance_low`. `coverage` is 2, 5 or 10 and controls sampled 1,000-bp reads.

The returned dictionary contains `catalogue`, `reads`, `truth_intervals`,
`source_truth_intervals`, `source_truth_bp`, `read_truth_bp`, `read_bases`,
`genome_size`, `condition`, `coverage`, `source_hash` and `reads_hash`.
`truth_intervals` use `(read_id, start, end, family)` half-open read coordinates;
source intervals use `(start, end, family)` coordinates. Truth is created from
the source genome masks before any analysis algorithm runs. Inserted bases in a
repeat array are masked as repeat; the 300-bp interruption is explicitly
unmasked. The 24-bp/80-bp homology background is not complete-family truth.

This fixture is for development negative controls and mechanical validation of
accuracy-gate plumbing. It is not a future holdout, benchmark result or claim
about real sequencing data. The module returns data in memory; callers should
avoid writing large files or selecting development parameters from these cases.
