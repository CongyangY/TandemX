# Frozen graph route on B1 v3 input-only bundle

The graph route source was frozen at commit `3e8ae23` after engineered B1 v2
development. It was run without opening B1 v3 truth or changing its exact-tile
adapter. All 39 B1 v3 cases returned `COMPONENT_UNRESOLVED`: 36
`tile_not_uniquely_exact_in_catalogue` and 3
`partial_or_nontiled_array_interval`. Native Col-CEN tile diversity is not
captured by the supplied operational three-monomer catalogue. This is a
method failure to produce decisions under the declared realistic development
input; performance must include all 39 abstentions.

Run: `python -m benchmarks.m2_routes.graph.run_common_development --inputs
benchmarks/controlled_collapse/v3/development_bundle/inputs.jsonl --predictions
benchmarks/m2_routes/graph/v3_predictions.jsonl --audit
benchmarks/m2_routes/graph/v3_input_only_audit.json`

The graph output does not establish whether the real array has a particular
HOR, order, or original-read support. B1 v3 read sequences are exact synthetic
copies of each source interval; original native reads are not paired.
