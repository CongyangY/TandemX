# TandemX assembly-only localization with supplied catalogue

On the frozen 36,371 bp assembly control used by TideCluster/CENdetectHOR,
`tandemx locate` was run with the control's `candidate_monomers.fa` as input.
This catalogue was derived from the engineered truth, so it gives TandemX a
strong prior that de novo TideCluster did not receive. It is an explicit
**catalogue-supplied condition**, not a fair cross-tool winner claim.

Command:

```bash
tandemx locate \
  --assembly benchmarks/assembly_hor_comparator_v1/control/control.chr1.fasta \
  --catalogue benchmarks/assembly_hor_comparator_v1/control/candidate_monomers.fa \
  --outdir /private/tmp/tandemx_assembly_control_catalogue_supplied
```

The command emitted 199 monomer interval calls versus 201 truth calls. Its
called base union covers 33,936/34,371 truth array bases (98.73% recall) with
no called bases outside the engineered array (100% base precision). A simple
truth-call check finds 176/201 with the same label and both boundaries within
10 bp of at least one call; this is not a one-to-one assignment. Native BED,
run config/log and a reproducible scorer are archived under
`benchmarks/assembly_hor_comparator_v1/`.

This route does not report HOR organization. The comparison with TideCluster
cannot rank algorithm accuracy because TandemX received truth-derived monomers,
and CENdetectHOR did not complete. A common no-prior HOR endpoint remains
unmeasured. The control consists of one engineered array, so even the
catalogue-supplied result is not a cross-genome accuracy estimate.
