# CENdetectHOR on the frozen flanked 171-bp control (2026-09-17)

The same 36,371-bp corrected-header FASTA used in the prior CENdetectHOR and
TideCluster assembly comparison was replayed without supplying monomers or the
designed 171-bp period. Its SHA-256 is
`f2a968383dd9bef0af38ee5891ce97b7228f67cbe71495e36340adc20d956e8e`.
The designed array is `[1000,35371)`, but no truth file was passed to the
native pipeline. The complete commands, frozen setting rationale, process
traces, stderr, native intermediate files and checksums are in
[`cendetecthor_171_precondition_20260917`](../../benchmarks/evidence/cendetecthor_171_precondition_20260917/).

The local upstream checkout is CENdetectHOR pipeline commit
`969f4ef515fadd8b8a035357294971ff57d4618f` with Snakemake 8.14.0.
The pinned source `Snakefile` SHA-256 is
`7cd4a1ba66224751b37f316fb7f3a33e7a33c18c087407cb30ce126f524ebfe1`.
We preregistered two conditions before running: the archived actual YAML
defaults, `windowSize=1000, kmerSize=10`, and one no-prior diagnostic with
`windowSize=5000, kmerSize=8`. Upstream `config/README.md` lists 8 as its
documented k-mer size; the 5-kb window had already been tried in the earlier
exploratory run. Both use one requested Snakemake core, `consFile=no`, and
`expMonSize` false from the YAML.

An initial invocation passed `expMonSize=no` explicitly. Snakemake treated this
CLI value as a truthy string, so the upstream `int()` cast failed during
Snakefile parsing, before any rule. That failed manifest and receipt remain in
the archive. The addendum records the compatibility correction: omit the
redundant override so YAML `no` parses false; rerun each predeclared setting
once in a fresh workdir. No biological setting or input was changed.

| Native condition | Last completed step | Failure | Wall seconds | Sampled process-tree RSS MiB |
| --- | --- | --- | ---: | ---: |
| Default replay, 1 kb / 10-mer | Periodicity BED, empty | `windowsFiltering.py` indexed empty list | 3.623 | 130.2 |
| Documented 5 kb / 8-mer | Filtered BED `[0,36371)` with nominated period 1710; empty consensus FASTA | `Extract5mon.py` referenced undefined `cons` | 9.510 | 214.0 |

The native `periodicityScript_full_fin.py` evaluates k-mer positions at fixed
10-bp offsets. The 171-bp designed phase does not align to that sampling grid
until 1710 bp; this explains why 1710 is a plausible alias, but does not prove
that it is the only cause of failure. The native no-consensus extractor uses
the nominated period as `monLen` and requires more than about 70% of expected
adjacent monomer intervals before writing a consensus. The 5-kb / 8-mer run
therefore reached an empty consensus file, then the downstream script failed
while trying to read it. Both runs lack the declared final decomposition and
HOR tree. We did not patch upstream period detection, inject the true period,
or adjust the control sequence.

**Endpoint:** CENdetectHOR accuracy on this flanked 171-bp control remains
**technical N/A**. The filtered BED in the diagnostic run is an intermediate
periodicity window, not a completed array or HOR call; scoring it as a final
prediction would change the endpoint. These runs establish a reproducible
precondition failure under the tested native no-prior settings, not a
biological false negative or general failure on plant centromeres. The earlier
5,700-bp shared smoke remains installation/interface evidence only. Resource
figures describe the failed complete-stage attempts, sampled at 0.2 s; they
are not successful-tool performance comparisons.
