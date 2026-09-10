# Formal SRF unified-endpoint protocol, v1

Status: design frozen before generation of the new comparator datasets. Execution
starts only after a committed protocol and tested runner are recorded in the run
receipt. This is a final comparator experiment for the frozen production method,
not an A3 held-out evaluation. No parameter tuning or algorithm changes follow
inspection of these results. A3 and TXF000695 are excluded.

## Evidence and scope

Existing four SRF workflow pilots and 32 development calls, including the guarded
empty-k-mer rerun, are inventoried in `srf_development_inventory_20260910.md`.
They remain development evidence; they will not be rerun or relabeled final.
Historical Ey15/Macadamia workflow evidence is separately inventoried in
`submission_workflow_evidence_matrix_20260910.md`. This new experiment evaluates
controlled planted read arrays, not physical genomic copy number or historical
assembly prioritization. Applicable historical SRF and competitive-mapping
comparisons remain `not_run`; unsupported endpoints are `N/A`.

## Inputs, before outcome inspection

Three seeds: 2026091001, 2026091002, 2026091003. Each has the same six declared
conditions below: 18 datasets. Each dataset has 100 reads of approximately
5,000 bp, three unrelated families, 12 tandem units per positive read, one array
per positive read, and 70% positive reads unless specified. Observed coordinates
after insertion/deletion are truth. Every truth family is in the abundance-error
denominator, including families not recovered by a method. A zero-truth family
would be reported separately, not divided by zero. Recovery also uses only truth families with at least one
planted array; the potential three-family catalogue and zero-truth IDs/counts
remain explicit.

| Condition | Change from clean 171-bp units |
| --- | --- |
| clean | None |
| substitution | 1% independent per-source-base substitution |
| indel | 0.1% insertion AND 0.1% deletion (0.2% combined nominal events) |
| divergence | 2% independent substitutions per tandem unit |
| low_abundance | 10% positive reads |
| shared_fragment | 120-bp units; negative reads contain first 100 bp of family 1, separated by independent random 30–80-bp gaps after a 200-bp random prefix |

Shared-fragment negatives have no planted full tandem array. Their partial
homology is a deliberately constructed challenge, not an estimate of genome-wide
background frequency. The historical simulator source remains byte-identical. A benchmark-only
negative-read overlay uses Random(seed*1000+zero-based-read-index), preserving
positive arrays and recorded strand while replacing random negative sequences.
Neutral read IDs and input FASTA are the only read inputs to native tools; truth
files are used solely for post-run scoring. Manifest SHA-256 binds all inputs.

## Frozen native methods

Production TandemX scientific source is frozen at
`81827c3e3fcd04bddbeba26c8165ac680f1b232e`; the comparator commit must verify
`tandemx/` is unchanged. Use the existing explicit analysis configuration:
`discover --discovery-method cascade --clustering-method sequence`, periods
30–1000 bp, minimum span 100 bp, support 1 read, minimum read length 1 bp,
Rust k-mer backend, one thread. This is an explicitly configured workflow,
not a claim about out-of-box defaults. Quantify its own discovered catalogue
with k=21, Rust backend and haploid depth=1. Quantify has no thread CLI flag;
its single input file uses the serial selected-counter path (the Rust batch
loops sequentially), rather than the multiple-file worker path. Genome size is total sampled bases
(required but inactive for explicit-depth normalization). Native `estimated_bp`
therefore estimates repeat bases in this sampled read collection.

The formal receipt binds the clean, consistent KMC rebuild and canonical-count
verification receipt, native binary hashes, and the previous upstream source/build
provenance. The earlier mixed-standard incremental KMC binary is excluded.

SRF primary: KMC k=151, minimum count=20; sensitivity: k=101, count=20. Both
are always reported; no selection of the better result. Use the existing native
seven-stage KMC/dump/SRF/enlong/minimap2/paf2bed/bed2abun workflow and its native
elongation and filtering defaults, one thread, KMC 2-GB cap. Preserve empty
eligible-k-mer and empty-catalogue biological/output states separately from
process failure. Native retained bp from bed2abun is the abundance endpoint.

Ordinary competitive-mapping occupancy baseline: use the **TandemX-discovered**
catalogue, not truth sequences. De novo recovery is N/A (shared catalogue).
Repeat each native unit to a 10,000-bp template. Run minimap2 `-x map-hifi -c
-N1000000 -f1000 -r100,100 -t1`; retain alignments of block length >=100 bp
and nmatch/block >=0.9, including secondary alignments. Union within native
family; bases mapped to multiple native families are ambiguous and excluded
from family abundance before any truth correspondence. All qualifying intervals,
including ambiguous ones, remain in global interval precision. Report ambiguous
mass. Report mapping-only cost and shared-discovery-plus-mapping cost separately.

## Common scoring and limits

Independent edlib global edit correspondence, >=0.90 identity, considers circular
rotations/reverse complements and pure integer repeats of a unit up to 64 times
(within 10% expected length). Exactly one passing truth family is a unique match;
multiple passing families are ambiguous, with no best-match forced assignment.
Mixed/unmatched/out-of-scope motifs retain their native quantities and separate
mass. No abundance is divided by the number of units in an SRF HOR: native bp is
already a base amount. This correspondence does not validate complex HOR structure.

Report recovery as unique truth families recovered / all positive-truth planted families for
independent discovery methods. Native abundance mapped to the same truth family
is summed, with unknown mass retained separately. Missing families receive zero
estimated bp. MARE is mean(abs(estimate-truth)/truth) over positive-truth families;
report each family as well as dataset means. Excess estimated bp is separate
from localization false attribution and is not labeled a binary FPR.

Global base-union recall/precision scores all native intervals, including
unmatched and ambiguous predictions; same-base duplicate calls are counted once.
Background false attribution is predicted union bp on reads with no planted
array, divided by total negative-read bp when reporting a fraction. Family-level
interval scores withhold conflicting family assignments. TandemX discovery
intervals do not directly measure where its k-mer quantifier allocated abundance;
its abundance excess is therefore a distinct endpoint. None of these quantities
is directly comparable to historical held-out binary FPR=16/972.

Every dataset runs once per method/configuration. Three seeds are independent
simulations, not technical timing repetitions. Record stage and end-to-end wall
time, user/system CPU time, maximum measured child RSS, throughput, output hashes,
and retained output bytes. Retained output bytes are not peak temporary disk;
peak temporary disk is `not_measured`. Run native processes serially, record host
context, and do not claim large-scale speedup or RAM reduction from these results.
The deferred old/new >=3-repeat ~100-Mb/~1-Gb scaling study is not part of this run.

Failures/timeouts retain logs and fates, not zero scores. Empty valid native output
can score zero recall/abundance with precision undefined. Stage timeout 180 seconds;
stop launching work after 45 minutes, retaining unrun cells as `not_run` and
allowing an already started bounded stage to exit. Resumption requires matching
input/source/protocol hashes; no overwriting completed attempts. Mechanical
adapter corrections must preserve original attempts and record the correction.

TideHunter, TRF, TRASH/TRASH2 and TideCluster existing results may be cited in a
separate matched-endpoint evidence inventory. They are not entered into this new
18-dataset table without actually running the same inputs, and are never scored
as failing an endpoint outside their design. Real historical-assembly comparisons
and manuscript artwork remain separate submission readiness checks.
