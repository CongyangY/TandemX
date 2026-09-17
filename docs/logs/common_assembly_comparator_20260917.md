# Equal-input assembly comparator follow-up (2026-09-17)

## Common endpoint that actually completed

We replayed the already archived deterministic 5,700-bp `((A+B)×9+A)×10`
assembly control. A and B are designed 30-bp units; every input base is in the
engineered tandem array. The exact input SHA-256 is
`0acf542594ff333c69ffa5bcb4850c8b3a75fc4c0a828a4d6b6765a7215616c5`.
Both tools received byte-identical FASTA and **no supplied monomer or read
catalogue**. This is a development interface control with no flanking negative
sequence and no independent biological truth. It provides one shared assembly
array-localization endpoint, not a performance or plant-accuracy benchmark.

The replay used the already pinned CENdetectHOR pipeline
`969f4ef515fadd8b8a035357294971ff57d4618f` and library
`bdc3a9df5b60adafdd17eaded1b994edd22aae03`, plus TideCluster source
`8131547e75f73a1e61298a89ff4e861fbd613b7e` and the source-built
TideHunter 1.4.3. CENdetectHOR used its complete 11-job Snakemake pipeline
with `windowSize=1000`, `kmerSize=10`, `consFile=no`, `expMonSize=no`, one
requested core. TideCluster used its `tidehunter` and `clustering` stages with
one requested core, TideHunter period scope 30–100 bp and minimum cluster
array length 100 bp. The stage commands, original executable paths, workdirs,
stdout/stderr and sampled resources are archived under
[`common_assembly_5700_v1`](../../benchmarks/evidence/common_assembly_5700_v1/).
This replay did not add TAREAN/KITE to TideCluster or use HiCAT-human, which
requires an eligible human alpha-satellite input.

| Native output / endpoint | CENdetectHOR | TideHunter → TideCluster |
| --- | --- | --- |
| Stage completion | 11/11 Snakemake jobs, exit 0 | Both stages exit 0 |
| Selected array | BED `[0,5700)` | GFF3 `1–5700`, normalized `[0,5700)` |
| Engineered array base recall / precision | 1.0 / 1.0 | 1.0 / 1.0 |
| Reported repeat unit | One 60-bp A+B monomer consensus | TideHunter 60-bp repeat consensus |
| Designed separate 30-bp A/B units and AB order | Not recovered as designed | Unit labels/order unavailable from these stages |

The two tools therefore agree on a whole-input tandem array and a 60-bp
repeat period. Because this input is all array, base precision has no negative
flank challenge. A 60-bp consensus cannot be credited as recovery of separate
A/B monomers or HOR label order. The `assessment.json` file checks the exact
input hashes, successful stages and native coordinates before recording these
numbers. Native BED, GFF3, monomer FASTA, decomposition and HOR tree were copied
from the fresh run; `SHA256SUMS.txt` covers the compact archive.

The successful process-table-enabled replay took 9.329 s for the complete
CENdetectHOR pipeline, 0.248 s for TideHunter and 11.722 s for TideCluster
clustering. Sampled process-tree RSS peaks were 2,328.5 and 7,861.9 MiB for
CENdetectHOR and clustering. RSS sums live processes and may double-count
shared pages. TideHunter's 0.234-MiB sampled peak missed its short native
process and is **invalid as a memory measurement**. These developmental times
include command startup and cannot support a cross-tool speed claim. A first
sandboxed replay also completed but returned zero process-table RSS for all
stages; its logs and traces are retained under `profiles/` and are not used for
memory inference.

## Why the flanked 171-bp control remains unresolved

The separate, previously frozen 36,371-bp assembly control is the meaningful
array-localization challenge because it contains two 1,000-bp negative flanks.
Its primary CENdetectHOR run produced an empty periodicity BED; upstream
`windowsFiltering.py` then indexed an empty list and exited. An exploratory
5-kb-window no-prior run nominated a 1,710-bp periodic interval but emitted an
empty consensus FASTA; the following `Extract5mon.py` raised `NameError` while
reading it. A further truth-length (`expMonSize=171`) exploratory condition
also stopped before final predictions. The upstream periodicity implementation
searches k-mer positions at 10-bp steps, which is relevant to diagnosing its
1,710-bp nomination on a 171-bp synthetic design; this is a code-level
inference, not a demonstrated root-cause repair. The header-format compatibility
issue was already fixed before the frozen primary run and both tools then used
the corrected identical FASTA. We did not change CENdetectHOR's period/window
algorithm or keep searching parameters against known truth.

Thus CENdetectHOR's 171-bp result remains a **technical incomplete run**, never
a false negative or zero accuracy. TideCluster localized that flanked array in
the earlier frozen run; its 855-bp reported period is a canonical HOR-length
signal, not a five-label HOR-order call. TandemX `locate` on that input received
truth-derived monomers, whereas these no-prior tools did not. No common
same-prior TandemX/CENdetectHOR/TideCluster HOR endpoint is available from the
current methods. The pinned CENdetectHOR source tree contains README input
instructions but no bundled FASTA/BED example to use as a separate official
completed benchmark. Its previously completed 30-bp full-pipeline toy and this
fresh equal-input replay establish that the installation can finish on an
eligible interface input; neither resolves the flanked 171-bp case.
