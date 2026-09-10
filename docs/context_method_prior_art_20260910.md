# Prior-art audit: periodic-context family abundance

## Question answered

This audit asks whether a proposed **context-constrained family abundance**
estimator is meaningfully different from known seed-and-chain, minimizer, or
satellite-reconstruction methods. It is a primary-source/official-code review,
not an implementation proposal or novelty claim.

The proposed idea is only potentially distinct if it is defined as an
**estimator eligibility rule**, not as “use seeds and chains”: for a pre-existing
family and its already-defined diagnostic anchors, count a read interval toward
family abundance only if it contains repeated, orientation-consistent family
anchors in collinear order with inter-anchor distances compatible with the
family period (within a predeclared tolerance). Isolated or non-periodic
anchor hits are retained as rejected evidence, not counted. The denominator and
normalization remain the existing abundance contract.

If the proposal merely finds exact seeds, chains them, or calls a tandem period,
it is established prior art and is not a new algorithm.

## Closest verified precedents

| Source | Concrete existing mechanism | What it does **not** establish about the proposed estimator |
| --- | --- | --- |
| [TideHunter paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6612900/) and [official code](https://github.com/yangao07/TideHunter) | A seed-and-chain algorithm recognizes a repeat pattern within an individual noisy/RCA long read, partitions it into repeat units, and applies POA to form a consensus. | This is read-local tandem detection/consensus calling. The paper/code does not define a cross-read, fixed-family abundance estimator that gates diagnostic-anchor counts by periodic context before normalization. Reusing its seed-chain idea alone is therefore not novel. |
| [minimap2 paper/source methods](https://github.com/lh3/minimap2/blob/master/tex/minimap2.tex) and [official code](https://github.com/lh3/minimap2) | Query minimizers seed exact reference anchors; dynamic programming forms collinear chains with a gap cost, then alignment extends chains. Frequent seeds may be skipped. | Collinearity is a standard mapping primitive, and its gap cost is not a repeat-period eligibility rule. It maps a query to a reference; it neither defines a family diagnostic set nor estimates family abundance from periodic support in unassembled reads. |
| [SRF paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760446/) and [official workflow](https://github.com/lh3/srf) | High-occurrence k-mers form a de Bruijn graph; SRF greedily follows the highest-occurrence successor at branches to reconstruct cycles. For abundance, elongated reconstructed units are mapped back to source reads/contigs, then competing mappings are resolved heuristically before summing retained mapped bases. | SRF reconstructs motifs/HORs and estimates mapped-bp abundance after contig mapping. Its documented abundance step does not require repeated preassigned family-diagnostic anchors at period-consistent spacing within every accepted source interval. Conversely, the proposed rule must not claim to reconstruct HORs or replace SRF's graph reconstruction. |
| [Winnowmap weighted-minimizer paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7355284/) and [official code](https://github.com/marbl/Winnowmap) | Weighted minimizer sampling down-weights highly repetitive k-mers to reduce excessive false matches while retaining a window guarantee for mapping; the method targets long-read mapping in repetitive sequence. | Frequency weighting/sampling is already known and is not family assignment. A context estimator that merely reweights or drops repetitive seeds is not distinguishable. Its testable contribution would have to be periodic, multi-anchor family support for *counting*, including explicit accepted/rejected interval evidence. |
| [syncmer conservation paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7869670/) and [reference implementation](https://github.com/rcedgar/syncmer) | Syncmers are a seed-sampling scheme intended to retain more conserved k-mers under mutation than minimizers. | Choosing syncmers instead of minimizers changes seed sampling, not the biological claim of family-specific periodic support or the abundance estimator. It is an implementation option only and cannot itself be presented as novelty. |

## Falsifiable distinction from prior art

The candidate method has a defensible separate contribution only if all of the
following are specified before evaluation:

1. **Fixed inputs.** Family representatives and diagnostic-anchor ownership are
   fixed by the current frozen catalogue; no truth, post hoc assembly sequence,
   or read-context result may change membership.
2. **Counting predicate.** An accepted interval requires at least a declared
   number of owned anchors, one declared orientation rule, collinear ordering,
   a declared supported span, and anchor-spacing residuals compatible with an
   integer multiple of that family's period. The predicate must say how indels,
   reverse complements, overlapping anchors, and ambiguous/shared anchors are
   handled.
3. **Estimator difference.** The output is abundance based on accepted
   context-supported bases (plus a separately reported rejected-anchor burden),
   not a new tandem detector, mapper, or graph assembler.
4. **Discriminating tests.** On the same fixed reads, it must reject decoy
   contexts containing isolated family anchors or non-periodic dispersed copies
   while retaining genuinely periodic arrays; then it must improve or expose a
   trade-off relative to the existing unconstrained diagnostic-k-mer estimate.
   A seed-chain implementation that produces the same accepted counts, or a
   gain only from retuning family membership, falsifies a distinct-method claim.
5. **Independent boundary.** The method must not describe a context-supported
   interval as a validated HOR, genomic locus, or physical copy number. It
   estimates family-associated read occupancy under its stated rule.

## What would be new versus what is already known

| Candidate component | Prior-art status | Required evidence before calling it a contribution |
| --- | --- | --- |
| Exact k-mer/minimizer/syncmer seeds | Established. | None; cite the seed scheme. |
| Collinear seed chaining | Established (minimap2 and many aligners). | None; do not describe as a new algorithm. |
| Read-local tandem period/consensus detection | Established (including TideHunter). | None; do not relabel it as a new detector. |
| High-occurrence-k-mer satellite graph reconstruction | Established (SRF). | None; do not claim a replacement. |
| Weighted suppression of repetitive anchors | Established (Winnowmap). | Demonstrate a different objective if used. |
| Fixed-family, periodic-context eligibility applied before abundance counting, with accepted and rejected occupancy reported | Not identified in these closest sources as their stated abundance endpoint. | A written rule, an independent implementation/audit, held-out decoy versus periodic tests, and comparison to the unchanged unconstrained estimator. This is a hypothesis, not a novelty claim yet. |

## Decision boundary

Proceed only if the proposed method is evaluated as a **bounded abundance
estimator** and beats, or clearly calibrates, the existing estimator on data not
used to define its predicate. Otherwise retain it as an explanatory diagnostic.
Do not open a new detector, HOR reconstruction, minimizer, syncmer, or generic
seed-chain development line: each is covered by established methods above.
