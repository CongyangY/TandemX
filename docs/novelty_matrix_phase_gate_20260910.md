# Phase-gate novelty matrix: family-conditioned periodic evidence

## Scope

This is a mechanical primary-source/official-implementation audit for the
accuracy-branch candidate described as family-conditioned, monomer-phase-aware,
indel-tolerant evidence scoring. It makes no conclusion that the candidate is
novel or should enter TandemX. “Not found” below means only that an equivalent
endpoint was not located among the selected primary sources; it is not evidence
of global absence.

The audit separates (a) seed construction, occurrence-context hashing,
mapping/chaining, and read-local repeat detection from (b) an abundance score
conditioned on a pre-existing repeat family.

## Verified prior-art matrix

| Concept in the gate | Primary source / official implementation | Verified mechanism | Relation to family-conditioned abundance scoring |
| --- | --- | --- | --- |
| Interval hashing / `mm2-ivh` | [Suzuki *et al.* 2025](https://doi.org/10.1093/bioinformatics/btaf648); [official `mm2-ivh`](https://github.com/ocxtal/mm2-ivh) | For an identical k-mer, encodes intervals to `W` same-k-mer occurrences on each side into an interval hash; the k-mer plus hash distinguishes repeat contexts for high-identity long-read overlap/mapping. The implementation states effectiveness requires roughly >=99.5% identity. | This is genuine multi-occurrence local context and cannot be renamed as new. It is an occurrence-context key for mapping, not an estimator that accepts/rejects evidence for a fixed biological family by monomer phase or reports family abundance. Any proposed method must compare directly to an interval-hash context baseline if it uses occurrence spacings as its principal discriminator. |
| Multi-context seeds (MCS) | [Tolstoganov *et al.* 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13059148/); [official strobealign](https://github.com/ksahlin/strobealign) | Strobemer hashes distribute bits over component strobes so a prefix lookup can query full or partial multi-strobe contexts at several resolutions. It is a general read-mapping seed/index scheme. | Multiple context scales, partial context lookup, and strobemer construction are prior art. The paper does not define a repeat-family abundance score, periodic monomer-phase residual, or a counted-versus-rejected family occupancy output. |
| Syncmers | [Edgar 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7869670/); [reference code](https://github.com/rcedgar/syncmer) | Selects k-mers by an internal s-mer rule to improve conservation under mutation relative to minimizers. | A syncmer changes which seeds survive errors. It neither establishes repeated copies/phase coherence nor assigns abundance to a repeat family. Calling a syncmer substitution phase-aware or context-aware abundance would be incorrect. |
| Strobemers | [Sahlin 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8554264/) | Links multiple k-mers spaced within a window, improving sensitivity under mutations/indels for sequence matching. | Linked seeds and error-tolerant spaced context are established. A new claim cannot rest on combining two or more seeds; it would need a different, family-conditioned scoring target. |
| Weighted minimizers | [Jain *et al.* 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7355284/); [official Winnowmap](https://github.com/marbl/Winnowmap) | Down-weights high-frequency k-mers during minimizer selection to improve repetitive-region long-read mapping while maintaining a window guarantee. | Specificity weights for seeds are known. A family-specific weight is not a new method unless its definition and effect are separate from generic repetitive-seed suppression. |
| Collinear seed chaining | [Li 2018](https://github.com/lh3/minimap2/blob/master/tex/minimap2.tex); [official minimap2](https://github.com/lh3/minimap2) | Exact minimizer matches become anchors; dynamic programming chains collinear anchors using a gap cost and then extends them by alignment. | Ordered anchors and gap penalties are standard mapping. A bounded phase-drift DP resembles a chain unless it has a distinct state space/objective tied to known repeat period and fixed family ownership. |
| Period-aware/error-tolerant read-local repeat recognition | [TideHunter](https://pmc.ncbi.nlm.nih.gov/articles/PMC6612900/); [official code](https://github.com/yangao07/TideHunter) | A tailored seed-and-chain step recognizes a repeat pattern in noisy/RCA long reads, partitions units, then SIMD POA produces a consensus; reported tolerance reaches high read error rates. | Period recognition, seed chains, unit partitioning and error tolerance are established. TideHunter is not a fixed-family, cross-read abundance estimator, so its mechanism is not automatically an equivalent abundance endpoint. |
| Error-tolerant periodic matching | [EquiRep](https://pmc.ncbi.nlm.nih.gov/articles/PMC11580891/) | Uses seed chaining for coarse candidate selection, then diagonal-free self-alignment to infer a repeating unit from error-prone sequence. | Self-alignment and indel-tolerant periodic detection are prior art. A local phase reset without full alignment may be an implementation trade-off, not novelty, unless it is shown to create a distinct family-conditioned abundance score and outcome. |
| Circular repeat/HOR reconstruction and mapped-bp abundance | [SRF](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760446/); [official workflow](https://github.com/lh3/srf) | Builds high-occurrence k-mer graph cycles, greedily follows high-occurrence branches, reconstructs satellite units/HORs, then maps elongated units back to reads/contigs and resolves competing mappings before reporting abundance. | Circular repeat representation, HORs and mapping-based abundance are prior art. SRF does not, in its stated workflow, score each preassigned family seed by a local monomer-phase residual or report a rejected non-periodic diagnostic-anchor burden. |

## Terminology checks

| Term | Audit result |
| --- | --- |
| “mm2-ivh interval hashing” | Verified: a minimap2 fork for high-identity repeat-aware overlap/mapping; it hashes intervals among repeated occurrences of the same k-mer. |
| “multi-context seeds” | Verified: a named strobemer-based multi-resolution seed/index construction, not a generic synonym for multiple nearby k-mers. |
| “period-aware seeds” | No single canonical primary method with this exact label was located in the selected sources. Period-aware repeat recognition itself is verified for TideHunter and EquiRep; do not treat the wording as a new category. |
| “phase-aware tandem-repeat analysis” | No selected primary source was found that uses this exact phrase for a fixed-family abundance estimator. Cyclic rotation/phase and repeat-unit ordering are inherent in circular/HOR and tandem-repeat work, so the phrase alone is unverified as a distinct method name and should not be used as a novelty shortcut. |
| “circular/error-tolerant periodic matching” | Circular repeat reconstruction is verified for SRF; error-tolerant periodic matching is verified for TideHunter/EquiRep. No selected source combined them into the specific proposed abundance endpoint. |
| “repeat abundance” | Verified for SRF as mapped retained bases after contig-to-source mapping; this must not be relabelled as a new read-occupancy concept. |

## Mechanical equivalence tests required before naming the method

The candidate should be considered equivalent to known seed/mapping machinery
if its score can be reproduced by one of: (1) generic collinear-chain score,
(2) interval-hash match/no-match, (3) multi-context seed lookup, (4) seed
frequency weighting, or (5) local tandem detector output. To show a distinct
endpoint, the design must predeclare all of the following:

1. A frozen family anchor set, circular monomer phase for each anchor, and
   specificity weight. These cannot be updated from held-out reads.
2. A local score whose state explicitly includes family identity and phase
   offset; define phase residual, copy recurrence, phase coverage, order rule,
   reset/indel penalty, and ambiguous-family output.
3. A counting rule that turns only accepted context-supported spans into the
   abundance numerator and separately reports rejected shared/non-periodic
   hits. This differs materially from merely outputting a chain or an alignment.
4. Discriminating data: shared-fragment background with non-periodic hits,
   true arrays with substitutions/indels, interruptions, related families, and
   heterogeneous members. The same inputs must be scored by unconstrained
   diagnostic k-mers, ordinary competitive occupancy, strict periodic context,
   and the proposed indel-tolerant score.
5. Negative conditions: if phase scoring does not reduce background attribution
   beyond generic occurrence-context/competitive mapping, or if resets simply
   reproduce local chain/alignment occupancy, the distinct-method hypothesis is
   falsified. If it gains only by reassigning family membership, it is not an
   abundance-scoring result.

## Boundary for the phase gate

The reviewed sources establish all generic building blocks: circularity,
periodicity, error tolerance, local occurrence context, linked/multi-context
seeds, seed weighting, and collinear chaining. The only candidate contribution
left for testing is the **combination as a fixed-family abundance eligibility
and uncertainty score**. Whether that combination is algorithmically distinct
depends on the predeclared implementation and held-out results; this audit does
not decide it.
