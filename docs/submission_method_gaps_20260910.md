# Submission method-evidence gaps — 2026-09-10

## Scope and status

This is a readiness audit, not a prediction of editorial outcome or a claim that
any named journal will accept the work.  It separates what the current evidence
can support from statements that it cannot support.  It does not authorize a
new algorithm, comparator, simulation, or real-data computation.

The frozen production baseline is `81827c3`.  The formal SRF unified comparison is now complete: 72/72 cells, zero process
failures, with three native k151 no_catalogue results retained. Protocol commit
7abafc3 preceded generation. See `srf_formal_results_20260910.md` and the compact
`paper/evidence/srf_formal_unified_v1` archive. The historical-evidence matrix
remains authoritative for what has and has not been run on Ey15 and Macadamia.

## Evidence that is already usable

The existing record supports a reproducible, bounded software and evidence
story:

* TandemX has donor-matched Ey15 and Macadamia discovery, read-abundance and
  assembly-representation products, with confidence and source-eligibility
  states preserved.  The historical primary denominators contain 19 and 43
  source-eligible families, respectively.
* Independent read evidence supports residual under-representation for two
  selected Ey15 families and one selected Macadamia family; three selected
  candidates remain unresolved.  This supports family-specific, qualified
  read--assembly interpretation.
* The offline report retains source hashes, commands, output states and
  candidate relationship warnings.  It is useful as an evidence/provenance
  artifact, not as new biological validation.
* The guarded SRF workflow has observed `ok`, `no_catalogue`,
  `no_eligible_kmers`, and failure fates.  Its prior runs are development
  evidence, not a final common-endpoint comparison.

These statements are documented in
[`submission_workflow_evidence_matrix_20260910.md`](submission_workflow_evidence_matrix_20260910.md),
[`release_program.md`](release_program.md), and the frozen formal protocol.

## Blocks on specific scientific claims

| Claim that must not be made yet | Concrete evidence gap | What the current record can say | What resolves the gap without changing the algorithm |
| --- | --- | --- | --- |
| TandemX is superior to SRF or ordinary competitive mapping on historical real reads | No matched SRF run exists for Ey15 or Macadamia.  Ey15 has no family-level competitive occupancy run; Macadamia has only a direction-only TandemX-template ONT occupancy endpoint. | Historical evidence supports TandemX-specific family prioritization and bounded orthogonal interpretation. | A complete, predeclared matched historical comparison would be required for this claim.  The fresh planted-read study, even if complete, is not a historical real-read comparison. |
| TandemX measures physical repeat copy number or the exact number of missing assembly bases | Read `estimated_bp`, SRF retained mapped bp, and occupancy are different abundance-like quantities.  Newer assemblies are reference proxies; biological replication and physical-copy truth are unavailable. | Report estimated read abundance, assembly representation, and read--assembly abundance deficit/possible under-representation. | Independent physical or otherwise calibrated copy-number truth would be needed; no software change can turn the present proxies into that truth. |
| Population-level collapse sensitivity/specificity, or broad species generality | The historical primary denominators are small (19 Ey15 and 43 Macadamia source-eligible families), and only three of six selected deficit candidates have orthogonal support. | Report the exact denominators, per-family states, and unresolved candidates.  Do not promote a direction-level result to a population estimate. | Additional independent biological material and a preregistered denominator would be needed, not a threshold change. |
| Complete higher-order-repeat reconstruction, validated HOR order, or HOR-aware superiority | `family_hierarchy.tsv` records candidate period multiples and unresolved relations only.  The current formal design permits pure integer-repeat correspondence but explicitly does not validate complex HOR structure. | Say that the tool preserves candidate relationships and does not validate HOR order or reconstruct complete arrays. | A separately validated family/HOR correspondence and structural truth tier are needed before such claims. |
| A new context/phase method is an accuracy improvement over ordinary mapping | The A3 development candidate improved an indel-recall failure relative to A2 but did not show sufficient overall advantage over A1: its final development MARE remained higher and it cost more.  No held-out benefit exists. | Describe this as an investigated development candidate without sufficient overall advantage; do not present it as the released method's novelty. | No additional method work is authorized in this release phase.  Any future novelty claim would require a separately governed, matched validation. |

The last row is a claim boundary, rather than an invitation to restart method
development.  The feature freeze remains in force.

## Limits of the fresh formal comparator

The completed formal study improves the controlled comparator record: it has three seeds, six declared
conditions, hidden truth for post-run scoring, native SRF abundance output and
retained empty/failure fates.  Its design is still deliberately narrow:

* Each dataset has 100 approximately 5-kb reads, three unrelated families, one
  array per positive read, and 12 units per array (70% positive except the 10%
  low-abundance condition).
* The shared-fragment condition is a constructed partial-homology stress test,
  not a genome-wide background-frequency estimate.
* It does not model complex multi-array read architectures, related/mosaic
  families, long or nested HOR structures, array-length heterogeneity, or
  calibrated physical copy number.
* Ordinary mapping receives the TandemX-discovered catalogue.  Its de novo
  family recovery endpoint is therefore `N/A`, rather than a loss or a failed
  discovery result.

Accordingly, even a completed formal study can support controlled
method-behavior statements within these inputs.  It cannot establish general
real-genome superiority, physical-copy accuracy, or complete HOR handling.

## Manuscript figures and package: what remains necessary

The following are package/evidence requirements.  They do not require a new
algorithm, but omitting them would make the evidence boundary hard to assess.

| Package item | Current position | Needed before using it as a central manuscript figure or table |
| --- | --- | --- |
| Formal comparator figure/table | Protocol and runner are frozen; the actual comparison is ongoing. | Preserve all planned cells, including failures, empty outputs and `not_run` states; show both SRF k settings, mapping `N/A` discovery, correspondence ambiguity/unassigned mass, per-family abundance errors, and resource scope.  Do not select a winning SRF setting after inspection. |
| Real-data deficit/orthogonal panel | Three supported and three unresolved selected candidates are documented. | Give every panel its source-data table, exact family denominator, assay/platform, source-eligibility state, and statement that newer assemblies are proxies.  Separate binary interpretation from deficit magnitude. |
| Architecture/HOR panel | Candidate hierarchy and report graphics exist. | Label edges as candidate period multiples or unresolved; do not depict a rooted hierarchy, validated HOR order, or complete array reconstruction. |
| Comparator and limitations table | Historical matrix already distinguishes `not_run` from `N/A`. | Carry those states into the supplement and methods.  In particular, do not fill historical SRF or Ey15 competitive-mapping cells with zeros, nor call the Macadamia ONT direction audit a three-tool comparison. |
| Reproducibility/source-data package | Compact provenance/report products exist. | Include the frozen commit, commands/configuration, input and output hashes, tool build provenance, raw or compact permitted source data, and retained failed/empty receipts.  Verify cross-references, bibliography/metadata, figure legends and table field definitions as a package. |

## Limitations that should remain explicit, not repaired by scope creep

These are limitations to state in the abstract, results, discussion and figure
legends where relevant.  They are not defects that can be solved by adding a
new classifier, threshold or algorithm during feature freeze.

1. Read-first abundance is an estimate with normalization assumptions; assembly
   representation is not physical copy truth.
2. Orthogonal evidence has a family-specific direction/abundance role and does
   not provide biological replication or an exact whole-array ground truth.
3. Candidate period-multiple links are not HOR order, ancestry, or complete
   repeat-array reconstruction.
4. The controlled comparator has small, simplified planted architectures; its
   results must remain scoped to that stress-test distribution.
5. Real historical cross-tool agreement or superiority is unavailable until a
   matched comparison actually exists.  `not_run` and `N/A` are informative
   evidence states, not missing values to impute.

## Readiness conclusion

The present package can be developed into a careful research-software and
evidence manuscript if it is framed around traceable read-derived families,
estimated abundance, assembly representation, and bounded orthogonal support.
It is not ready to claim exact physical copy number, validated HOR
reconstruction, broad population performance, or superiority over SRF/ordinary
mapping on the real historical inputs.  Completing the already-authorized
formal study can close only its controlled-comparator gap; it does not erase
the real-data, denominator, physical-truth, or complex-architecture limits.

## Observed formal comparison changes the claim boundary

Ordinary competitive mapping has lower positive-family MARE in every one of the
six declared conditions. TandemX recovers all planted families but emits 66–82
native families in shared-fragment backgrounds and predicts 27,096–35,717 bp on
negative reads (SRF k151: 0 bp). SRF k101 improves divergence recovery over k151.
These results rule out an unqualified algorithmic accuracy advantage. The real
submission question is whether the integrated evidence workflow adds demonstrated
value to the family audit beyond existing SRF/ordinary-mapping analysis. That
real-input comparative value remains unmeasured. No algorithm change is inferred
or authorized by this observation.
