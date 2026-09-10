# Formal SRF unified comparison: observed results

Protocol/source frozen at commit 7abafc3 before generation. All 72 cells
completed in 174.336 s controller elapsed; no process failures or unrun cells.
TandemX production and Rust source remain baseline 81827c3. No A3, third
algorithm revision, new detector, 1-Gbase timing or release action was performed.

## Existing development evidence

The four pilots, original 32-call development suite and guarded 32-call rerun
remain development evidence, separately inventoried in
`srf_development_inventory_20260910.md`. None was repeated for this comparison.
The new dataset seeds and protocol are in `srf_formal_unified_protocol_20260910.md`.
The A3/history FPR definitions differ, and the absolute .02 gate is withdrawn;
see `fpr_scope_reconciliation_20260910.md`.

## Formal abundance comparison

MARE below is the arithmetic mean across three new simulated datasets per
condition, with all three positive-truth families per dataset included (nine
family rows per condition). A missed true family contributes zero estimated bp.
No truth family had zero observations in the realized datasets. These are
sampled-read abundance errors with depth=1, not calibrated genome copy errors.

| Condition | TandemX | SRF k=151 | SRF k=101 | Competitive mapping |
| --- | ---: | ---: | ---: | ---: |
| Clean | 0.008361 | 0.045623 | 0.049840 | 0.000359 |
| 1% substitutions | 0.193880 | 0.057697 | 0.052819 | 0.000261 |
| 0.1% insertion + 0.1% deletion | 0.047374 | 0.055464 | 0.049672 | 0.000251 |
| 2% unit divergence | 0.349075 | 1.000000 | 0.162351 | 0.000127 |
| 10% positive reads | 0.003704 | 0.035118 | 0.047998 | 0.000394 |
| Shared-fragment background | 0.013280 | 0.058794 | 0.046549 | 0.000437 |

MARE alone does not measure false catalogue discovery. Unmatched native abundance
is explicitly separate and must not be omitted when interpreting this table.
Ordinary mapping uses TandemX's discovered catalogue; its de novo recovery is N/A.

## Recovery, background attribution and retained mass

TandemX recovered 54/54 positive truth-family conditions. SRF k151 recovered 45/54:
all nine missing family conditions were in the three 2% divergence datasets,
where native SRF exited successfully with `no_catalogue`. SRF k101 recovered
53/54, including 8/9 divergence family conditions. No result-dependent choice
between SRF k values is permitted. SRF's clean/indel interval recall is measured
on its native filtered intervals, not its unfiltered mapping output.

In the three shared-fragment datasets:

| Endpoint (mean per dataset) | TandemX | SRF k151 | SRF k101 | Competitive mapping |
| --- | ---: | ---: | ---: | ---: |
| Negative-read predicted bp | 30,886.67 | 0 | 227.33 | 11,165.67 |
| Negative-read predicted fraction | 0.205911 | 0 | 0.001516 | 0.074438 |
| Unassigned native abundance bp | 10,528 | 0 | 0 | 11,008 |
| Native catalogue size | 71.33 | 3 | 3 | 71.33 (shared catalogue) |

TandemX native catalogue sizes were 66/82/66; global base-union precision was
0.770218/0.737413/0.786916. Competitive mapping retained the same catalogues and
its precision was 0.901546/0.837550/0.970463. In particular, very low MARE for
known positive families does not mean this baseline has no false attribution.
The background is deliberately constructed partial homology, not population
specificity or an assertion that every such motif lacks biological interest.
Discovery-interval negative bp must not be described as k-mer-quantifier
localization. Full native quantities and family correspondences are retained.

## Resource measurements

`condition_summary.tsv` reports wall time, user/system CPU, peak child RSS,
throughput and observed ranges. Three independent seeds are not technical
repeats. As a check on accounting, clean mean wall times are TandemX 0.572 s,
SRF k151 2.233 s, SRF k101 1.895 s, and mapping-only 2.060 s. The mapping table
also provides shared-discovery-plus-mapping costs; the discovery is not free.
No GB-scale speedup or memory-reduction claim follows from these ~0.5-Mb inputs.
Peak temporary disk was not measured. Retained output byte counts are separate.

## Interpretation and remaining gap

The comparison does not support a general accuracy advantage of the production
k-mer estimator over ordinary competitive occupancy. It also exposes substantial
partial-homology catalogue/interval over-attribution and substitution/divergence
sensitivity in the tested uncorrected FASTA path. Conversely, the configured
TandemX detector retains family recovery where SRF k151 has no catalogue, and
SRF k101 substantially reduces that gap. These results support conditional
trade-offs and an integrated family-auditing workflow claim, not a uniquely new
measurement or universally best repeat detector.

The formal **historical real-assembly prioritization** SRF/ordinary-mapping
comparison remains `not_run`. Existing Ey15/Macadamia validation belongs to
TandemX, with a newer-assembly proxy and three orthogonally supported selected
families; it does not establish added real-data value over these comparators.
That remains the main comparative methodological gap for a Genome Research
submission. New algorithm development is not reopened. The manuscript also
requires its final figure/source-data package and later formal scaling study.

## Evidence

- Compact, hash-manifested evidence: `paper/evidence/srf_formal_unified_v1/`.
- Full raw/native/source archive: `/Volumes/T7/Codex/TandemX/results/srf_formal_unified_v1_20260910`.
- Independent arithmetic/native-abundance audit: `srf_formal_independent_audit_20260910.md`.
- Workflow capabilities and missing endpoints: `submission_workflow_evidence_matrix_20260910.md`.
- Revised manuscript: `paper/0910/manuscript_v4.md`.
