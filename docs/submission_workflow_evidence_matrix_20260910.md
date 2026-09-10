# Submission workflow and historical-evidence matrix

## Scope

This is a read-only capability and evidence audit for the planned formal SRF
comparison: three fresh seeds × six planted-read conditions (`clean171`, 1%
substitution, 0.1% indel, 2% unit divergence, 10% low abundance, and a 100-bp
shared background for a 120-bp unit). The TandemX arm is fixed to production
baseline `81827c3`, cascade plus sequence, Rust and one thread. SRF reports
k=151/ci20 as primary and k=101/ci20 as sensitivity. No benchmark, real input,
or historical re-analysis was run for this audit.

“Not run” means no matching execution artifact was found. “N/A” means an
endpoint is outside the stated arm, rather than a missing output.

## Historical evidence cannot provide a common tool-prioritization endpoint

Existing TandemX evidence can rank its own discovered families by frozen
read--assembly comparison and retain interpretation states. There is no matched
SRF execution on Ey15 or Macadamia and no SRF/TandemX family-or-HOR
correspondence table. It therefore cannot support a cross-tool historical
prioritization comparison without new execution and a predeclared common score.

| Historical material / endpoint | TandemX already executed | Ordinary mapping already executed | SRF already executed | Common historical endpoint |
| --- | --- | --- | --- | --- |
| Ey15 donor-matched catalogue, abundance and old/new arrays | Yes. `paper/evidence/ey15_donor_matched_collapse_v1/results/run/` retains discovery, primary/depth-107 copy-number and old/new arrays. Primary metrics retain 2,133 rows, 19 source-eligible rows and `not_source_eligible` states. | Not run as a family-level competitive read-occupancy table. The retained PAF is old-versus-new assembly context, not read-to-family occupancy. | Not run on Ey15 reads; available SRF products are synthetic development only. | Not available. |
| Macadamia donor-matched catalogue, abundance and old/new arrays | Yes. `paper/evidence/macadamia_jansenii_donor_matched_collapse_v1/results/run/` retains discovery, primary/reported-depth copy-number and old/new arrays; primary metrics retain 1,227 rows and 43 source-eligible rows. | Partly, for a separate secondary endpoint. `orthogonal_abundance_validation_v1/ont/ont_occupancy.tsv` covers 43 eligible TandemX representatives; the direction table evaluates three preselected candidates only. | Not run on Macadamia reads. | Not available. The ONT output is a TandemX-template, direction-only audit, not a three-tool comparison. |
| Final orthogonal interpretation | Yes. `final_evidence_summary.json` records two Ey15 and one Macadamia family with orthogonal support, and three unresolved candidates. | Macadamia ONT gives direction only; Ey15 has orthogonal Illumina k21/k31 metrics, not mapping occupancy. | Not run. | Not available; do not relabel it as cross-tool validation. |
| De novo family recovery | Historical TandemX discovery outputs exist. | N/A: the planned mapper receives TandemX representatives, so it is an occupancy comparator. | Not run historically. | N/A for ordinary mapping; unmeasured for SRF. |

The reusable mapping implementation is
`benchmarks/scripts/evaluate_macadamia_ont_occupancy.py:summarize_paf`: it
unions accepted primary query intervals by family and excludes whole reads with
multiple accepted families. Its mapping-efficiency correction is explicitly
direction-only, not an SRF abundance estimator or physical copy-number truth.

## Fresh-comparison endpoint readiness

| Endpoint | TandemX arm | SRF arm | Ordinary mapper | Remaining link |
| --- | --- | --- | --- | --- |
| De novo sequence/family recovery | Operational monomer families and representatives. | Native `srf.fa` motifs, possibly HORs. | N/A by design. | Lock family/HOR-aware correspondence before scoring. |
| Read interval/base recovery | Candidate/read evidence needs a validated candidate-to-family crosswalk and interval union. | `srf.paf` to native BED plus normalized predictions are emitted. | Query-interval union has a reusable implementation. | Use one shared scorer and frozen mapping/ambiguity policy. |
| Abundance-like quantity | Frozen diagnostic-k-mer estimated bp with normalization warnings. | `bed2abun` mapped-bp fraction divided by observed FASTA bases. | Library-fraction occupancy after mapping. | No shared calibrated physical-copy endpoint exists. |
| Array/negative controls | Challenge interfaces score hidden truth after execution. | In-scope 30--1000-bp, >=100-bp arrays and empty states are preserved. | Not a de novo recovery arm. | Generate six conditions with truth withheld from commands. |
| Resources/failures | Formal protocol must retain command/resource receipts. | Seven native-stage commands, stderr, hashes, wall/RSS and statuses are already retained. | Needs same receipt discipline. | Existing SRF timings are one-run development observations, not paired repetitions. |

`run_srf_pilot.py:workflow()` accepts and receipts an explicit k; existing
development runs establish k151/ci20 execution. The CLI `run()` still calls the
default k151, and no completed k101/ci20 result was located. k101 is therefore
planned sensitivity output, not prior evidence.

## Workflow capability evidence

| Capability | Evidence | Supported behavior | Boundary |
| --- | --- | --- | --- |
| Operational families and sequence membership | `tandemx/discover/clustering.py`; historical `families.tsv`, `monomers.fa`, `monomer_membership.tsv`, and family-audit products. | Traceable read-derived operational clusters and representatives. | Not ancestral or automatically homologous families. |
| Candidate relationships and unresolved state | `tandemx/discover/hierarchy.py`, `family_hierarchy.tsv`, `tandemx/report_data.py`. | Candidate period-multiple and unresolved relationships can be displayed. | Not validated HOR order or complete array reconstruction. |
| Abundance, assembly and confidence | Frozen `copy_number.tsv`, `arrays.bed`, `assembly_vs_read_cn.tsv`; report preserves stage confidence/warnings. | Estimated read abundance, assembly representation, possible under-representation and unavailable state. | Deficit is not missing-base truth; failed/stale stages are excluded when records are supplied. |
| Orthogonal interpretation | `orthogonal_abundance_validation_v1/final_evidence_summary.json` and candidate tables. | Three supported and three unresolved final states with input hashes. | No automatic upgrade of unresolved candidates; magnitude and binary state remain separate. |
| Offline evidence report/provenance | `report_data.py`, `report_html.py`, `report_figures.py`, `reporting.py`; report and pipeline tests. | Offline HTML, TSV/JSON/GraphML, SVG/PDF/PNG/source/receipt products, source hashes and commands. | Reporting adapter does not re-cluster, re-estimate or mark recovery partial as success. |
| Compact-source report reuse | `benchmarks/scripts/render_existing_family_report.py`. | Hash-verified compact copy, source-reuse manifest and report without reopening raw reads. | Not a new benchmark; recovery needs a complete hash-valid candidate lock. |
| SRF workflow and states | `run_srf_pilot.py`, `docs/srf_workflow.md`, and T7 `srf_workflow_{pilot,development,empty_guard}_v1_20260906`. | KMC/SRF/enlong/minimap2/paf2bed/bed2abun, native product retention, and `ok`, `no_catalogue`, `no_eligible_kmers`, `failed`. | Existing inputs are s1101 development data and motifs may be HORs. |
| SRF observed development | `srf_development_inventory_20260910.md`; guarded 32-run receipt has 17 `ok`, 4 `no_catalogue`, 11 `no_eligible_kmers`, zero process failures. | Operational workflow and honest empty-state handling. | Neither a real-read historical run nor a final TandemX-versus-SRF result. |

## Submission-safe conclusion

Before the authorized fresh comparison, the donor-matched records support only
TandemX-specific prioritization and bounded orthogonal interpretation. They do
not establish historical agreement, disagreement or superiority versus SRF or
ordinary mapping. The formal fresh study should retain every tool status and
raw native product, score the locked common endpoints, keep ordinary mapping
de-novo recovery as N/A, and preserve SRF higher-order motifs rather than force
them into monomer labels.
