# Frozen orthogonal abundance validation

This compact archive is the final method-science gate for the frozen Ey15-2
and *Macadamia jansenii* family catalogues. It contains all family-level
Illumina results, the Macadamia ONT direction audit, input-QC summaries and the
final cross-platform interpretations. Raw reads, KMC databases, mapping PAF
files and scratch data remain under `/Volumes/T7/Codex/TandemX` and are not
versioned.

## Independent evidence

- Ey15-2: PCR-free Illumina run `ERR8666067`, independently extracted from the
  same ground tissue pool as the HMW DNA. The converted paired FASTQ contains
  157,774,340 reads and 23,666,151,000 bases.
- Macadamia: Illumina run `SRR11191912` and PromethION run `SRR11191910` from
  the earlier accession-1005 study. The update paper establishes same-sample
  provenance at paper level, but differing archival BioSamples do not prove
  the same DNA extraction for the later HiFi data.
- The two Macadamia inputs are official SRA Lite objects. They passed official
  MD5, `vdb-validate`, and converted record/base reconciliation. SRA Lite lacks
  original quality scores; KMC counts sequence only, so qualities are not used
  in this validation.

## Frozen analysis and result

Illumina family abundance used the preregistered family-exclusive diagnostic
k-mers at k=21 and k=31, normalized by the frozen 20,000 empirical controls.
Macadamia ONT supplied direction-only occupancy after the pre-inspection HiFi
mapping-efficiency calibration. Family definitions, eligibility, the 0.6 ratio
threshold and interpretation rules were not retuned.

Among the six pre-orthogonal deficit candidates, three have stable orthogonal
support for residual tandem-repeat under-representation in the newer assembly:
Ey15 `TXF000002` and `TXF000154`, and Macadamia `TXF000496`. Ey15
`TXF001517` and Macadamia `TXF000563`/`TXF000695` remain `unresolved` because
k values or platforms disagree. No family satisfies a stable cross-k
quantification-bias rule, so no TandemX method repair is triggered.

This is a small, deliberately selected set and is not a population accuracy
estimate. Continuous values are read--assembly abundance deficits or estimated
under-representation, not physical missing-base truth. Binary interpretation
is evaluated separately from magnitude.

## Files

- `final_evidence_summary.json`: combined decision and exact input counts.
- `results/*_final_candidates.tsv`: the six candidate-level final decisions.
- `results/*_k{21,31}_family_metrics.tsv`: all catalogue rows, including
  source-ineligible and eligible negative/context rows.
- `results/*_summary.json`: control diagnostics, candidate counts and hashes.
- `ont/`: Macadamia mapping calibration, occupancy and direction-only results.
- `qc/`: official-object integrity and converted FASTQ reconciliation.
- `figures_v1/`: retained failed visual layout; its legend obscured evidence.
- `figures_v2/`: accepted editable SVG/PDF/PNG, complete 72-row panel source,
  legend, hashes and independent PDF/direct-PNG visual review.
- `run_fates.tsv`: retained failed/excluded attempts and the final accepted
  replacement path.
- `provenance.md`: versions, frozen workflow boundary and large-file locations.
- `manifest.tsv`: SHA-256 and byte count for every compact archive member except
  the manifest itself.

The source/metadata selection preceding this run is archived separately at
`paper/evidence/orthogonal_abundance_source_audit_v1`.
