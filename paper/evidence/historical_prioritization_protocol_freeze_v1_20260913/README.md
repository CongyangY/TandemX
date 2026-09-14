# Historical prioritization comparator protocol freeze

This directory records the pre-execution validation of the historical-specific
SRF/ordinary-mapping protocol frozen on 2026-09-13.

`binding_audit.json` is the runner's read-only full-file audit. It recomputed
SHA-256 for the complete Ey15-2 FASTQ, both complete Macadamia FASTQs in the
frozen order, both complete TandemX catalogues, both primary endpoint tables,
both source configurations, the historical and upstream simulation protocols,
the benchmark implementation sources, and all six native tools/scripts.

`validation_receipt.json` records the structural plan and focused tests. No SRF,
KMC, minimap2 native comparator cell was launched. No pre-existing PAF or native
abundance output was imported or scored. Formal results therefore remain
unavailable pending main-task review and the explicit execution acknowledgement
defined by the protocol.
