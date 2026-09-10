#!/usr/bin/env python3
"""Remap locked comparison spans against candidate-excluded source reads."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from tandemx.io.sequences import read_sequence_records
from tandemx.recovery.evidence import read_paf
from tandemx.recovery.poc import run_mapping, sha256, write_json, write_table


def remap(proxy: Path, reads: Path, minimap2: Path, outdir: Path) -> None:
    completion = json.loads((proxy / "proxy_completion.json").read_text())
    receipt = json.loads((reads / "receipt.json").read_text())
    if completion.get("status") != "posthoc_proxy_complete" or receipt.get("complete") is not True:
        raise ValueError("Completed proxy and independent-read receipts required")
    spans, fasta = proxy / "comparison_spans.fa", reads / "reads.fa"
    if (sha256(spans) != completion["outputs"]["comparison_spans.fa"]
            or sha256(fasta) != receipt["outputs"]["reads.fa"]):
        raise ValueError("Remapping input receipt/hash mismatch")
    if outdir.exists():
        raise ValueError("Read remapping requires a fresh output directory")
    ids = {r.id for r in read_sequence_records(fasta)}
    if receipt["candidate_read_id_excluded"] in ids or len(ids) != receipt["extracted_independent_read_ids"]:
        raise ValueError("Independent read membership is invalid")
    outdir.mkdir(parents=True)
    write_json(outdir / "input_manifest.json", {"proxy": str(proxy), "reads": str(reads),
               "proxy_completion_sha256": sha256(proxy / "proxy_completion.json"),
               "read_receipt_sha256": sha256(reads / "receipt.json"),
               "source_sha256": sha256(Path(__file__)), "minimap2_sha256": sha256(minimap2),
               "query_coverage_min": .9, "alignment_identity_min": .98,
               "scope": "same_library_source_read_excluded_alignment_support_not_copy_truth"})
    paf = outdir / "spans_to_independent_reads.paf"
    run_mapping([str(minimap2), "-x", "asm5", "-c", "-k", "15", "-w", "5", "-N", "1000",
                 "-p", "0.5", "--secondary=yes", "-t", "4", str(fasta), str(spans)], paf)
    support = defaultdict(set)
    for hit in read_paf(paf):
        if hit.target not in ids:
            raise ValueError("Remapping PAF contains a read outside the enrolled subset")
        if hit.query_coverage >= .9 and hit.identity >= .98:
            support[hit.query].add(hit.target)
    rows = [dict(sequence_id=r.id, span_bp=len(r.sequence), independently_supporting_reads=len(support[r.id]),
                 supplied_independent_reads=len(ids), query_coverage_min=.9, identity_min=.98,
                 warning="same_library_alignment_support_not_independent_copy_number_truth")
            for r in read_sequence_records(spans)]
    write_table(outdir / "read_remapping.tsv", rows)
    write_json(outdir / "remapping_completion.json", {"status": "read_remapping_complete",
               "outputs": {name: sha256(outdir/name) for name in ("read_remapping.tsv", "spans_to_independent_reads.paf")}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("proxy", "reads", "minimap2", "outdir"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    remap(args.proxy, args.reads, args.minimap2, args.outdir)
