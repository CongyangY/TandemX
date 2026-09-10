"""Bounded streaming collection of recovery recruitment evidence."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from tandemx.recovery.evidence import Alignment, read_paf


RECRUITMENT_FIELDS = (
    "read_id", "family_id", "locus_id", "evidence_type", "target",
    "read_length_bp", "read_start", "read_end", "strand", "aligned_bp",
    "identity", "mapq",
)


def collect_recruitment(
    paf: Path,
    anchors: dict[str, dict],
    output_tsv: Path,
    selected_families: set[str],
    maximum_reads: int,
    maximum_rows: int,
) -> tuple[
    dict[str, set[str]],
    dict[str, dict[str, list[Alignment]]],
    dict[str, int],
    int,
]:
    """Filter a recruitment PAF while streaming TSV evidence to a partial file.

    The final TSV is installed only after every accepted alignment is processed.
    On malformed or over-limit input, the partial TSV is intentionally retained
    and the caller receives an exception rather than a silently truncated set.
    """
    if maximum_reads < 1 or maximum_rows < 1:
        raise ValueError("Recruitment limits must be positive")
    if output_tsv.exists():
        raise ValueError(f"Recruitment output already exists: {output_tsv}")
    partial = output_tsv.with_suffix(output_tsv.suffix + ".partial")
    if partial.exists():
        raise ValueError(f"Recruitment partial output exists: {partial}")

    repeat_reads: dict[str, set[str]] = defaultdict(set)
    anchor_hits: dict[str, dict[str, list[Alignment]]] = defaultdict(lambda: defaultdict(list))
    read_lengths: dict[str, int] = {}
    row_count = 0

    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECRUITMENT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for hit in read_paf(paf):
            family: str
            locus: str
            evidence_type: str
            retain_anchor = False
            if hit.target.startswith("repeat_"):
                family = hit.target.removeprefix("repeat_")
                if family not in selected_families:
                    raise ValueError(f"Recruitment PAF names an unselected repeat family: {family}")
                if hit.block_length < 500 or hit.identity < 0.85:
                    continue
                locus, evidence_type = "NA", "repeat_support_only"
            elif hit.target in anchors:
                if hit.target_coverage < 0.9 or hit.identity < 0.98:
                    continue
                anchor = anchors[hit.target]
                family = anchor["family_id"]
                if family not in selected_families:
                    raise ValueError(f"Recruitment anchor has an unselected family: {family}")
                locus, evidence_type, retain_anchor = anchor["locus_id"], "flank_anchored", True
            else:
                raise ValueError(f"Unexpected recruitment target: {hit.target}")

            if row_count >= maximum_rows:
                raise ValueError("Bounded PoC recruitment row limit exceeded; partial TSV retained")
            if hit.query not in read_lengths and len(read_lengths) >= maximum_reads:
                raise ValueError("Bounded PoC recruitment read-ID limit exceeded; partial TSV retained")

            read_lengths[hit.query] = hit.query_length
            if retain_anchor:
                anchor_hits[locus][hit.query].append(hit)
            else:
                repeat_reads[family].add(hit.query)
            writer.writerow({
                "read_id": hit.query,
                "family_id": family,
                "locus_id": locus,
                "evidence_type": evidence_type,
                "target": hit.target,
                "read_length_bp": hit.query_length,
                "read_start": hit.query_start,
                "read_end": hit.query_end,
                "strand": hit.strand,
                "aligned_bp": hit.block_length,
                "identity": hit.identity,
                "mapq": hit.mapq,
            })
            row_count += 1
    partial.replace(output_tsv)
    return dict(repeat_reads), {locus: dict(reads) for locus, reads in anchor_hits.items()}, read_lengths, row_count
