"""Coordinate and evidence rules for recovery without a newer assembly."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable

from tandemx.io.sequences import read_fasta_chunks


@dataclass(frozen=True)
class Alignment:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int

    @property
    def identity(self) -> float:
        return self.matches / self.block_length

    @property
    def query_coverage(self) -> float:
        return (self.query_end - self.query_start) / self.query_length

    @property
    def target_coverage(self) -> float:
        return (self.target_end - self.target_start) / self.target_length


def read_paf(path: Path) -> Iterable[Alignment]:
    """Reject malformed mappings rather than silently dropping evidence."""
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF at {path}:{number}")
            try:
                hit = Alignment(fields[0], *(int(x) for x in fields[1:4]), fields[4],
                                fields[5], *(int(x) for x in fields[6:12]))
            except ValueError as exc:
                raise ValueError(f"Invalid PAF number at {path}:{number}") from exc
            if (not hit.query or not hit.target or hit.strand not in {"+", "-"}
                    or not 0 <= hit.query_start < hit.query_end <= hit.query_length
                    or not 0 <= hit.target_start < hit.target_end <= hit.target_length
                    or not 0 <= hit.matches <= hit.block_length or hit.block_length <= 0
                    or not 0 <= hit.mapq <= 255):
                raise ValueError(f"Invalid PAF coordinates at {path}:{number}")
            yield hit


def historical_loci(path: Path, families: set[str]) -> list[dict]:
    """Keep all old-assembly loci; coalesce overlapping same-family intervals."""
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 4:
                raise ValueError(f"Malformed BED at line {number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at line {number}")
            if fields[3] in families:
                grouped.setdefault((fields[3], fields[0]), []).append((start, end))
    loci = []
    for (family, chrom), intervals in sorted(grouped.items()):
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        for start, end in merged:
            loci.append(dict(locus_id=f"L{len(loci)+1:04d}", family_id=family,
                             chromosome=chrom, start=start, end=end))
    return loci


def extract_intervals(assembly: Path, intervals: list[dict]) -> dict[str, str]:
    """Extract only requested windows, streaming even chromosome-scale records."""
    grouped: dict[str, list[dict]] = {}
    parts: dict[str, list[str]] = {}
    for row in intervals:
        if row["start"] < 0 or row["end"] <= row["start"]:
            continue
        grouped.setdefault(row["chromosome"], []).append(row)
        parts[row["anchor_id"]] = []
    for chunk in read_fasta_chunks(assembly):
        for row in grouped.get(chunk.id, []):
            start, end = max(row["start"], chunk.start), min(row["end"], chunk.start + len(chunk.sequence))
            if start < end:
                parts[row["anchor_id"]].append(chunk.sequence[start-chunk.start:end-chunk.start])
    return {key: "".join(value) for key, value in parts.items()}


def make_flanks(loci: list[dict], length: int = 2000,
                offsets: tuple[int, ...] = (0, 2000, 5000, 10000)) -> list[dict]:
    if length < 1 or any(offset < 0 for offset in offsets):
        raise ValueError("Flank length/offset invalid")
    rows = []
    for locus in loci:
        for side in ("left", "right"):
            for offset in offsets:
                start = locus["start"] - offset - length if side == "left" else locus["end"] + offset
                rows.append(dict(anchor_id=f'{locus["locus_id"]}_{side}_{offset}',
                                 locus_id=locus["locus_id"], family_id=locus["family_id"],
                                 chromosome=locus["chromosome"], start=start, end=start+length,
                                 side=side, offset=offset))
    return rows


def anchor_uniqueness(anchor: dict, sequence: str, hits: list[Alignment]) -> dict:
    """Mapping-based uniqueness relative to the old assembly, not genomic truth."""
    result = dict(anchor, sequence_length=len(sequence), eligible=False,
                  qualifying_hits=0, uniqueness="unresolved", reason="")
    if len(sequence) != anchor["end"]-anchor["start"] or anchor["start"] < 0:
        return dict(result, reason="incomplete_or_out_of_contig_flank")
    if set(sequence.upper()) - set("ACGT"):
        return dict(result, reason="ambiguous_bases_in_flank")
    qualifying = [h for h in hits if h.query_coverage >= .90 and h.identity >= .95]
    result["qualifying_hits"] = len(qualifying)
    expected = [h for h in qualifying if h.target == anchor["chromosome"]
                and h.strand == "+" and h.identity >= .98
                and abs(h.target_start - anchor["start"]) <= 200
                and abs(h.target_end - anchor["end"]) <= 200]
    if len(qualifying) == 1 and len(expected) == 1:
        return dict(result, eligible=True, uniqueness="unique_in_old_assembly",
                    reason="assembly_relative_mapping_uniqueness_not_copy_truth")
    return dict(result, uniqueness="nonunique" if len(qualifying) > 1 else "unresolved",
                reason="competing_old_assembly_hits" if len(qualifying) > 1 else "expected_self_hit_missing")


def spanning_interval(left: Alignment, right: Alignment) -> tuple[int, int, str] | None:
    """Return forward-old-orientation read interval bounded by anchor interiors."""
    if left.query != right.query or left.strand != right.strand:
        return None
    if left.target_coverage < .90 or right.target_coverage < .90:
        return None
    if left.identity < .98 or right.identity < .98:
        return None
    if left.strand == "+" and left.query_end < right.query_start:
        return left.query_end, right.query_start, "+"
    if left.strand == "-" and right.query_end < left.query_start:
        return right.query_end, left.query_start, "-"
    return None


def select_spanning_candidate(spans: list[dict]) -> tuple[dict | None, str]:
    """Select an observed median span only; never synthesize copy-count consensus."""
    by_read: dict[str, dict] = {}
    for span in spans:
        if span["end"] <= span["start"]:
            raise ValueError("Invalid spanning read coordinates")
        prior = by_read.get(span["read_id"])
        if prior and prior != span:
            return None, "unresolved_conflicting_paths"
        by_read[span["read_id"]] = span
    values = list(by_read.values())
    if len(values) < 3:
        return None, "insufficient_read_support"
    lengths = [v["end"]-v["start"] for v in values]
    center = median(lengths)
    if max(lengths)-min(lengths) > max(100, .01*center):
        return None, "unresolved_conflicting_paths"
    chosen = min(values, key=lambda row: (abs(row["end"]-row["start"]-center), row["read_id"]))
    return chosen, "partially_resolved"
