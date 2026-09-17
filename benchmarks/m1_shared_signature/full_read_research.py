"""Research-only full-read competitive occupancy with conservative abstention.

Reads are streamed one at a time. Overlapping context and disjoint core tiles
are both classified; core bp are assigned only when both agree uniquely.
This is not a calibrated copy-number backend.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from tandemx.io.sequences import SequenceRecord, read_sequence_records_many

from .occupancy_research import CompetitiveConfig, classify_window, load_catalogue


@dataclass(frozen=True)
class FullReadConfig:
    core_bp: int = 40
    flank_bp: int = 40
    max_edit_fraction: float = 0.225
    min_margin_edits: int = 2
    max_read_bp: int = 1_000_000
    max_families: int = 8


def classify_read(sequence: str, catalogue: dict[str, str], config: FullReadConfig
                  ) -> list[dict]:
    """Partition a read exactly; each row owns disjoint 0-based half-open bp."""
    sequence = sequence.upper()
    if not sequence or len(sequence) > config.max_read_bp:
        raise ValueError("Empty read or research read-length cap exceeded")
    if config.core_bp < 1 or config.flank_bp < 0 or not 0 <= config.max_edit_fraction < 0.5:
        raise ValueError("Invalid full-read tile settings")
    if config.min_margin_edits < 1 or not catalogue or len(catalogue) > config.max_families:
        raise ValueError("Invalid full-read catalogue or margin")
    probe = CompetitiveConfig(
        reads=Path("unused"), monomers=Path("unused"), genome_size=1,
        outdir=Path("unused"), min_tail_bp=1,
        max_edit_fraction=config.max_edit_fraction,
        min_margin_edits=config.min_margin_edits,
    )
    segments = []
    for start in range(0, len(sequence), config.core_bp):
        end = min(start + config.core_bp, len(sequence))
        left = max(0, start - config.flank_bp)
        right = min(len(sequence), end + config.flank_bp)
        core_status, core_id = classify_window(sequence[start:end], catalogue, probe)
        context_status, context_id = classify_window(sequence[left:right], catalogue, probe)
        if core_status == context_status == "assigned" and core_id == context_id:
            status, identity = "assigned", core_id
        elif core_status == "unknown" or context_status == "unknown":
            status, identity = "unknown", None
        else:
            status = "ambiguous"
            identities = set((core_id or "").split("+")) | set((context_id or "").split("+"))
            identity = "+".join(sorted(identities - {""})) or None
        segments.append(dict(start=start, end=end, status=status, identity=identity))
    assert sum(row["end"] - row["start"] for row in segments) == len(sequence)
    return segments


def summarize_records(records: Iterable[SequenceRecord], catalogue: dict[str, str],
                      config: FullReadConfig) -> dict:
    """Only counters and the current read are retained across the stream."""
    counts: Counter[str] = Counter()
    family_bp: Counter[str] = Counter({family: 0 for family in catalogue})
    groups: Counter[str] = Counter()
    read_count = 0
    total_bp = 0
    for record in records:
        read_count += 1
        total_bp += len(record.sequence)
        for segment in classify_read(record.sequence, catalogue, config):
            bp = segment["end"] - segment["start"]
            counts[segment["status"]] += bp
            if segment["status"] == "assigned":
                family_bp[segment["identity"]] += bp
            elif segment["status"] == "ambiguous":
                groups[segment["identity"] or "unresolved"] += bp
    if read_count == 0:
        raise ValueError("No reads")
    if sum(counts.values()) != total_bp:
        raise AssertionError("Full-read mass was not conserved")
    return dict(read_count=read_count, total_read_bp=total_bp,
                assigned_read_bp=counts["assigned"], ambiguous_read_bp=counts["ambiguous"],
                unknown_read_bp=counts["unknown"], assigned_by_family=dict(family_bp),
                ambiguity_groups=dict(sorted(groups.items())),
                boundary="Research occupancy only; no physical-copy calibration")


def stream_files(reads: Path, monomers: Path, config: FullReadConfig) -> dict:
    catalogue = load_catalogue(monomers, CompetitiveConfig(
        reads=reads, monomers=monomers, genome_size=1, outdir=Path("unused"),
        max_catalogue_families=config.max_families))
    return summarize_records(read_sequence_records_many([reads]), catalogue, config)
