"""Toy-scale tandem repeat discovery MVP."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Iterable, Sequence

from tandemx.simulate.toy import reverse_complement, wrap_sequence


@dataclass(frozen=True)
class FastaRecord:
    read_id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class CandidateRepeat:
    read_id: str
    candidate_id: str
    sequence: str
    read_start: int
    read_end: int
    strand: str
    period_bp: int
    repeat_span_bp: int
    unit_count: float
    score: float
    low_complexity_flag: bool
    confidence: str
    warning: str


@dataclass(frozen=True)
class RepeatFamily:
    family_id: str
    monomer_id: str
    monomer_sequence: str
    monomer_length_bp: int
    support_read_count: int
    support_span_bp: int
    mean_identity: float
    low_complexity_flag: bool
    confidence: str
    warning: str


@dataclass(frozen=True)
class DiscoverConfig:
    reads: Path
    outdir: Path
    min_monomer_len: int
    max_monomer_len: int
    min_support_reads: int
    min_repeat_span: int


def discover_toy_repeats(config: DiscoverConfig) -> tuple[list[CandidateRepeat], list[RepeatFamily]]:
    """Discover toy-scale simple tandem repeat families from FASTA reads."""
    records = list(read_fasta(config.reads))
    if not records:
        raise ValueError("No FASTA records found. discover MVP currently supports FASTA input only.")

    candidates = []
    for record in records:
        candidate = find_best_periodic_candidate(
            record,
            min_period=config.min_monomer_len,
            max_period=config.max_monomer_len,
            min_repeat_span=config.min_repeat_span,
            candidate_index=len(candidates) + 1,
        )
        if candidate is not None:
            candidates.append(candidate)

    families = cluster_candidates(candidates, config.min_support_reads)
    write_candidate_reads(config.outdir / "candidate_reads.tsv", candidates)
    write_monomers(config.outdir / "monomers.fa", families)
    write_families(config.outdir / "families.tsv", families)
    return candidates, families


def read_fasta(path: Path) -> Iterable[FastaRecord]:
    current_header: str | None = None
    sequence_parts: list[str] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                yield make_fasta_record(current_header, sequence_parts)
            current_header = line[1:]
            sequence_parts = []
        elif current_header is None:
            raise ValueError(f"Invalid FASTA: sequence found before header at line {line_number}")
        else:
            sequence = line.upper()
            if any(base not in "ACGTN" for base in sequence):
                raise ValueError(f"Invalid FASTA: unsupported base at line {line_number}")
            sequence_parts.append(sequence)
    if current_header is not None:
        yield make_fasta_record(current_header, sequence_parts)


def make_fasta_record(header: str, sequence_parts: Sequence[str]) -> FastaRecord:
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"Invalid FASTA: empty sequence for record {header}")
    read_id = header.split()[0].split(";")[0]
    return FastaRecord(read_id=read_id, description=header, sequence=sequence)


def find_best_periodic_candidate(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
) -> CandidateRepeat | None:
    sequence = record.sequence.replace("N", "")
    if len(sequence) < max(min_repeat_span, min_period * 2):
        return None
    max_test_period = min(max_period, len(sequence) // 2)
    if max_test_period < min_period:
        return None

    best_period = 0
    best_score = 0.0
    for period in range(min_period, max_test_period + 1):
        score = periodicity_score(sequence, period)
        if score > best_score:
            best_score = score
            best_period = period

    if best_score < 0.75:
        return None

    repeat_span = len(sequence)
    warning = ""
    low_complexity = is_low_complexity(sequence[:best_period])
    if low_complexity:
        warning = "low_complexity_candidate"
    confidence = "high" if best_score >= 0.9 else "medium"
    return CandidateRepeat(
        read_id=record.read_id,
        candidate_id=f"TXC{candidate_index:06d}",
        sequence=sequence[:best_period],
        read_start=0,
        read_end=len(sequence),
        strand=parse_strand(record.description),
        period_bp=best_period,
        repeat_span_bp=repeat_span,
        unit_count=repeat_span / best_period,
        score=best_score,
        low_complexity_flag=low_complexity,
        confidence=confidence,
        warning=warning,
    )


def periodicity_score(sequence: str, period: int) -> float:
    compared = len(sequence) - period
    if compared <= 0:
        return 0.0
    matches = 0
    valid = 0
    for index in range(compared):
        left = sequence[index]
        right = sequence[index + period]
        if left == "N" or right == "N":
            continue
        valid += 1
        if left == right:
            matches += 1
    if valid == 0:
        return 0.0
    return matches / valid


def parse_strand(description: str) -> str:
    for part in description.split(";"):
        if part.startswith("strand="):
            value = part.split("=", 1)[1]
            if value in {"+", "-"}:
                return value
    return "."


def is_low_complexity(sequence: str) -> bool:
    if not sequence:
        return True
    counts = {base: sequence.count(base) for base in "ACGT"}
    return max(counts.values()) / len(sequence) >= 0.8


def cluster_candidates(
    candidates: Sequence[CandidateRepeat],
    min_support_reads: int,
    period_tolerance_bp: int = 5,
) -> list[RepeatFamily]:
    clusters: list[list[CandidateRepeat]] = []
    for candidate in sorted(candidates, key=lambda item: item.period_bp):
        placed = False
        for cluster in clusters:
            center = round(mean(item.period_bp for item in cluster))
            if abs(candidate.period_bp - center) <= period_tolerance_bp:
                cluster.append(candidate)
                placed = True
                break
        if not placed:
            clusters.append([candidate])

    supported = [cluster for cluster in clusters if len({item.read_id for item in cluster}) >= min_support_reads]
    supported.sort(key=lambda cluster: (-len({item.read_id for item in cluster}), round(mean(item.period_bp for item in cluster))))

    families = []
    for index, cluster in enumerate(supported, start=1):
        representative = max(cluster, key=lambda item: (item.score, item.repeat_span_bp))
        monomer_sequence = orient_monomer(representative.sequence)
        mean_identity = mean(item.score for item in cluster)
        low_complexity = any(item.low_complexity_flag for item in cluster)
        warning = "low_complexity_family" if low_complexity else ""
        confidence = "high" if len({item.read_id for item in cluster}) >= max(3, min_support_reads) and mean_identity >= 0.9 else "medium"
        families.append(
            RepeatFamily(
                family_id=f"TXF{index:06d}",
                monomer_id=f"TXM{index:06d}",
                monomer_sequence=monomer_sequence,
                monomer_length_bp=representative.period_bp,
                support_read_count=len({item.read_id for item in cluster}),
                support_span_bp=sum(item.repeat_span_bp for item in cluster),
                mean_identity=mean_identity,
                low_complexity_flag=low_complexity,
                confidence=confidence,
                warning=warning,
            )
        )
    return families


def orient_monomer(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def write_candidate_reads(path: Path, candidates: Sequence[CandidateRepeat]) -> None:
    lines = [
        (
            "read_id\tcandidate_id\tread_start\tread_end\tstrand\tperiod_bp\t"
            "repeat_span_bp\tunit_count\tscore\tlow_complexity_flag\tconfidence\twarning"
        )
    ]
    for candidate in candidates:
        lines.append(
            "\t".join(
                [
                    candidate.read_id,
                    candidate.candidate_id,
                    str(candidate.read_start),
                    str(candidate.read_end),
                    candidate.strand,
                    str(candidate.period_bp),
                    str(candidate.repeat_span_bp),
                    f"{candidate.unit_count:.3f}",
                    f"{candidate.score:.4f}",
                    str(candidate.low_complexity_flag).lower(),
                    candidate.confidence,
                    candidate.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_monomers(path: Path, families: Sequence[RepeatFamily]) -> None:
    lines = []
    for family in families:
        lines.append(
            (
                f">family_id={family.family_id};monomer_id={family.monomer_id};"
                f"length_bp={family.monomer_length_bp};confidence={family.confidence}"
            )
        )
        lines.append(wrap_sequence(family.monomer_sequence))
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_families(path: Path, families: Sequence[RepeatFamily]) -> None:
    lines = [
        (
            "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\t"
            "support_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\t"
            "confidence\twarning"
        )
    ]
    for family in families:
        gc_fraction = (family.monomer_sequence.count("G") + family.monomer_sequence.count("C")) / len(family.monomer_sequence)
        lines.append(
            "\t".join(
                [
                    family.family_id,
                    family.monomer_id,
                    str(family.monomer_length_bp),
                    sequence_md5(family.monomer_sequence),
                    f"{gc_fraction:.4f}",
                    str(family.support_read_count),
                    str(family.support_span_bp),
                    f"{family.mean_identity:.4f}",
                    str(family.low_complexity_flag).lower(),
                    family.confidence,
                    family.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
