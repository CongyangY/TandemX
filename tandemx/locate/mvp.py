"""Toy-scale assembly localization MVP."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from tandemx.discover.mvp import FastaRecord, read_fasta
from tandemx.quantify.mvp import MonomerRecord, canonical_kmer, is_low_complexity_kmer, read_monomer_fasta
from tandemx.compare.mvp import (
    AssemblyReadComparison,
    classify_assembly_read_ratio,
    compare_assembly_to_reads,
    read_copy_number,
    write_comparisons,
)


@dataclass(frozen=True)
class LocateConfig:
    assembly: Path
    monomers: Path
    copy_number: Path | None
    outdir: Path
    window_size: int
    step_size: int
    k: int


@dataclass(frozen=True)
class ArrayHit:
    chrom: str
    start: int
    end: int
    family_id: str
    score: int
    strand: str
    confidence: str
    warning: str


@dataclass(frozen=True)
class DensityWindow:
    chrom: str
    start: int
    end: int
    score: float


def locate_toy_arrays(config: LocateConfig) -> tuple[list[DensityWindow], list[ArrayHit], list[AssemblyReadComparison]]:
    validate_locate_config(config)
    assembly_records = list(read_fasta(config.assembly))
    monomers = list(read_monomer_fasta(config.monomers))
    if not assembly_records:
        raise ValueError("No assembly FASTA records found")
    if not monomers:
        raise ValueError("No monomers found for locate")
    if all(len(record.sequence) < config.k for record in assembly_records):
        raise ValueError("--k is greater than all assembly sequence lengths")
    if all(len(monomer.sequence) < config.k for monomer in monomers):
        raise ValueError("--k is greater than all monomer lengths in the catalogue")

    arrays: list[ArrayHit] = []
    per_chrom_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for record in assembly_records:
        for monomer in monomers:
            intervals = family_hit_intervals(record.sequence, monomer, config.k)
            merged = merge_intervals(intervals, max_gap=config.k * 2)
            for start, end, hit_count in merged:
                min_array_bp = max(len(monomer.sequence) // 2, config.k * 2)
                if end - start < min_array_bp:
                    continue
                score = min(1000, int(1000 * hit_count / max(1, (end - start - config.k + 1))))
                confidence = "high" if score >= 500 else "medium"
                arrays.append(
                    ArrayHit(
                        chrom=record.read_id,
                        start=start,
                        end=end,
                        family_id=monomer.family_id,
                        score=score,
                        strand=".",
                        confidence=confidence,
                        warning="",
                    )
                )
                per_chrom_intervals[record.read_id].append((start, end))

    density = []
    for record in assembly_records:
        density.extend(window_density(record, per_chrom_intervals.get(record.read_id, []), config.window_size, config.step_size))

    comparisons = compare_assembly_to_reads(arrays, config.copy_number)
    write_density(config.outdir / "repeat_density.bedgraph", density)
    write_arrays(config.outdir / "arrays.bed", arrays)
    write_comparisons(config.outdir / "assembly_vs_read_cn.tsv", comparisons)
    return density, arrays, comparisons


def validate_locate_config(config: LocateConfig) -> None:
    if config.window_size <= 0:
        raise ValueError("--window-size must be positive")
    if config.step_size <= 0:
        raise ValueError("--step-size must be positive")
    if config.k <= 0:
        raise ValueError("--k must be positive")


def family_hit_intervals(sequence: str, monomer: MonomerRecord, k: int) -> list[tuple[int, int]]:
    family_kmers = {
        kmer
        for kmer in (canonical_kmer(monomer.sequence[index : index + k]) for index in range(0, len(monomer.sequence) - k + 1))
        if not is_low_complexity_kmer(kmer)
    }
    intervals = []
    sequence = sequence.upper()
    for index in range(0, len(sequence) - k + 1):
        kmer = sequence[index : index + k]
        if "N" in kmer:
            continue
        if canonical_kmer(kmer) in family_kmers:
            intervals.append((index, index + k))
    return intervals


def merge_intervals(intervals: Sequence[tuple[int, int]], max_gap: int = 0) -> list[tuple[int, int, int]]:
    if not intervals:
        return []
    sorted_intervals = sorted(intervals)
    merged = []
    current_start, current_end = sorted_intervals[0]
    hit_count = 1
    for start, end in sorted_intervals[1:]:
        if start <= current_end + max_gap:
            current_end = max(current_end, end)
            hit_count += 1
        else:
            merged.append((current_start, current_end, hit_count))
            current_start, current_end = start, end
            hit_count = 1
    merged.append((current_start, current_end, hit_count))
    return merged


def covered_bp(intervals: Sequence[tuple[int, int]], start: int, end: int) -> int:
    total = 0
    for interval_start, interval_end in intervals:
        overlap_start = max(start, interval_start)
        overlap_end = min(end, interval_end)
        if overlap_end > overlap_start:
            total += overlap_end - overlap_start
    return total


def window_density(record: FastaRecord, intervals: Sequence[tuple[int, int]], window_size: int, step_size: int) -> list[DensityWindow]:
    windows = []
    length = len(record.sequence)
    for start in range(0, length, step_size):
        end = min(length, start + window_size)
        if start >= end:
            continue
        score = covered_bp(intervals, start, end) / (end - start)
        windows.append(DensityWindow(chrom=record.read_id, start=start, end=end, score=score))
        if end == length:
            break
    return windows


def write_density(path: Path, density: Sequence[DensityWindow]) -> None:
    lines = [f"{row.chrom}\t{row.start}\t{row.end}\t{row.score:.4f}" for row in density]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_arrays(path: Path, arrays: Sequence[ArrayHit]) -> None:
    lines = [
        "\t".join(
            [
                array.chrom,
                str(array.start),
                str(array.end),
                array.family_id,
                str(array.score),
                array.strand,
                array.confidence,
                array.warning,
            ]
        )
        for array in arrays
    ]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")

