"""Streaming assembly localization using a shared diagnostic k-mer index."""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from tandemx.compare.mvp import (
    AssemblyReadComparison,
    classify_assembly_read_ratio,
    compare_assembly_to_reads,
    read_copy_number,
    write_comparisons,
)
from tandemx.discover.mvp import FastaRecord
from tandemx.io.sequences import read_fasta_chunks
from tandemx.quantify.mvp import MonomerRecord, read_monomer_fasta
from tandemx.utils.kmers import (
    canonical_kmer_code,
    is_low_complexity_kmer,
    iter_canonical_kmer_codes,
    iter_circular_canonical_kmers,
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
    min_identity: float = 0.8


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
    """Locate all families in one streaming pass over each assembly contig."""
    validate_locate_config(config)
    monomers = list(read_monomer_fasta(config.monomers))
    if not monomers:
        raise ValueError("No monomers found for locate")
    kmer_index, indexed_families, shared_kmers = build_family_kmer_index(monomers, config.k)
    if not kmer_index:
        raise ValueError("No informative monomer k-mers remain after complexity and family-specificity filtering")

    config.outdir.mkdir(parents=True, exist_ok=True)
    arrays: list[ArrayHit] = []
    density: list[DensityWindow] = []
    monomer_lengths = {monomer.family_id: len(monomer.sequence) for monomer in monomers}
    assembly_record_count = 0
    current_id: str | None = None
    current_length = 0
    tail = ""
    active: dict[str, list[int]] = {}
    current_arrays: list[ArrayHit] = []

    def finish_family(family_id: str) -> None:
        state = active.pop(family_id)
        array = array_from_state(
            current_id or "",
            family_id,
            state,
            monomer_lengths[family_id],
            config.k,
            config.min_identity,
            bool(shared_kmers),
        )
        if array is not None:
            current_arrays.append(array)

    def finish_record() -> None:
        nonlocal current_arrays
        if current_id is None:
            return
        for family_id in list(active):
            finish_family(family_id)
        arrays.extend(current_arrays)
        union_intervals = merge_coverage_intervals(
            (item.start, item.end) for item in current_arrays
        )
        density.extend(
            window_density_for_length(
                current_id,
                current_length,
                union_intervals,
                config.window_size,
                config.step_size,
            )
        )
        current_arrays = []

    for chunk in read_fasta_chunks(config.assembly):
        if chunk.id != current_id:
            finish_record()
            assembly_record_count += 1
            current_id = chunk.id
            current_length = 0
            tail = ""
            active.clear()
        combined = tail + chunk.sequence
        combined_start = chunk.start - len(tail)
        first_new_position = max(0, chunk.start - config.k + 1)
        for local_position, code in iter_canonical_kmer_codes(
            combined,
            config.k,
            filter_low_complexity=True,
        ):
            position = combined_start + local_position
            if position < first_new_position:
                continue
            family_id = kmer_index.get(code)
            if family_id is None:
                continue
            state = active.get(family_id)
            if state is None:
                active[family_id] = [position, position + config.k, 1]
            elif position <= state[1] + config.k * 2:
                state[1] = position + config.k
                state[2] += 1
            else:
                finish_family(family_id)
                active[family_id] = [position, position + config.k, 1]
        current_length = chunk.start + len(chunk.sequence)
        tail = combined[-(config.k - 1) :] if config.k > 1 else ""
    finish_record()
    if assembly_record_count == 0:
        raise ValueError("No assembly FASTA records found")

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
    if not 1 <= config.k <= 31:
        raise ValueError("--k must be in 1..31 for bounded-memory rolling k-mer scans")
    if not 0 < config.min_identity <= 1:
        raise ValueError("--min-identity must be in (0, 1]")


def build_family_kmer_index(
    monomers: Sequence[MonomerRecord],
    k: int,
) -> tuple[dict[int, str], set[str], set[int]]:
    """Build one code-to-family index, excluding k-mers shared by multiple families."""
    owners: dict[int, set[str]] = defaultdict(set)
    family_codes: dict[str, set[int]] = defaultdict(set)
    for monomer in monomers:
        for kmer in set(iter_circular_canonical_kmers(monomer.sequence, k)):
            if is_low_complexity_kmer(kmer):
                continue
            code = canonical_kmer_code(kmer)
            owners[code].add(monomer.family_id)
            family_codes[monomer.family_id].add(code)
    shared = {code for code, families in owners.items() if len(families) > 1}
    index = {
        code: next(iter(families))
        for code, families in owners.items()
        if len(families) == 1
    }
    indexed_families = {
        family_id
        for family_id, codes in family_codes.items()
        if any(code not in shared for code in codes)
    }
    return index, indexed_families, shared


def locate_record_arrays(
    record: FastaRecord,
    monomers: Sequence[MonomerRecord],
    kmer_index: dict[int, str],
    indexed_families: set[str],
    shared_kmers: set[int],
    k: int,
    min_identity: float,
) -> list[ArrayHit]:
    hits: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for position, code in iter_canonical_kmer_codes(record.sequence, k, filter_low_complexity=True):
        family_id = kmer_index.get(code)
        if family_id is not None:
            hits[family_id].append((position, position + k))

    monomer_lengths = {monomer.family_id: len(monomer.sequence) for monomer in monomers}
    arrays: list[ArrayHit] = []
    for family_id in sorted(indexed_families):
        intervals = merge_intervals(hits.get(family_id, ()), max_gap=k * 2)
        for start, end, hit_count in intervals:
            array = array_from_state(
                record.read_id,
                family_id,
                [start, end, hit_count],
                monomer_lengths[family_id],
                k,
                min_identity,
                bool(shared_kmers),
            )
            if array is not None:
                arrays.append(array)
    return arrays


def array_from_state(
    chrom: str,
    family_id: str,
    state: Sequence[int],
    monomer_length: int,
    k: int,
    min_identity: float,
    shared_kmers_excluded: bool,
) -> ArrayHit | None:
    start, end, hit_count = state
    min_array_bp = max(monomer_length // 2, k * 2)
    if end - start < min_array_bp:
        return None
    possible_kmers = max(1, end - start - k + 1)
    identity_proxy = min(1.0, hit_count / possible_kmers)
    if identity_proxy < min_identity:
        return None
    warnings = ["exact_kmer_identity_proxy"]
    if shared_kmers_excluded:
        warnings.append("shared_family_kmers_excluded")
    return ArrayHit(
        chrom=chrom,
        start=start,
        end=end,
        family_id=family_id,
        score=min(1000, round(1000 * identity_proxy)),
        strand=".",
        confidence="medium",
        warning=";".join(warnings),
    )


def family_hit_intervals(sequence: str, monomer: MonomerRecord, k: int) -> list[tuple[int, int]]:
    """Compatibility helper using circular monomer k-mers and a rolling sequence scan."""
    codes = {
        canonical_kmer_code(kmer)
        for kmer in iter_circular_canonical_kmers(monomer.sequence, k)
        if not is_low_complexity_kmer(kmer)
    }
    return [
        (position, position + k)
        for position, code in iter_canonical_kmer_codes(sequence, k, filter_low_complexity=True)
        if code in codes
    ]


def merge_intervals(intervals: Sequence[tuple[int, int]], max_gap: int = 0) -> list[tuple[int, int, int]]:
    if not intervals:
        return []
    sorted_intervals = sorted(intervals)
    merged: list[tuple[int, int, int]] = []
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


def merge_coverage_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return []
    merged: list[tuple[int, int]] = []
    start, end = sorted_intervals[0]
    for next_start, next_end in sorted_intervals[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def covered_bp(intervals: Sequence[tuple[int, int]], start: int, end: int) -> int:
    """Return union coverage, so overlapping family arrays are never double counted."""
    total = 0
    for interval_start, interval_end in merge_coverage_intervals(intervals):
        overlap_start = max(start, interval_start)
        overlap_end = min(end, interval_end)
        if overlap_end > overlap_start:
            total += overlap_end - overlap_start
    return total


def window_density(record: FastaRecord, intervals: Sequence[tuple[int, int]], window_size: int, step_size: int) -> list[DensityWindow]:
    """Compute window union coverage with binary search over merged intervals."""
    return window_density_for_length(
        record.read_id,
        len(record.sequence),
        intervals,
        window_size,
        step_size,
    )


def window_density_for_length(
    chrom: str,
    length: int,
    intervals: Sequence[tuple[int, int]],
    window_size: int,
    step_size: int,
) -> list[DensityWindow]:
    merged = merge_coverage_intervals(intervals)
    starts = [start for start, _ in merged]
    prefix = [0]
    for start, end in merged:
        prefix.append(prefix[-1] + end - start)

    def coverage_before(position: int) -> int:
        count = bisect_left(starts, position)
        total = prefix[count]
        if count and merged[count - 1][1] > position:
            total -= merged[count - 1][1] - position
        return total

    def fast_covered(window_start: int, window_end: int) -> int:
        return max(0, coverage_before(window_end) - coverage_before(window_start))

    windows: list[DensityWindow] = []
    for start in range(0, length, step_size):
        end = min(length, start + window_size)
        if start >= end:
            continue
        score = fast_covered(start, end) / (end - start)
        windows.append(DensityWindow(chrom=chrom, start=start, end=end, score=score))
        if end == length:
            break
    return windows


def write_density(path: Path, density: Sequence[DensityWindow]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        for row in density:
            handle.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.score:.4f}\n")


def write_arrays(path: Path, arrays: Sequence[ArrayHit]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        for array in arrays:
            handle.write(
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
                + "\n"
            )
