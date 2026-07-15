"""Bounded-memory FISH probe ranking with one indexed assembly scan."""

from __future__ import annotations

import math
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from tandemx.io.sequences import read_fasta_chunks
from tandemx.locate.mvp import merge_coverage_intervals, merge_intervals
from tandemx.quantify.mvp import canonical_kmer, is_low_complexity_kmer, read_monomer_fasta
from tandemx.simulate.toy import wrap_sequence
from tandemx.utils.kmers import canonical_kmer_code, iter_canonical_kmer_codes


@dataclass(frozen=True)
class ProbeConfig:
    monomers: Path
    assembly: Path
    copy_number: Path
    arrays: Path
    outdir: Path
    min_len: int
    max_len: int
    sodium_millimolar: float = 50.0
    formamide_percent: float = 0.0


@dataclass(frozen=True)
class ArrayRegion:
    chrom: str
    start: int
    end: int
    family_id: str


@dataclass(frozen=True)
class ProbeCandidate:
    probe_id: str
    family_id: str
    sequence: str
    sequence_length: int
    gc_content: float
    tm: float
    estimated_copy_number: float
    arrayiness_score: float
    specificity_score: float
    off_target_hits: int
    predicted_regions: str
    probe_score: float
    confidence: str
    warning: str


@dataclass(frozen=True)
class FishSignal:
    probe_id: str
    chrom: str
    start: int
    end: int
    predicted_signal: float
    confidence: str
    warning: str


@dataclass(frozen=True)
class _ProbeSeed:
    probe_id: str
    family_id: str
    sequence: str


def rank_toy_probes(config: ProbeConfig) -> tuple[list[ProbeCandidate], list[FishSignal]]:
    validate_probe_config(config)
    config.outdir.mkdir(parents=True, exist_ok=True)
    monomers = list(read_monomer_fasta(config.monomers))
    if not monomers:
        raise ValueError("No monomers found for probe ranking")
    copy_numbers = read_copy_numbers(config.copy_number)
    arrays = read_arrays(config.arrays)
    array_index = build_array_index(arrays)

    seeds: list[_ProbeSeed] = []
    for monomer in monomers:
        for sequence in candidate_probe_sequences(monomer.sequence, config.min_len, config.max_len):
            if low_complexity_ratio(sequence) >= 0.8:
                continue
            seeds.append(
                _ProbeSeed(
                    probe_id=f"TXP{len(seeds) + 1:06d}",
                    family_id=monomer.family_id,
                    sequence=sequence,
                )
            )
    if not seeds:
        write_probes_fasta(config.outdir / "probes.fa", ())
        write_probe_ranks(config.outdir / "probes.rank.tsv", ())
        write_fish_signals(config.outdir / "in_silico_fish.tsv", ())
        return [], []

    k = min(21, max(8, config.min_len // 4))
    regions_by_probe = indexed_probe_regions(seeds, config.assembly, k)
    candidates: list[ProbeCandidate] = []
    signals: list[FishSignal] = []
    for index, seed in enumerate(seeds):
        hit_regions = regions_by_probe[index]
        target_regions = [
            region
            for region in hit_regions
            if indexed_region_overlap(region, seed.family_id, array_index)
        ]
        off_target_hits = len(hit_regions) - len(target_regions)
        arrayiness = len(target_regions) / len(hit_regions) if hit_regions else 0.0
        specificity = 1.0 / (1.0 + off_target_hits)
        copy_number = copy_numbers.get(seed.family_id, 0.0)
        gc = gc_content(seed.sequence)
        tm = long_oligo_tm(
            seed.sequence,
            sodium_millimolar=config.sodium_millimolar,
            formamide_percent=config.formamide_percent,
        )
        gc_balance = max(0.0, 1.0 - abs(gc - 0.5))
        copy_component = min(1.0, copy_number / 10.0)
        probe_score = copy_component * specificity * max(arrayiness, 0.1) * gc_balance
        predicted_regions = ";".join(f"{chrom}:{start}-{end}" for chrom, start, end in target_regions)
        warnings = ["heuristic_probe_score_not_experimentally_calibrated"]
        if not target_regions:
            warnings.append("no_predicted_array_region")
        if off_target_hits:
            warnings.append("predicted_off_target_regions")
        confidence = "medium" if target_regions and off_target_hits == 0 else "low"
        candidate = ProbeCandidate(
            probe_id=seed.probe_id,
            family_id=seed.family_id,
            sequence=seed.sequence,
            sequence_length=len(seed.sequence),
            gc_content=gc,
            tm=tm,
            estimated_copy_number=copy_number,
            arrayiness_score=arrayiness,
            specificity_score=specificity,
            off_target_hits=off_target_hits,
            predicted_regions=predicted_regions,
            probe_score=probe_score,
            confidence=confidence,
            warning=";".join(warnings),
        )
        candidates.append(candidate)
        for chrom, start, end in target_regions:
            signals.append(
                FishSignal(
                    probe_id=seed.probe_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    predicted_signal=probe_score,
                    confidence=confidence,
                    warning="heuristic_signal_prediction",
                )
            )

    candidates.sort(key=lambda item: (-item.probe_score, item.family_id, item.probe_id))
    write_probes_fasta(config.outdir / "probes.fa", candidates)
    write_probe_ranks(config.outdir / "probes.rank.tsv", candidates)
    write_fish_signals(config.outdir / "in_silico_fish.tsv", signals)
    return candidates, signals


def validate_probe_config(config: ProbeConfig) -> None:
    if config.min_len <= 0:
        raise ValueError("--min-len must be positive")
    if config.max_len < config.min_len:
        raise ValueError("--max-len must be greater than or equal to --min-len")
    if config.sodium_millimolar <= 0:
        raise ValueError("--sodium-millimolar must be positive")
    if not 0 <= config.formamide_percent <= 100:
        raise ValueError("--formamide-percent must be in [0, 100]")


def candidate_probe_sequences(sequence: str, min_len: int, max_len: int) -> list[str]:
    """Generate deterministic probes, using tandem context for short monomers."""
    normalized = sequence.upper()
    if not normalized:
        return []
    probe_len = min(max_len, max(min_len, min(len(normalized), max_len)))
    if len(normalized) < probe_len:
        repeats = (probe_len + len(normalized) - 2) // len(normalized) + 1
        source = normalized * repeats
        raw = [source[phase : phase + probe_len] for phase in range(len(normalized))]
    else:
        source = normalized
        step = max(1, probe_len // 2)
        raw = [source[start : start + probe_len] for start in range(0, len(source) - probe_len + 1, step)]
        if source[-probe_len:] not in raw:
            raw.append(source[-probe_len:])
    return list(dict.fromkeys(raw))


def gc_content(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def simple_tm(sequence: str) -> float:
    """Legacy Wallace-rule helper retained for short-primer compatibility."""
    gc = sequence.count("G") + sequence.count("C")
    at = sequence.count("A") + sequence.count("T")
    return 4.0 * gc + 2.0 * at


def long_oligo_tm(
    sequence: str,
    *,
    sodium_millimolar: float = 50.0,
    formamide_percent: float = 0.0,
) -> float:
    """Salt/formamide-adjusted approximation suitable for ranking long DNA probes."""
    if len(sequence) < 14:
        return simple_tm(sequence)
    gc_percent = 100.0 * gc_content(sequence)
    sodium_molar = sodium_millimolar / 1000.0
    return 81.5 + 16.6 * math.log10(sodium_molar) + 0.41 * gc_percent - 600.0 / len(sequence) - 0.62 * formamide_percent


def low_complexity_ratio(sequence: str) -> float:
    if not sequence:
        return 1.0
    counts = {base: sequence.count(base) for base in "ACGT"}
    return max(counts.values()) / len(sequence)


def read_copy_numbers(path: Path) -> dict[str, float]:
    with path.open("rt", encoding="utf-8") as handle:
        header_line = handle.readline().rstrip("\n\r")
        if not header_line:
            raise ValueError(f"Copy-number table is empty: {path}")
        header = header_line.split("\t")
        required = {"family_id", "estimated_copy_number"}
        missing = sorted(required - set(header))
        if missing:
            raise ValueError(f"{path} is missing required field(s): {', '.join(missing)}")
        family_index = header.index("family_id")
        copy_index = header.index("estimated_copy_number")
        values: dict[str, float] = {}
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            try:
                values[parts[family_index]] = float(parts[copy_index])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid copy-number value in {path} at line {line_number}") from exc
    if not values:
        raise ValueError(f"Copy-number table has no records: {path}")
    return values


def read_arrays(path: Path) -> list[ArrayRegion]:
    arrays: list[ArrayRegion] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 4:
                raise ValueError(f"Invalid arrays BED record in {path} at line {line_number}")
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as exc:
                raise ValueError(f"Invalid arrays BED coordinates in {path} at line {line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid 0-based half-open interval in {path} at line {line_number}")
            arrays.append(ArrayRegion(chrom=parts[0], start=start, end=end, family_id=parts[3]))
    if not arrays:
        raise ValueError(f"Array BED file has no records: {path}")
    return arrays


def build_array_index(arrays: Sequence[ArrayRegion]) -> dict[tuple[str, str], tuple[list[tuple[int, int]], list[int]]]:
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for array in arrays:
        grouped[(array.family_id, array.chrom)].append((array.start, array.end))
    return {
        key: (merged, [start for start, _ in merged])
        for key, values in grouped.items()
        for merged in [merge_coverage_intervals(values)]
    }


def indexed_region_overlap(
    region: tuple[str, int, int],
    family_id: str,
    array_index: dict[tuple[str, str], tuple[list[tuple[int, int]], list[int]]],
) -> bool:
    chrom, start, end = region
    entry = array_index.get((family_id, chrom))
    if entry is None:
        return False
    intervals, starts = entry
    position = bisect_left(starts, end)
    if position == 0:
        return False
    array_start, array_end = intervals[position - 1]
    return array_end > start and array_start < end


def indexed_probe_regions(
    seeds: Sequence[_ProbeSeed],
    assembly_path: Path,
    k: int,
) -> list[list[tuple[str, int, int]]]:
    code_to_probes: dict[int, list[int]] = defaultdict(list)
    kmer_counts: list[int] = []
    for probe_index, seed in enumerate(seeds):
        codes = {
            canonical_kmer_code(canonical_kmer(seed.sequence[index : index + k]))
            for index in range(0, len(seed.sequence) - k + 1)
            if not is_low_complexity_kmer(seed.sequence[index : index + k])
        }
        kmer_counts.append(len(codes))
        for code in codes:
            code_to_probes[code].append(probe_index)

    regions: list[list[tuple[str, int, int]]] = [[] for _ in seeds]
    assembly_count = 0
    current_id: str | None = None
    tail = ""
    active: dict[int, list[int]] = {}

    def finish(probe_index: int) -> None:
        state = active.pop(probe_index)
        threshold = max(2, kmer_counts[probe_index] // 10)
        if state[2] >= threshold:
            regions[probe_index].append((current_id or "", state[0], state[1]))

    def finish_record() -> None:
        for probe_index in list(active):
            finish(probe_index)

    for chunk in read_fasta_chunks(assembly_path):
        if chunk.id != current_id:
            finish_record()
            current_id = chunk.id
            tail = ""
            assembly_count += 1
        combined = tail + chunk.sequence
        combined_start = chunk.start - len(tail)
        first_new_position = max(0, chunk.start - k + 1)
        for local_position, code in iter_canonical_kmer_codes(
            combined,
            k,
            filter_low_complexity=True,
        ):
            position = combined_start + local_position
            if position < first_new_position:
                continue
            for probe_index in code_to_probes.get(code, ()):
                state = active.get(probe_index)
                if state is None:
                    active[probe_index] = [position, position + k, 1]
                elif position <= state[1] + k * 2:
                    state[1] = position + k
                    state[2] += 1
                else:
                    finish(probe_index)
                    active[probe_index] = [position, position + k, 1]
        tail = combined[-(k - 1) :] if k > 1 else ""
    finish_record()
    if assembly_count == 0:
        raise ValueError("No assembly records found for probe ranking")
    return regions


def find_probe_regions(sequence: str, assembly_records: Sequence) -> list[tuple[str, int, int]]:
    """Compatibility helper for already-loaded toy assembly records."""
    k = min(21, max(8, len(sequence) // 4))
    kmers = {
        canonical_kmer(sequence[index : index + k])
        for index in range(0, len(sequence) - k + 1)
        if not is_low_complexity_kmer(sequence[index : index + k])
    }
    if not kmers:
        return []
    regions: list[tuple[str, int, int]] = []
    for record in assembly_records:
        intervals = []
        for index in range(0, len(record.sequence) - k + 1):
            kmer = record.sequence[index : index + k]
            if "N" not in kmer and canonical_kmer(kmer) in kmers:
                intervals.append((index, index + k))
        for start, end, hit_count in merge_intervals(intervals, max_gap=k * 2):
            if hit_count >= max(2, len(kmers) // 10):
                regions.append((record.read_id, start, end))
    return regions


def overlaps_family_array(region: tuple[str, int, int], family_id: str, arrays: Sequence[ArrayRegion]) -> bool:
    return indexed_region_overlap(region, family_id, build_array_index(arrays))


def write_probes_fasta(path: Path, candidates: Sequence[ProbeCandidate]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        for candidate in candidates:
            handle.write(
                f">probe_id={candidate.probe_id};family_id={candidate.family_id};"
                f"length_bp={candidate.sequence_length};probe_score={candidate.probe_score:.4f};"
                f"confidence={candidate.confidence}\n"
            )
            handle.write(wrap_sequence(candidate.sequence) + "\n")


def write_probe_ranks(path: Path, candidates: Sequence[ProbeCandidate]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        handle.write(
            "probe_id\tfamily_id\tsequence_length\tgc_content\ttm\testimated_copy_number\t"
            "arrayiness_score\tspecificity_score\toff_target_hits\tpredicted_regions\tprobe_score\tconfidence\twarning\n"
        )
        for candidate in candidates:
            handle.write(
                "\t".join(
                    [
                        candidate.probe_id,
                        candidate.family_id,
                        str(candidate.sequence_length),
                        f"{candidate.gc_content:.4f}",
                        f"{candidate.tm:.2f}",
                        f"{candidate.estimated_copy_number:.4f}",
                        f"{candidate.arrayiness_score:.4f}",
                        f"{candidate.specificity_score:.4f}",
                        str(candidate.off_target_hits),
                        candidate.predicted_regions,
                        f"{candidate.probe_score:.4f}",
                        candidate.confidence,
                        candidate.warning,
                    ]
                )
                + "\n"
            )


def write_fish_signals(path: Path, signals: Sequence[FishSignal]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        handle.write("probe_id\tchrom\tstart\tend\tpredicted_signal\tconfidence\twarning\n")
        for signal in signals:
            handle.write(
                "\t".join(
                    [
                        signal.probe_id,
                        signal.chrom,
                        str(signal.start),
                        str(signal.end),
                        f"{signal.predicted_signal:.4f}",
                        signal.confidence,
                        signal.warning,
                    ]
                )
                + "\n"
            )
