"""Toy-scale FISH probe ranking MVP."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from tandemx.discover.mvp import read_fasta
from tandemx.locate.mvp import merge_intervals
from tandemx.quantify.mvp import MonomerRecord, canonical_kmer, is_low_complexity_kmer, read_monomer_fasta
from tandemx.simulate.toy import reverse_complement, wrap_sequence


@dataclass(frozen=True)
class ProbeConfig:
    monomers: Path
    assembly: Path
    copy_number: Path
    arrays: Path
    outdir: Path
    min_len: int
    max_len: int


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
    warning: str


@dataclass(frozen=True)
class FishSignal:
    probe_id: str
    chrom: str
    start: int
    end: int
    predicted_signal: float
    confidence: str


def rank_toy_probes(config: ProbeConfig) -> tuple[list[ProbeCandidate], list[FishSignal]]:
    validate_probe_config(config)
    monomers = list(read_monomer_fasta(config.monomers))
    assembly = list(read_fasta(config.assembly))
    if not monomers:
        raise ValueError("No monomers found for probe ranking")
    if not assembly:
        raise ValueError("No assembly records found for probe ranking")
    copy_numbers = read_copy_numbers(config.copy_number)
    arrays = read_arrays(config.arrays)

    candidates: list[ProbeCandidate] = []
    signals: list[FishSignal] = []
    for monomer in monomers:
        for sequence in candidate_probe_sequences(monomer.sequence, config.min_len, config.max_len):
            if low_complexity_ratio(sequence) >= 0.8:
                continue
            probe_id = f"TXP{len(candidates) + 1:06d}"
            hit_regions = find_probe_regions(sequence, assembly)
            target_regions = [region for region in hit_regions if overlaps_family_array(region, monomer.family_id, arrays)]
            off_target_hits = max(0, len(hit_regions) - len(target_regions))
            arrayiness = len(target_regions) / len(hit_regions) if hit_regions else 0.0
            specificity = 1.0 / (1.0 + off_target_hits)
            copy_number = copy_numbers.get(monomer.family_id, 0.0)
            gc = gc_content(sequence)
            tm = simple_tm(sequence)
            gc_balance = max(0.0, 1.0 - abs(gc - 0.5))
            copy_component = min(1.0, copy_number / 10.0)
            probe_score = copy_component * specificity * max(arrayiness, 0.1) * gc_balance
            predicted_regions = ";".join(f"{chrom}:{start}-{end}" for chrom, start, end in target_regions)
            warning = "" if target_regions else "no_predicted_array_region"
            candidate = ProbeCandidate(
                probe_id=probe_id,
                family_id=monomer.family_id,
                sequence=sequence,
                sequence_length=len(sequence),
                gc_content=gc,
                tm=tm,
                estimated_copy_number=copy_number,
                arrayiness_score=arrayiness,
                specificity_score=specificity,
                off_target_hits=off_target_hits,
                predicted_regions=predicted_regions,
                probe_score=probe_score,
                warning=warning,
            )
            candidates.append(candidate)
            for chrom, start, end in target_regions:
                signals.append(
                    FishSignal(
                        probe_id=probe_id,
                        chrom=chrom,
                        start=start,
                        end=end,
                        predicted_signal=probe_score,
                        confidence="medium" if warning else "high",
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


def candidate_probe_sequences(sequence: str, min_len: int, max_len: int) -> list[str]:
    if len(sequence) < min_len:
        return []
    probe_len = min(max_len, len(sequence))
    step = max(1, probe_len // 2)
    candidates = []
    for start in range(0, len(sequence) - probe_len + 1, step):
        candidates.append(sequence[start : start + probe_len])
    if sequence[-probe_len:] not in candidates:
        candidates.append(sequence[-probe_len:])
    return candidates


def gc_content(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def simple_tm(sequence: str) -> float:
    gc = sequence.count("G") + sequence.count("C")
    at = sequence.count("A") + sequence.count("T")
    return 4.0 * gc + 2.0 * at


def low_complexity_ratio(sequence: str) -> float:
    if not sequence:
        return 1.0
    counts = {base: sequence.count(base) for base in "ACGT"}
    return max(counts.values()) / len(sequence)


def read_copy_numbers(path: Path) -> dict[str, float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return {}
    header = lines[0].split("\t")
    family_index = header.index("family_id")
    copy_index = header.index("estimated_copy_number")
    values = {}
    for line in lines[1:]:
        if not line:
            continue
        parts = line.split("\t")
        values[parts[family_index]] = float(parts[copy_index])
    return values


def read_arrays(path: Path) -> list[ArrayRegion]:
    arrays = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        parts = line.split("\t")
        arrays.append(ArrayRegion(chrom=parts[0], start=int(parts[1]), end=int(parts[2]), family_id=parts[3]))
    return arrays


def find_probe_regions(sequence: str, assembly_records: Sequence) -> list[tuple[str, int, int]]:
    k = min(21, max(8, len(sequence) // 4))
    kmers = {
        canonical_kmer(sequence[index : index + k])
        for index in range(0, len(sequence) - k + 1)
        if not is_low_complexity_kmer(sequence[index : index + k])
    }
    if not kmers:
        return []
    regions = []
    for record in assembly_records:
        intervals = []
        assembly = record.sequence.upper()
        for index in range(0, len(assembly) - k + 1):
            kmer = assembly[index : index + k]
            if "N" in kmer:
                continue
            if canonical_kmer(kmer) in kmers:
                intervals.append((index, index + k))
        for start, end, hit_count in merge_intervals(intervals, max_gap=k * 2):
            if hit_count >= max(2, len(kmers) // 10):
                regions.append((record.read_id, start, end))
    return regions


def overlaps_family_array(region: tuple[str, int, int], family_id: str, arrays: Sequence[ArrayRegion]) -> bool:
    chrom, start, end = region
    for array in arrays:
        if array.family_id != family_id or array.chrom != chrom:
            continue
        if min(end, array.end) > max(start, array.start):
            return True
    return False


def write_probes_fasta(path: Path, candidates: Sequence[ProbeCandidate]) -> None:
    lines = []
    for candidate in candidates:
        lines.append(
            (
                f">probe_id={candidate.probe_id};family_id={candidate.family_id};"
                f"length_bp={candidate.sequence_length};probe_score={candidate.probe_score:.4f}"
            )
        )
        lines.append(wrap_sequence(candidate.sequence))
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_probe_ranks(path: Path, candidates: Sequence[ProbeCandidate]) -> None:
    lines = [
        (
            "probe_id\tfamily_id\tsequence_length\tgc_content\ttm\testimated_copy_number\t"
            "arrayiness_score\tspecificity_score\toff_target_hits\tpredicted_regions\tprobe_score\twarning"
        )
    ]
    for candidate in candidates:
        lines.append(
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
                    candidate.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_fish_signals(path: Path, signals: Sequence[FishSignal]) -> None:
    lines = ["probe_id\tchrom\tstart\tend\tpredicted_signal\tconfidence"]
    for signal in signals:
        lines.append(
            "\t".join(
                [
                    signal.probe_id,
                    signal.chrom,
                    str(signal.start),
                    str(signal.end),
                    f"{signal.predicted_signal:.4f}",
                    signal.confidence,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
