"""Toy-scale read-based copy-number quantification MVP."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

from tandemx.discover.mvp import FastaRecord, read_fasta
from tandemx.simulate.toy import reverse_complement


@dataclass(frozen=True)
class MonomerRecord:
    family_id: str
    sequence: str


@dataclass(frozen=True)
class QuantifyConfig:
    reads: Path
    monomers: Path
    genome_size: int
    outdir: Path
    k: int
    haploid_depth: float | None


@dataclass(frozen=True)
class CopyNumberEstimate:
    family_id: str
    monomer_length: int
    diagnostic_kmer_count: int
    median_kmer_depth: float
    haploid_depth: float
    estimated_copy_number: float
    estimated_bp: float
    confidence: str
    warning: str


def quantify_toy_copy_number(config: QuantifyConfig) -> list[CopyNumberEstimate]:
    validate_quantify_config(config)
    monomers = list(read_monomer_fasta(config.monomers))
    if not monomers:
        raise ValueError("No monomers found for quantify")
    if all(len(monomer.sequence) < config.k for monomer in monomers):
        raise ValueError("--k is greater than all monomer lengths in the catalogue")

    read_kmers, total_read_bases, read_count, max_read_len = count_read_kmers_and_bases(config.reads, config.k)
    if read_count == 0:
        raise ValueError("No reads found for quantify")
    if max_read_len < config.k:
        raise ValueError("--k is greater than all read lengths")
    haploid_depth = (
        config.haploid_depth
        if config.haploid_depth is not None
        else total_read_bases / config.genome_size
    )
    shared_map = family_kmer_membership(monomers, config.k)

    estimates = []
    for monomer in monomers:
        monomer_counts = monomer_kmer_counts(monomer.sequence, config.k)
        diagnostic = {
            kmer: multiplicity
            for kmer, multiplicity in monomer_counts.items()
            if len(shared_map[kmer]) == 1 and not is_low_complexity_kmer(kmer)
        }
        corrected_depths = [
            read_kmers.get(kmer, 0) / multiplicity
            for kmer, multiplicity in diagnostic.items()
            if multiplicity > 0
        ]
        median_depth = float(median(corrected_depths)) if corrected_depths else 0.0
        estimated_copy_number = median_depth / haploid_depth if haploid_depth > 0 else 0.0
        estimated_bp = estimated_copy_number * len(monomer.sequence)
        warning_parts = []
        if config.haploid_depth is None:
            warning_parts.append("haploid_depth_estimated_from_total_read_bases_and_genome_size")
        if not diagnostic:
            warning_parts.append("no_diagnostic_kmers")
        confidence = "high"
        if len(diagnostic) < 10 or config.haploid_depth is None:
            confidence = "medium"
        if not diagnostic or haploid_depth <= 0:
            confidence = "low"
        estimates.append(
            CopyNumberEstimate(
                family_id=monomer.family_id,
                monomer_length=len(monomer.sequence),
                diagnostic_kmer_count=len(diagnostic),
                median_kmer_depth=median_depth,
                haploid_depth=haploid_depth,
                estimated_copy_number=estimated_copy_number,
                estimated_bp=estimated_bp,
                confidence=confidence,
                warning=";".join(warning_parts),
            )
        )

    write_copy_number(config.outdir / "copy_number.tsv", estimates)
    return estimates


def validate_quantify_config(config: QuantifyConfig) -> None:
    if config.genome_size <= 0:
        raise ValueError("--genome-size must be positive")
    if config.k <= 0:
        raise ValueError("--k must be positive")
    if config.haploid_depth is not None and config.haploid_depth <= 0:
        raise ValueError("--haploid-depth must be positive when provided")


def read_monomer_fasta(path: Path) -> Iterable[MonomerRecord]:
    for record in read_fasta(path):
        family_id = parse_family_id(record.description)
        yield MonomerRecord(family_id=family_id, sequence=record.sequence.replace("N", ""))


def parse_family_id(header: str) -> str:
    for part in header.split(";"):
        if part.startswith("family_id="):
            return part.split("=", 1)[1]
    return header.split()[0]


def canonical_kmer(kmer: str) -> str:
    reverse = reverse_complement(kmer)
    return min(kmer.upper(), reverse)


def iter_kmers(sequence: str, k: int) -> Iterable[str]:
    sequence = sequence.upper()
    for index in range(0, len(sequence) - k + 1):
        kmer = sequence[index : index + k]
        if "N" not in kmer:
            yield canonical_kmer(kmer)


def count_read_kmers(reads: Sequence[FastaRecord], k: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    for read in reads:
        counts.update(iter_kmers(read.sequence, k))
    return counts


def count_read_kmers_and_bases(path: Path, k: int) -> tuple[Counter[str], int, int, int]:
    counts: Counter[str] = Counter()
    total_bases = 0
    read_count = 0
    max_read_len = 0
    for read in read_fasta(path):
        read_count += 1
        total_bases += len(read.sequence)
        max_read_len = max(max_read_len, len(read.sequence))
        counts.update(iter_kmers(read.sequence, k))
    return counts, total_bases, read_count, max_read_len


def monomer_kmer_counts(sequence: str, k: int) -> Counter[str]:
    return Counter(iter_kmers(sequence, k))


def family_kmer_membership(monomers: Sequence[MonomerRecord], k: int) -> dict[str, set[str]]:
    membership: dict[str, set[str]] = defaultdict(set)
    for monomer in monomers:
        for kmer in monomer_kmer_counts(monomer.sequence, k):
            membership[kmer].add(monomer.family_id)
    return membership


def is_low_complexity_kmer(kmer: str) -> bool:
    if not kmer:
        return True
    counts = Counter(kmer)
    if max(counts.values()) / len(kmer) >= 0.8:
        return True
    dinucleotides = {kmer[index : index + 2] for index in range(0, len(kmer) - 1, 2)}
    return len(dinucleotides) <= 1


def write_copy_number(path: Path, estimates: Sequence[CopyNumberEstimate]) -> None:
    lines = [
        (
            "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
            "haploid_depth\testimated_copy_number\testimated_bp\tconfidence\twarning"
        )
    ]
    for estimate in estimates:
        lines.append(
            "\t".join(
                [
                    estimate.family_id,
                    str(estimate.monomer_length),
                    str(estimate.diagnostic_kmer_count),
                    f"{estimate.median_kmer_depth:.4f}",
                    f"{estimate.haploid_depth:.4f}",
                    f"{estimate.estimated_copy_number:.4f}",
                    f"{estimate.estimated_bp:.4f}",
                    estimate.confidence,
                    estimate.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
