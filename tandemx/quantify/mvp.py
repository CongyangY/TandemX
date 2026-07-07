"""Toy-scale read-based copy-number quantification MVP."""

from __future__ import annotations

import logging
import os
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

from tandemx.discover.mvp import FastaRecord, read_fasta, read_fasta_many
from tandemx.io.sequences import normalize_sequence_paths
from tandemx.discover.rust_backend import RustDiagnosticKmerCounter
from tandemx.simulate.toy import reverse_complement
from tandemx.utils.threads import discover_thread_limit
from tandemx.utils.progress import ProgressSnapshot, TerminalProgress


@dataclass(frozen=True)
class MonomerRecord:
    family_id: str
    sequence: str


@dataclass(frozen=True)
class QuantifyConfig:
    reads: Path | Sequence[Path]
    monomers: Path
    genome_size: int
    outdir: Path
    k: int
    haploid_depth: float | None
    kmer_backend: str = "python"
    max_reads: int | None = None
    max_read_bases: int | None = None
    progress_every: int = 1000


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


def quantify_toy_copy_number(
    config: QuantifyConfig,
    logger: logging.Logger | None = None,
    progress: TerminalProgress | None = None,
) -> list[CopyNumberEstimate]:
    validate_quantify_config(config)
    logger = logger or logging.getLogger("tandemx.quantify")
    update_quantify_terminal_progress(progress, "load_catalog", 0, 0, config)
    monomers = list(read_monomer_fasta(config.monomers))
    if not monomers:
        raise ValueError("No monomers found for quantify")
    if all(len(monomer.sequence) < config.k for monomer in monomers):
        raise ValueError("--k is greater than all monomer lengths in the catalogue")

    update_quantify_terminal_progress(progress, "build_diagnostic_kmers", 0, 0, config)
    shared_map = family_kmer_membership(monomers, config.k)
    diagnostic_by_family = {
        monomer.family_id: {
            kmer: multiplicity
            for kmer, multiplicity in monomer_kmer_counts(monomer.sequence, config.k).items()
            if len(shared_map[kmer]) == 1 and not is_low_complexity_kmer(kmer)
        }
        for monomer in monomers
    }
    target_kmers = {
        kmer
        for diagnostic in diagnostic_by_family.values()
        for kmer in diagnostic
    }
    read_kmers, total_read_bases, read_count, max_read_len = count_selected_read_kmers_and_bases(
        config.reads,
        config.k,
        target_kmers,
        config.kmer_backend,
        max_reads=config.max_reads,
        max_read_bases=config.max_read_bases,
        progress_every=config.progress_every,
        logger=logger,
        progress=progress,
    )
    if read_count == 0:
        raise ValueError("No reads found for quantify")
    if max_read_len < config.k:
        raise ValueError("--k is greater than all read lengths")
    haploid_depth = (
        config.haploid_depth
        if config.haploid_depth is not None
        else total_read_bases / config.genome_size
    )

    update_quantify_terminal_progress(progress, "estimate_copy_number", read_count, total_read_bases, config)
    estimates = []
    for monomer in monomers:
        diagnostic = diagnostic_by_family[monomer.family_id]
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

    update_quantify_terminal_progress(progress, "write_outputs", read_count, total_read_bases, config)
    write_copy_number(config.outdir / "copy_number.tsv", estimates)
    return estimates


def validate_quantify_config(config: QuantifyConfig) -> None:
    if config.genome_size <= 0:
        raise ValueError("--genome-size must be positive")
    if config.k <= 0:
        raise ValueError("--k must be positive")
    if config.haploid_depth is not None and config.haploid_depth <= 0:
        raise ValueError("--haploid-depth must be positive when provided")
    if config.kmer_backend not in {"python", "rust"}:
        raise ValueError("--kmer-backend must be python or rust")
    if config.kmer_backend == "rust" and config.k > 31:
        raise ValueError("Rust backend requires --k at most 31")
    if config.max_reads is not None and config.max_reads <= 0:
        raise ValueError("--max-reads must be positive when provided")
    if config.max_read_bases is not None and config.max_read_bases <= 0:
        raise ValueError("--max-read-bases must be positive when provided")
    if config.progress_every <= 0:
        raise ValueError("--progress-every must be positive")


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


def count_selected_read_kmers_and_bases(
    path: Path | Sequence[Path],
    k: int,
    targets: set[str],
    backend: str,
    *,
    max_reads: int | None = None,
    max_read_bases: int | None = None,
    progress_every: int = 1000,
    logger: logging.Logger | None = None,
    progress: TerminalProgress | None = None,
) -> tuple[Counter[str], int, int, int]:
    sequence_paths = normalize_sequence_paths(path)
    if (
        len(sequence_paths) > 1
        and max_reads is None
        and max_read_bases is None
    ):
        return count_selected_read_kmers_and_bases_parallel_files(
            sequence_paths,
            k,
            targets,
            backend,
            progress_every=progress_every,
            logger=logger,
            progress=progress,
        )

    counts: Counter[str] = Counter()
    rust_counter = RustDiagnosticKmerCounter(k, targets) if backend == "rust" else None
    total_bases = 0
    read_count = 0
    max_read_len = 0
    started = time.perf_counter()
    logger = logger or logging.getLogger("tandemx.quantify")
    update_quantify_read_progress(
        progress,
        read_count,
        total_bases,
        max_reads,
        max_read_bases,
    )
    for read in read_fasta_many(sequence_paths):
        if max_reads is not None and read_count >= max_reads:
            break
        if max_read_bases is not None and total_bases + len(read.sequence) > max_read_bases:
            logger.info(
                "limit_reached=max_read_bases configured_bases=%s next_read_bases=%s",
                max_read_bases,
                len(read.sequence),
            )
            break
        read_count += 1
        total_bases += len(read.sequence)
        max_read_len = max(max_read_len, len(read.sequence))
        if rust_counter is not None:
            rust_counter.count_sequence(read.sequence)
        else:
            counts.update(kmer for kmer in iter_kmers(read.sequence, k) if kmer in targets)
        if read_count % progress_every == 0:
            log_quantify_progress(logger, read_count, total_bases, started, max_reads, max_read_bases)
            update_quantify_read_progress(
                progress,
                read_count,
                total_bases,
                max_reads,
                max_read_bases,
            )
    log_quantify_progress(logger, read_count, total_bases, started, max_reads, max_read_bases)
    update_quantify_read_progress(
        progress,
        read_count,
        total_bases,
        max_reads,
        max_read_bases,
    )
    if rust_counter is not None:
        counts.update(rust_counter.counts())
    return counts, total_bases, read_count, max_read_len


def count_selected_read_kmers_and_bases_parallel_files(
    paths: Sequence[Path],
    k: int,
    targets: set[str],
    backend: str,
    *,
    progress_every: int = 1000,
    logger: logging.Logger | None = None,
    progress: TerminalProgress | None = None,
) -> tuple[Counter[str], int, int, int]:
    logger = logger or logging.getLogger("tandemx.quantify")
    started = time.perf_counter()
    workers = min(len(paths), max(1, min(discover_thread_limit(), os.cpu_count() or 1)))
    logger.info(
        "parallel_file_count enabled=true read_files=%s workers=%s backend=%s",
        len(paths),
        workers,
        backend,
    )
    update_quantify_read_progress(progress, 0, 0, None, None)
    partial_results: dict[Path, tuple[Counter[str], int, int, int]] = {}
    completed_reads = 0
    completed_bases = 0
    completed_max_read_len = 0
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_path = {
            executor.submit(
                count_selected_read_kmers_and_bases_one_file,
                path,
                k,
                targets,
                backend,
            ): path
            for path in paths
        }
        for completed_index, future in enumerate(as_completed(future_to_path), start=1):
            path = future_to_path[future]
            file_counts, file_bases, file_reads, file_max_len = future.result()
            partial_results[path] = (file_counts, file_bases, file_reads, file_max_len)
            completed_reads += file_reads
            completed_bases += file_bases
            completed_max_read_len = max(completed_max_read_len, file_max_len)
            if completed_index == len(paths) or completed_reads % progress_every == 0:
                log_quantify_progress(logger, completed_reads, completed_bases, started, None, None)
                update_quantify_read_progress(progress, completed_reads, completed_bases, None, None)

    merged_counts: Counter[str] = Counter()
    total_bases = 0
    read_count = 0
    max_read_len = 0
    for path in paths:
        file_counts, file_bases, file_reads, file_max_len = partial_results[path]
        merged_counts.update(file_counts)
        total_bases += file_bases
        read_count += file_reads
        max_read_len = max(max_read_len, file_max_len)
    return merged_counts, total_bases, read_count, max_read_len


def count_selected_read_kmers_and_bases_one_file(
    path: Path,
    k: int,
    targets: set[str],
    backend: str,
) -> tuple[Counter[str], int, int, int]:
    counts: Counter[str] = Counter()
    rust_counter = RustDiagnosticKmerCounter(k, targets) if backend == "rust" else None
    total_bases = 0
    read_count = 0
    max_read_len = 0
    for read in read_fasta(path):
        read_count += 1
        total_bases += len(read.sequence)
        max_read_len = max(max_read_len, len(read.sequence))
        if rust_counter is not None:
            rust_counter.count_sequence(read.sequence)
        else:
            counts.update(kmer for kmer in iter_kmers(read.sequence, k) if kmer in targets)
    if rust_counter is not None:
        counts.update(rust_counter.counts())
    return counts, total_bases, read_count, max_read_len


def update_quantify_terminal_progress(
    progress: TerminalProgress | None,
    step: str,
    processed_reads: int,
    processed_bases: int,
    config: QuantifyConfig,
) -> None:
    if progress is None:
        return
    progress.update(
        ProgressSnapshot(
            command="quantify",
            step=step,
            processed_reads=processed_reads,
            processed_bases=processed_bases,
            total_reads=config.max_reads,
            total_bases=config.max_read_bases,
        )
    )


def update_quantify_read_progress(
    progress: TerminalProgress | None,
    processed_reads: int,
    processed_bases: int,
    max_reads: int | None,
    max_read_bases: int | None,
) -> None:
    if progress is None:
        return
    progress.update(
        ProgressSnapshot(
            command="quantify",
            step="scan_reads",
            processed_reads=processed_reads,
            processed_bases=processed_bases,
            total_reads=max_reads,
            total_bases=max_read_bases,
        )
    )


def log_quantify_progress(
    logger: logging.Logger,
    processed_reads: int,
    processed_bases: int,
    started: float,
    max_reads: int | None,
    max_read_bases: int | None,
) -> None:
    elapsed = max(time.perf_counter() - started, 1e-9)
    reads_per_second = processed_reads / elapsed
    mb_per_second = (processed_bases / 1_000_000) / elapsed
    remaining_estimates: list[float] = []
    if max_reads is not None and reads_per_second > 0:
        remaining_estimates.append(max(0, max_reads - processed_reads) / reads_per_second)
    if max_read_bases is not None and mb_per_second > 0:
        remaining_mb = max(0, max_read_bases - processed_bases) / 1_000_000
        remaining_estimates.append(remaining_mb / mb_per_second)
    estimated_remaining = f"{min(remaining_estimates):.1f}" if remaining_estimates else "unknown"
    logger.info(
        "progress processed_reads=%s processed_bases=%s elapsed_seconds=%.3f "
        "reads_per_second=%.3f mb_per_second=%.3f estimated_remaining_seconds=%s",
        processed_reads,
        processed_bases,
        elapsed,
        reads_per_second,
        mb_per_second,
        estimated_remaining,
    )


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
