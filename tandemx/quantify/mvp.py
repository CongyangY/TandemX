"""Toy-scale read-based copy-number quantification MVP."""

from __future__ import annotations

import logging
import math
import os
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

from tandemx.discover.mvp import FastaRecord, read_fasta, read_fasta_many
from tandemx.io.sequences import (
    detect_sequence_format,
    normalize_sequence_paths,
    read_sequence_records_many,
)
from tandemx.discover.rust_backend import RustDiagnosticKmerCounter
from tandemx.utils.threads import discover_thread_limit
from tandemx.utils.progress import ProgressSnapshot, TerminalProgress
from tandemx.utils.kmers import (
    canonical_kmer,
    canonical_kmer_code,
    circular_kmer_counts,
    is_low_complexity_kmer,
    iter_canonical_kmer_codes,
    iter_linear_canonical_kmers,
)


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
    single_copy_kmers: Path | None = None
    single_copy_min_depth: float | None = None
    read_error_rate: float | None = None
    quality_correction_enabled: bool = True


@dataclass(frozen=True)
class CopyNumberEstimate:
    family_id: str
    monomer_length: int
    diagnostic_kmer_count: int
    median_kmer_depth: float
    haploid_depth: float
    estimated_copy_number: float
    estimated_bp: float
    depth_mad: float
    copy_number_interval_low: float
    copy_number_interval_high: float
    normalization_method: str
    raw_median_kmer_depth: float
    kmer_survival_probability: float
    single_copy_control_kmer_count: int
    single_copy_control_mean_depth: float | None
    single_copy_control_median_depth: float | None
    single_copy_control_depth_mad: float | None
    single_copy_control_zero_fraction: float | None
    quality_window_count: int
    confidence: str
    warning: str


@dataclass(frozen=True)
class KmerSurvivalStats:
    survival_sum: float = 0.0
    window_count: int = 0
    quality_window_count: int = 0
    assumed_error_window_count: int = 0

    @property
    def mean_survival(self) -> float:
        return self.survival_sum / self.window_count if self.window_count else 1.0


@dataclass(frozen=True)
class ControlDepthStats:
    mean: float
    median: float
    mad: float
    zero_fraction: float


def select_haploid_depth(
    config: QuantifyConfig,
    control_stats: ControlDepthStats | None,
    total_read_bases: int,
) -> tuple[float, str]:
    """Select a recorded normalization path, including an optional depth gate."""
    if config.haploid_depth is not None:
        return config.haploid_depth, "explicit_haploid_depth"
    if control_stats is not None:
        if (
            config.single_copy_min_depth is not None
            and control_stats.mean < config.single_copy_min_depth
        ):
            return (
                total_read_bases / config.genome_size,
                "total_read_bases_divided_by_genome_size_low_control_depth_fallback",
            )
        if control_stats.mean <= 0:
            raise ValueError("Single-copy controls have zero mean observed depth")
        return control_stats.mean, "empirical_single_copy_kmers_mean"
    return total_read_bases / config.genome_size, "total_read_bases_divided_by_genome_size"


def quantify_toy_copy_number(
    config: QuantifyConfig,
    logger: logging.Logger | None = None,
    progress: TerminalProgress | None = None,
) -> list[CopyNumberEstimate]:
    validate_quantify_config(config)
    config.outdir.mkdir(parents=True, exist_ok=True)
    logger = logger or logging.getLogger("tandemx.quantify")
    update_quantify_terminal_progress(progress, "load_catalog", 0, 0, config)
    monomers = list(read_monomer_fasta(config.monomers))
    if not monomers:
        raise ValueError("No monomers found for quantify")
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
    single_copy_controls = (
        read_single_copy_kmers(config.single_copy_kmers, config.k)
        if config.single_copy_kmers is not None
        else {}
    )
    target_kmers = {
        kmer
        for diagnostic in diagnostic_by_family.values()
        for kmer in diagnostic
    }
    overlap = target_kmers.intersection(single_copy_controls)
    if overlap:
        raise ValueError(
            f"Single-copy control k-mers overlap repeat diagnostic k-mers: {sorted(overlap)[:3]}"
        )
    target_kmers.update(single_copy_controls)
    sequence_paths = normalize_sequence_paths(config.reads)
    uses_fastq = any(detect_sequence_format(path) == "fastq" for path in sequence_paths)
    if config.quality_correction_enabled and (uses_fastq or config.read_error_rate is not None):
        read_kmers, total_read_bases, read_count, max_read_len, survival = (
            count_selected_read_kmers_quality_and_bases(
                sequence_paths,
                config.k,
                target_kmers,
                config.kmer_backend,
                max_reads=config.max_reads,
                max_read_bases=config.max_read_bases,
                progress_every=config.progress_every,
                logger=logger,
                progress=progress,
                assumed_error_rate=config.read_error_rate,
            )
        )
    else:
        read_kmers, total_read_bases, read_count, max_read_len = count_selected_read_kmers_and_bases(
            sequence_paths,
            config.k,
            target_kmers,
            config.kmer_backend,
            max_reads=config.max_reads,
            max_read_bases=config.max_read_bases,
            progress_every=config.progress_every,
            logger=logger,
            progress=progress,
        )
        survival = KmerSurvivalStats()
    if read_count == 0:
        raise ValueError("No reads found for quantify")
    if max_read_len < config.k:
        raise ValueError("--k is greater than all read lengths")
    survival_probability = survival.mean_survival
    if (
        config.quality_correction_enabled
        and (uses_fastq or config.read_error_rate is not None)
        and survival.window_count == 0
    ):
        raise ValueError("No valid ACGT k-mer windows available for survival correction")
    if not 0 < survival_probability <= 1:
        raise ValueError("Effective k-mer survival probability must be in (0,1]")
    corrected_read_kmers = {
        kmer: count / survival_probability for kmer, count in read_kmers.items()
    }
    control_depths = [
        corrected_read_kmers.get(kmer, 0.0) / expected_copy_number
        for kmer, expected_copy_number in single_copy_controls.items()
    ]
    control_stats = estimate_control_depth(control_depths) if control_depths else None
    # Include zero-observation controls. Their arithmetic mean is an unbiased
    # depth estimator under uniform independent sampling, whereas the median
    # becomes exactly zero around or below 1x coverage. The optional gate
    # preserves the observed control statistics while selecting total-bases
    # normalization when their mean support is below a declared threshold.
    haploid_depth, normalization_method = select_haploid_depth(
        config, control_stats, total_read_bases
    )

    update_quantify_terminal_progress(progress, "estimate_copy_number", read_count, total_read_bases, config)
    estimates = []
    for monomer in monomers:
        diagnostic = diagnostic_by_family[monomer.family_id]
        raw_depths = [
            read_kmers.get(kmer, 0) / multiplicity
            for kmer, multiplicity in diagnostic.items()
            if multiplicity > 0
        ]
        corrected_depths = [depth / survival_probability for depth in raw_depths]
        raw_median_depth = float(median(raw_depths)) if raw_depths else 0.0
        median_depth = float(median(corrected_depths)) if corrected_depths else 0.0
        depth_mad = (
            float(median(abs(value - median_depth) for value in corrected_depths))
            if corrected_depths
            else 0.0
        )
        estimated_copy_number = median_depth / haploid_depth if haploid_depth > 0 else 0.0
        estimated_bp = estimated_copy_number * len(monomer.sequence)
        interval_low_depth = empirical_quantile(corrected_depths, 0.10)
        interval_high_depth = empirical_quantile(corrected_depths, 0.90)
        interval_low = interval_low_depth / haploid_depth if haploid_depth > 0 else 0.0
        interval_high = interval_high_depth / haploid_depth if haploid_depth > 0 else 0.0
        warning_parts = []
        if normalization_method.startswith("total_read_bases_divided_by_genome_size"):
            warning_parts.append("haploid_depth_estimated_from_total_read_bases_and_genome_size")
            if normalization_method.endswith("low_control_depth_fallback"):
                warning_parts.append("single_copy_control_mean_depth_below_configured_minimum")
                warning_parts.append("single_copy_controls_reported_but_not_used")
        elif normalization_method == "empirical_single_copy_kmers_mean":
            warning_parts.append("single_copy_status_depends_on_user_supplied_control_provenance")
            if len(single_copy_controls) < 100:
                warning_parts.append("fewer_than_100_single_copy_control_kmers")
            if control_stats is not None and control_stats.zero_fraction > 0.5:
                warning_parts.append("more_than_half_single_copy_controls_unobserved")
        elif single_copy_controls:
            warning_parts.append("single_copy_controls_reported_but_explicit_depth_used")
        if survival_probability < 1:
            warning_parts.append("independent_base_error_survival_approximation")
            if survival.quality_window_count:
                warning_parts.append("phred_scores_treated_as_calibrated_probabilities")
            if (
                survival.assumed_error_window_count
                and config.read_error_rate is None
            ):
                warning_parts.append("fasta_windows_in_mixed_input_treated_as_error_free")
        if not diagnostic:
            warning_parts.append("no_diagnostic_kmers")
        else:
            warning_parts.append("genome_background_uniqueness_not_verified")
        # Catalogue-only diagnostic k-mers have not been checked against an
        # independent genome background, so this MVP must not label them high.
        confidence = "medium"
        relative_mad = depth_mad / median_depth if median_depth > 0 else 0.0
        if len(diagnostic) < 10 or config.haploid_depth is None or relative_mad > 0.25:
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
                depth_mad=depth_mad,
                copy_number_interval_low=interval_low,
                copy_number_interval_high=interval_high,
                normalization_method=normalization_method,
                raw_median_kmer_depth=raw_median_depth,
                kmer_survival_probability=survival_probability,
                single_copy_control_kmer_count=len(single_copy_controls),
                single_copy_control_mean_depth=control_stats.mean if control_stats else None,
                single_copy_control_median_depth=control_stats.median if control_stats else None,
                single_copy_control_depth_mad=control_stats.mad if control_stats else None,
                single_copy_control_zero_fraction=control_stats.zero_fraction if control_stats else None,
                quality_window_count=survival.quality_window_count,
                confidence=confidence,
                warning=";".join(warning_parts),
            )
        )

    update_quantify_terminal_progress(progress, "write_outputs", read_count, total_read_bases, config)
    write_copy_number(config.outdir / "copy_number.tsv", estimates)
    return estimates


def validate_quantify_config(config: QuantifyConfig) -> None:
    paths = normalize_sequence_paths(config.reads)
    if len(set(paths)) != len(paths):
        raise ValueError("--reads must not contain the same file path more than once")
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
    if config.single_copy_kmers is not None and not config.single_copy_kmers.is_file():
        raise ValueError("--single-copy-kmers must be an existing TSV file")
    if config.single_copy_min_depth is not None:
        if not math.isfinite(config.single_copy_min_depth) or config.single_copy_min_depth < 0:
            raise ValueError("--single-copy-min-depth must be finite and nonnegative")
        if config.single_copy_kmers is None:
            raise ValueError("--single-copy-min-depth requires --single-copy-kmers")
        if config.haploid_depth is not None:
            raise ValueError("--single-copy-min-depth cannot be combined with --haploid-depth")
    if config.read_error_rate is not None and not 0 <= config.read_error_rate < 1:
        raise ValueError("--read-error-rate must be in [0,1)")
    if type(config.quality_correction_enabled) is not bool:
        raise ValueError("quality correction flag must be boolean")
    if config.read_error_rate is not None and not config.quality_correction_enabled:
        raise ValueError(
            "--read-error-rate cannot be combined with --disable-quality-correction"
        )


def read_single_copy_kmers(path: Path, k: int) -> dict[str, float]:
    """Read canonical background controls with declared haploid copy number."""
    import csv

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Single-copy control TSV requires unique header names")
        if not {"kmer", "expected_copy_number"} <= set(reader.fieldnames):
            raise ValueError("Single-copy control TSV requires kmer and expected_copy_number columns")
        controls: dict[str, float] = {}
        for line_number, row in enumerate(reader, 2):
            raw_word = row.get("kmer")
            if not isinstance(raw_word, str):
                raise ValueError(f"Missing k-mer control at line {line_number}")
            word = raw_word.strip().upper()
            if len(word) != k or set(word) - set("ACGT"):
                raise ValueError(f"Invalid {k}-mer control at line {line_number}")
            if is_low_complexity_kmer(word):
                raise ValueError(f"Low-complexity control k-mer at line {line_number}")
            word = canonical_kmer(word)
            try:
                expected = float(row.get("expected_copy_number", ""))
            except (TypeError, ValueError) as error:
                raise ValueError(f"Invalid expected copy number at line {line_number}") from error
            if not math.isfinite(expected) or expected <= 0:
                raise ValueError(f"Expected copy number must be positive at line {line_number}")
            if word in controls:
                raise ValueError(f"Duplicate canonical control k-mer at line {line_number}: {word}")
            controls[word] = expected
    if not controls:
        raise ValueError("Single-copy control TSV contains no records")
    return controls


def estimate_control_depth(depths: Sequence[float]) -> ControlDepthStats:
    """Summarize normalized controls without dropping unobserved k-mers."""
    if not depths or any(not math.isfinite(value) or value < 0 for value in depths):
        raise ValueError("Control depths must be nonempty, finite and nonnegative")
    center = float(median(depths))
    return ControlDepthStats(
        mean=math.fsum(depths) / len(depths),
        median=center,
        mad=float(median(abs(value - center) for value in depths)),
        zero_fraction=sum(value == 0 for value in depths) / len(depths),
    )


def read_monomer_fasta(path: Path) -> Iterable[MonomerRecord]:
    for record in read_fasta(path):
        family_id = parse_family_id(record.description)
        yield MonomerRecord(family_id=family_id, sequence=record.sequence)


def parse_family_id(header: str) -> str:
    for part in header.split(";"):
        if part.startswith("family_id="):
            return part.split("=", 1)[1]
    return header.split()[0]


def iter_kmers(sequence: str, k: int) -> Iterable[str]:
    yield from iter_linear_canonical_kmers(sequence, k)


def selected_target_code_map(targets: set[str], k: int) -> dict[int, str] | None:
    """Map exactly matchable canonical target strings to rolling 2-bit codes.

    ``None`` preserves the legacy string iterator for k values and text that the
    encoder cannot represent.  Invalid/noncanonical target strings never match
    legacy canonical output, so omitting them is equivalent.
    """
    if not 1 <= k <= 31:
        return None
    result: dict[int, str] = {}
    for target in targets:
        if not isinstance(target, str) or len(target) != k or not target.isascii():
            continue
        normalized = target.upper()
        if set(normalized) - set("ACGT") or canonical_kmer(normalized) != target:
            continue
        result[canonical_kmer_code(normalized)] = target
    return result


def update_selected_python_kmer_counts(
    counts: Counter[str], sequence: str, k: int, targets: set[str], code_targets: dict[int, str] | None,
) -> None:
    """Count the same selected canonical words as the historical string path."""
    if code_targets is None or not sequence.isascii():
        counts.update(kmer for kmer in iter_kmers(sequence, k) if kmer in targets)
        return
    counts.update(code_targets[code] for _, code in iter_canonical_kmer_codes(sequence, k)
                  if code in code_targets)


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
    python_target_codes = None if rust_counter is not None else selected_target_code_map(targets, k)
    total_bases = 0
    read_count = 0
    max_read_len = 0
    started = time.perf_counter()
    rust_batch: list[str] = []
    rust_batch_bases = 0
    logger = logger or logging.getLogger("tandemx.quantify")
    update_quantify_read_progress(
        progress,
        read_count,
        total_bases,
        max_reads,
        max_read_bases,
    )
    for read in read_fasta_many(
        sequence_paths,
        check_duplicate_ids_across_files=False,
    ):
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
            rust_batch.append(read.sequence)
            rust_batch_bases += len(read.sequence)
            if len(rust_batch) >= 512 or rust_batch_bases >= 8_000_000:
                rust_counter.count_sequences(rust_batch)
                rust_batch.clear()
                rust_batch_bases = 0
        else:
            update_selected_python_kmer_counts(counts, read.sequence, k, targets, python_target_codes)
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
        rust_counter.count_sequences(rust_batch)
        counts.update(rust_counter.counts())
    return counts, total_bases, read_count, max_read_len


def quality_window_survival(
    sequence: str,
    quality: str | None,
    k: int,
    assumed_error_rate: float | None,
) -> tuple[float, int, int, int]:
    """Return survival sum, valid windows, quality windows and assumed windows."""
    sequence = sequence.upper()
    if quality is not None and len(quality) != len(sequence):
        raise ValueError("FASTQ sequence and quality lengths must match")
    if assumed_error_rate is not None and not 0 <= assumed_error_rate < 1:
        raise ValueError("Assumed error rate must be in [0,1)")
    if quality is not None and any(not 0 <= ord(character) - 33 <= 93 for character in quality):
        raise ValueError("FASTQ quality characters must encode Phred+33 values in 0..93")
    if quality is None:
        if not set(sequence).difference("ACGT"):
            valid_windows = max(0, len(sequence) - k + 1)
        else:
            valid_windows = 0
            run_length = 0
            for character in sequence:
                if character in "ACGT":
                    run_length += 1
                else:
                    valid_windows += max(0, run_length - k + 1)
                    run_length = 0
            valid_windows += max(0, run_length - k + 1)
        rate = assumed_error_rate or 0.0
        return (
            valid_windows * (1.0 - rate) ** k,
            valid_windows,
            0,
            valid_windows,
        )
    survival_sum = 0.0
    valid_windows = 0
    quality_windows = 0
    assumed_windows = 0
    start = 0
    while start < len(sequence):
        while start < len(sequence) and sequence[start] not in "ACGT":
            start += 1
        end = start
        while end < len(sequence) and sequence[end] in "ACGT":
            end += 1
        length = end - start
        windows = max(0, length - k + 1)
        if windows:
            valid_windows += windows
            if quality is None:
                rate = assumed_error_rate or 0.0
                survival_sum += windows * (1.0 - rate) ** k
                assumed_windows += windows
            else:
                logs: list[float] = []
                zeros: list[int] = []
                for character in quality[start:end]:
                    phred = ord(character) - 33
                    probability = 1.0 - 10.0 ** (-phred / 10.0)
                    logs.append(math.log(probability) if probability > 0 else 0.0)
                    zeros.append(int(probability == 0))
                log_sum = sum(logs[:k])
                zero_count = sum(zeros[:k])
                survival_sum += 0.0 if zero_count else math.exp(log_sum)
                for index in range(k, length):
                    log_sum += logs[index] - logs[index - k]
                    zero_count += zeros[index] - zeros[index - k]
                    survival_sum += 0.0 if zero_count else math.exp(log_sum)
                quality_windows += windows
        start = end + 1
    return survival_sum, valid_windows, quality_windows, assumed_windows


def count_selected_read_kmers_quality_and_bases(
    paths: Sequence[Path],
    k: int,
    targets: set[str],
    backend: str,
    *,
    max_reads: int | None,
    max_read_bases: int | None,
    progress_every: int,
    logger: logging.Logger,
    progress: TerminalProgress | None,
    assumed_error_rate: float | None,
) -> tuple[Counter[str], int, int, int, KmerSurvivalStats]:
    """Count target k-mers and their global error survival in one streaming pass."""
    counts: Counter[str] = Counter()
    rust_counter = RustDiagnosticKmerCounter(k, targets) if backend == "rust" else None
    python_target_codes = None if rust_counter is not None else selected_target_code_map(targets, k)
    total_bases = read_count = max_read_len = 0
    survival_sum = 0.0
    window_count = quality_windows = assumed_windows = 0
    rust_batch: list[str] = []
    rust_batch_bases = 0
    started = time.perf_counter()
    for record in read_sequence_records_many(
        paths, check_duplicate_ids_across_files=False
    ):
        if max_reads is not None and read_count >= max_reads:
            break
        if max_read_bases is not None and total_bases + len(record.sequence) > max_read_bases:
            logger.info(
                "limit_reached=max_read_bases configured_bases=%s next_read_bases=%s",
                max_read_bases,
                len(record.sequence),
            )
            break
        read_count += 1
        total_bases += len(record.sequence)
        max_read_len = max(max_read_len, len(record.sequence))
        observed = quality_window_survival(
            record.sequence, record.quality, k, assumed_error_rate
        )
        survival_sum += observed[0]
        window_count += observed[1]
        quality_windows += observed[2]
        assumed_windows += observed[3]
        if rust_counter is not None:
            rust_batch.append(record.sequence)
            rust_batch_bases += len(record.sequence)
            if len(rust_batch) >= 512 or rust_batch_bases >= 8_000_000:
                rust_counter.count_sequences(rust_batch)
                rust_batch.clear()
                rust_batch_bases = 0
        else:
            update_selected_python_kmer_counts(counts, record.sequence, k, targets, python_target_codes)
        if read_count % progress_every == 0:
            log_quantify_progress(
                logger, read_count, total_bases, started, max_reads, max_read_bases
            )
            update_quantify_read_progress(
                progress, read_count, total_bases, max_reads, max_read_bases
            )
    if rust_counter is not None:
        rust_counter.count_sequences(rust_batch)
        counts.update(rust_counter.counts())
    log_quantify_progress(
        logger, read_count, total_bases, started, max_reads, max_read_bases
    )
    update_quantify_read_progress(
        progress, read_count, total_bases, max_reads, max_read_bases
    )
    return counts, total_bases, read_count, max_read_len, KmerSurvivalStats(
        survival_sum, window_count, quality_windows, assumed_windows
    )


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
    python_target_codes = None if rust_counter is not None else selected_target_code_map(targets, k)
    total_bases = 0
    read_count = 0
    max_read_len = 0
    rust_batch: list[str] = []
    rust_batch_bases = 0
    for read in read_fasta(path):
        read_count += 1
        total_bases += len(read.sequence)
        max_read_len = max(max_read_len, len(read.sequence))
        if rust_counter is not None:
            rust_batch.append(read.sequence)
            rust_batch_bases += len(read.sequence)
            if len(rust_batch) >= 512 or rust_batch_bases >= 8_000_000:
                rust_counter.count_sequences(rust_batch)
                rust_batch.clear()
                rust_batch_bases = 0
        else:
            update_selected_python_kmer_counts(counts, read.sequence, k, targets, python_target_codes)
    if rust_counter is not None:
        rust_counter.count_sequences(rust_batch)
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
    return circular_kmer_counts(sequence, k)


def family_kmer_membership(monomers: Sequence[MonomerRecord], k: int) -> dict[str, set[str]]:
    membership: dict[str, set[str]] = defaultdict(set)
    for monomer in monomers:
        for kmer in monomer_kmer_counts(monomer.sequence, k):
            membership[kmer].add(monomer.family_id)
    return membership


def empirical_quantile(values: Sequence[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(len(ordered) - 1, lower + 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def write_copy_number(path: Path, estimates: Sequence[CopyNumberEstimate]) -> None:
    lines = [
        (
            "family_id\tmonomer_length\tdiagnostic_kmer_count\tmedian_kmer_depth\t"
            "haploid_depth\testimated_copy_number\testimated_bp\tdepth_mad\t"
            "copy_number_interval_low\tcopy_number_interval_high\t"
            "normalization_method\traw_median_kmer_depth\tkmer_survival_probability\t"
            "single_copy_control_kmer_count\tsingle_copy_control_mean_depth\t"
            "single_copy_control_median_depth\tsingle_copy_control_depth_mad\t"
            "single_copy_control_zero_fraction\t"
            "quality_window_count\tconfidence\twarning"
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
                    f"{estimate.depth_mad:.4f}",
                    f"{estimate.copy_number_interval_low:.4f}",
                    f"{estimate.copy_number_interval_high:.4f}",
                    estimate.normalization_method,
                    f"{estimate.raw_median_kmer_depth:.4f}",
                    f"{estimate.kmer_survival_probability:.8f}",
                    str(estimate.single_copy_control_kmer_count),
                    (
                        f"{estimate.single_copy_control_mean_depth:.4f}"
                        if estimate.single_copy_control_mean_depth is not None
                        else "NA"
                    ),
                    (
                        f"{estimate.single_copy_control_median_depth:.4f}"
                        if estimate.single_copy_control_median_depth is not None
                        else "NA"
                    ),
                    (
                        f"{estimate.single_copy_control_depth_mad:.4f}"
                        if estimate.single_copy_control_depth_mad is not None
                        else "NA"
                    ),
                    (
                        f"{estimate.single_copy_control_zero_fraction:.6f}"
                        if estimate.single_copy_control_zero_fraction is not None
                        else "NA"
                    ),
                    str(estimate.quality_window_count),
                    estimate.confidence,
                    estimate.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
