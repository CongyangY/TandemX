"""Toy-scale tandem repeat discovery MVP."""

from __future__ import annotations

from collections import Counter, deque
import hashlib
import logging
import random
import time
from contextlib import nullcontext
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

from tandemx.io.sequences import (
    DuplicateIdTracker,
    SequenceFormatError,
    SequenceStats,
    normalize_sequence_paths,
    read_sequence_records,
    read_sequence_records_many,
)
from tandemx.discover.spacing import (
    best_local_periodicity_score,
    build_spacing_histogram,
    canonical_kmer,
    extract_repeated_kmer_positions,
    is_low_complexity_kmer,
    modulo_periodicity_score,
    refine_candidate_period,
    refine_candidate_period_with_interval,
    select_candidate_periods,
)
from tandemx.discover.rust_backend import RustBackendUnavailable, scan_read_for_periods, scan_reads_for_periods
from tandemx.simulate.toy import reverse_complement, wrap_sequence
from tandemx.utils.progress import ProgressSnapshot, TerminalProgress
from tandemx.utils.kmers import iter_canonical_kmer_codes
from tandemx.utils.threads import effective_discover_threads


SHORT_PERIOD_SCAN_MAX = 19
SHORT_PERIOD_ACCEPTANCE_SCORE = 0.80
DEFAULT_AUTO_DISCOVERY_TRIGGER_BASES = 20_000_000_000
DEFAULT_AUTO_DISCOVERY_MAX_BASES = 40_000_000_000
DEFAULT_TARGET_DISCOVERY_COVERAGE = 10.0


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


@dataclass
class CandidateCluster:
    """Incremental candidate cluster statistics used by family clustering."""

    period_sum: int
    candidate_count: int
    read_ids: set[str]
    support_span_sum: int
    score_sum: float
    low_complexity_flag: bool
    representative: CandidateRepeat
    members: list[CandidateRepeat]
    representative_kmers: frozenset[int]

    @classmethod
    def from_candidate(
        cls,
        candidate: CandidateRepeat,
        candidate_kmers: frozenset[int] | None = None,
    ) -> CandidateCluster:
        return cls(
            period_sum=candidate.period_bp,
            candidate_count=1,
            read_ids={candidate.read_id},
            support_span_sum=candidate.repeat_span_bp,
            score_sum=candidate.score,
            low_complexity_flag=candidate.low_complexity_flag,
            representative=candidate,
            members=[candidate],
            representative_kmers=(
                candidate_kmers
                if candidate_kmers is not None
                else candidate_kmer_sketch(candidate.sequence)
            ),
        )

    @property
    def center(self) -> int:
        return round(self.period_sum / self.candidate_count)

    @property
    def support_read_count(self) -> int:
        return len(self.read_ids)

    @property
    def mean_identity(self) -> float:
        return self.score_sum / self.candidate_count

    def add(
        self,
        candidate: CandidateRepeat,
        candidate_kmers: frozenset[int] | None = None,
    ) -> None:
        self.period_sum += candidate.period_bp
        self.candidate_count += 1
        self.read_ids.add(candidate.read_id)
        self.support_span_sum += candidate.repeat_span_bp
        self.score_sum += candidate.score
        self.low_complexity_flag = self.low_complexity_flag or candidate.low_complexity_flag
        self.members.append(candidate)
        if (candidate.score, candidate.repeat_span_bp) > (
            self.representative.score,
            self.representative.repeat_span_bp,
        ):
            self.representative = candidate
            self.representative_kmers = (
                candidate_kmers
                if candidate_kmers is not None
                else candidate_kmer_sketch(candidate.sequence)
            )


@dataclass(frozen=True)
class FamilySimilarity:
    family_a: str
    family_b: str
    length_a_bp: int
    length_b_bp: int
    kmer_jaccard: float
    shared_kmer_fraction: float
    local_identity: float
    local_overlap_bp: int
    local_overlap_fraction_shorter: float
    length_ratio: float
    orientation: str
    relationship: str
    redundant_candidate: bool
    notes: str


@dataclass(frozen=True)
class FamilyCollapse:
    original_family_id: str
    retained_family_id: str
    action: str
    reason: str
    relationship: str
    similarity_metrics: str
    notes: str


@dataclass(frozen=True)
class DiscoverConfig:
    reads: Path | Sequence[Path]
    outdir: Path
    min_monomer_len: int
    max_monomer_len: int
    min_support_reads: int
    min_repeat_span: int
    min_read_length: int = 1
    kmer_size: int = 11
    top_periods: int = 5
    min_seed_occurrences: int = 2
    min_spacing_support: int = 2
    max_pairs_per_kmer: int = 100
    max_reads: int | None = None
    max_read_bases: int | None = None
    genome_size: int | None = None
    total_reads: int | None = None
    total_read_bases: int | None = None
    sample_rate: float = 1.0
    seed: int = 1
    progress_every: int = 1000
    chunk_size: int = 1000
    chunk_bases: int = 8_000_000
    threads: int = 8
    kmer_backend: str = "python"
    target_discovery_coverage: float = DEFAULT_TARGET_DISCOVERY_COVERAGE
    auto_discovery_trigger_bases: int = DEFAULT_AUTO_DISCOVERY_TRIGGER_BASES
    auto_discovery_max_bases: int = DEFAULT_AUTO_DISCOVERY_MAX_BASES
    enable_auto_discovery_budget: bool = False
    collapse_redundant_families: bool = False


@dataclass(frozen=True)
class ReadScanTask:
    record: FastaRecord


@dataclass(frozen=True)
class ReadScanResult:
    read_bases: int
    candidate: CandidateRepeat | None = None
    skipped_short_reads: int = 0
    skipped_short_kmer: int = 0
    skipped_low_complexity: int = 0
    seed_overflow_count: int = 0


@dataclass
class DiscoverProgressTotals:
    total_reads: int | None = None
    total_read_bases: int | None = None
    effective_max_reads: int | None = None
    effective_max_read_bases: int | None = None
    count_done: bool = False


def discover_toy_repeats(
    config: DiscoverConfig,
    logger: logging.Logger | None = None,
    progress: TerminalProgress | None = None,
    count_future: Future[SequenceStats] | None = None,
    progress_totals: DiscoverProgressTotals | None = None,
) -> tuple[list[CandidateRepeat], list[RepeatFamily]]:
    """Discover tandem repeat families with a bounded k-mer spacing prefilter."""
    validate_discover_config(config)
    config.outdir.mkdir(parents=True, exist_ok=True)
    logger = logger or logging.getLogger("tandemx.discover")
    candidate_path = config.outdir / "candidate_reads.tsv"
    candidate_tmp_path = candidate_path.with_suffix(candidate_path.suffix + ".tmp")
    for stale_path in (
        config.outdir / "monomers.fa",
        config.outdir / "families.tsv",
        config.outdir / "family_similarity.tsv",
        config.outdir / "collapsed_families.tsv",
        config.outdir / "collapsed_monomers.fa",
        config.outdir / "family_collapse.tsv",
    ):
        stale_path.unlink(missing_ok=True)

    candidates: list[CandidateRepeat] = []
    processed_reads = 0
    processed_bases = 0
    skipped_short_reads = 0
    skipped_short_kmer = 0
    skipped_low_complexity = 0
    seed_overflow_count = 0
    started = time.perf_counter()
    rng = random.Random(config.seed)
    sequence_paths = normalize_sequence_paths(config.reads)
    effective_max_reads = config.max_reads
    effective_max_read_bases = config.max_read_bases
    totals = progress_totals or DiscoverProgressTotals(
        total_reads=config.total_reads,
        total_read_bases=config.total_read_bases,
        effective_max_reads=effective_max_reads,
        effective_max_read_bases=effective_max_read_bases,
        count_done=config.total_reads is not None or config.total_read_bases is not None,
    )
    effective_threads = effective_discover_threads(config.threads)
    parallel_scan = config.kmer_backend == "rust" and effective_threads > 1
    parallel_file_scan = False
    effective_max_reads, effective_max_read_bases = maybe_enable_auto_discovery_budget(
        config,
        totals,
        logger,
        effective_max_reads,
        effective_max_read_bases,
    )
    logger.info(
        "scan_threads requested=%s effective=%s parallel=%s parallel_files=%s backend=%s",
        config.threads,
        effective_threads,
        str(parallel_scan).lower(),
        str(parallel_file_scan).lower(),
        config.kmer_backend,
    )
    if effective_threads < config.threads:
        logger.info(
            "threads_capped requested=%s effective=%s reason=host_thread_policy",
            config.threads,
            effective_threads,
        )
    if config.kmer_backend != "rust" and effective_threads > 1:
        logger.info("parallel_scan_disabled reason=python_backend_gil effective_scan_threads=1")
    if progress is not None:
        poll_discover_count_future(count_future, totals, logger, config)
        effective_max_reads, effective_max_read_bases = maybe_enable_auto_discovery_budget(
            config,
            totals,
            logger,
            effective_max_reads,
            effective_max_read_bases,
        )
        progress.update(
            ProgressSnapshot(
                command="discover",
                step="scan_reads",
                total_reads=discover_progress_total_reads(config, totals),
                total_bases=discover_progress_total_bases(config, totals),
                extra=discover_progress_extra(0, totals),
            )
        )

    candidate_tmp_path.write_text(candidate_reads_header() + "\n", encoding="utf-8")
    candidate_tmp_path.replace(candidate_path)
    with candidate_path.open("a", encoding="utf-8") as handle:
        chunk: list[ReadScanTask] = []
        chunk_bases = 0
        selected_reads = 0
        selected_bases = 0

        def handle_scan_result(result: ReadScanResult) -> None:
            nonlocal processed_reads
            nonlocal processed_bases
            nonlocal skipped_short_reads
            nonlocal skipped_short_kmer
            nonlocal skipped_low_complexity
            nonlocal seed_overflow_count

            processed_reads += 1
            processed_bases += result.read_bases
            skipped_short_reads += result.skipped_short_reads
            skipped_short_kmer += result.skipped_short_kmer
            skipped_low_complexity += result.skipped_low_complexity
            seed_overflow_count += result.seed_overflow_count
            if result.candidate is not None:
                candidate = renumber_candidate(
                    result.candidate,
                    f"TXC{len(candidates) + 1:06d}",
                )
                candidates.append(candidate)
                handle.write(format_candidate_read(candidate) + "\n")

            if processed_reads % config.progress_every == 0:
                poll_discover_count_future(count_future, totals, logger, config)
                log_discover_progress(
                    logger,
                    processed_reads,
                    processed_bases,
                    len(candidates),
                    started,
                    discover_progress_total_reads(config, totals),
                    discover_progress_total_bases(config, totals),
                )
                update_discover_terminal_progress(
                    progress,
                    "scan_reads",
                    processed_reads,
                    processed_bases,
                    len(candidates),
                    config,
                    totals,
                )

        def flush_chunk(executor: ThreadPoolExecutor | None) -> None:
            nonlocal chunk_bases
            if not chunk:
                return
            try:
                if executor is None:
                    results = scan_discover_chunk(chunk, config)
                else:
                    batches = split_read_scan_tasks(chunk, effective_threads)
                    results = (
                        result
                        for batch_results in executor.map(
                            scan_discover_chunk,
                            batches,
                            [config] * len(batches),
                        )
                        for result in batch_results
                    )
                for result in results:
                    handle_scan_result(result)
                handle.flush()
            except RustBackendUnavailable as exc:
                logger.error("rust_backend_unavailable=%s", exc)
                raise ValueError(str(exc)) from exc
            finally:
                chunk.clear()
                chunk_bases = 0

        executor_context = (
            ThreadPoolExecutor(max_workers=effective_threads)
            if parallel_scan
            else nullcontext(None)
        )
        with executor_context as executor:
            record_iter = (
                read_fasta_many_round_robin(sequence_paths)
                if auto_discovery_uses_round_robin(sequence_paths, effective_max_read_bases, config)
                else read_fasta_many(sequence_paths)
            )
            for record in record_iter:
                poll_discover_count_future(count_future, totals, logger, config)
                effective_max_reads, effective_max_read_bases = maybe_enable_auto_discovery_budget(
                    config,
                    totals,
                    logger,
                    effective_max_reads,
                    effective_max_read_bases,
                )
                if effective_max_reads is not None and selected_reads >= effective_max_reads:
                    break
                if config.sample_rate < 1.0 and rng.random() > config.sample_rate:
                    continue
                if (
                    effective_max_read_bases is not None
                    and selected_bases + len(record.sequence) > effective_max_read_bases
                ):
                    logger.info(
                        "limit_reached=max_read_bases configured_bases=%s next_read_bases=%s",
                        effective_max_read_bases,
                        len(record.sequence),
                    )
                    break

                selected_reads += 1
                selected_bases += len(record.sequence)
                chunk.append(ReadScanTask(record=record))
                chunk_bases += len(record.sequence)
                if len(chunk) >= config.chunk_size or chunk_bases >= config.chunk_bases:
                    flush_chunk(executor)
            flush_chunk(executor)

        poll_discover_count_future(count_future, totals, logger, config)
        log_discover_progress(
            logger,
            processed_reads,
            processed_bases,
            len(candidates),
            started,
            discover_progress_total_reads(config, totals),
            discover_progress_total_bases(config, totals),
        )
        update_discover_terminal_progress(
            progress,
            "scan_reads",
            processed_reads,
            processed_bases,
            len(candidates),
            config,
            totals,
        )
        logger.info(
            "filter_summary skipped_short_reads=%s skipped_short_kmer=%s "
            "skipped_low_complexity=%s seed_overflow_count=%s",
            skipped_short_reads,
            skipped_short_kmer,
            skipped_low_complexity,
            seed_overflow_count,
        )

    if not candidates:
        if processed_reads and skipped_short_kmer == processed_reads:
            raise ValueError(
                f"No reads were long enough for --kmer-size {config.kmer_size}; "
                f"processed {processed_reads} reads"
            )
        raise ValueError(
            "No tandem repeat candidates found after k-mer spacing prefilter; "
            f"processed {processed_reads} reads"
        )

    update_discover_terminal_progress(
        progress,
        "cluster_families",
        processed_reads,
        processed_bases,
        len(candidates),
        config,
        totals,
    )
    families = cluster_candidates(candidates, config.min_support_reads)
    if not families:
        raise ValueError(
            "No repeat families passed the minimum support threshold; "
            "try lowering --min-support-reads for toy data or check the reads"
        )
    update_discover_terminal_progress(
        progress,
        "compare_families",
        processed_reads,
        processed_bases,
        len(candidates),
        config,
        totals,
    )
    similarities = compare_families(families, k=config.kmer_size)
    families = annotate_family_redundancy(families, similarities)
    update_discover_terminal_progress(
        progress,
        "write_outputs",
        processed_reads,
        processed_bases,
        len(candidates),
        config,
        totals,
    )
    write_monomers(config.outdir / "monomers.fa", families)
    write_families(config.outdir / "families.tsv", families)
    write_family_similarity(config.outdir / "family_similarity.tsv", similarities)
    if config.collapse_redundant_families:
        collapsed_families, collapse_records = collapse_redundant_families(families, similarities)
        write_monomers(config.outdir / "collapsed_monomers.fa", collapsed_families)
        write_families(config.outdir / "collapsed_families.tsv", collapsed_families)
        write_family_collapse(config.outdir / "family_collapse.tsv", collapse_records)
    return candidates, families


def update_discover_terminal_progress(
    progress: TerminalProgress | None,
    step: str,
    processed_reads: int,
    processed_bases: int,
    candidate_reads: int,
    config: DiscoverConfig,
    totals: DiscoverProgressTotals,
) -> None:
    if progress is None:
        return
    progress.update(
        ProgressSnapshot(
            command="discover",
            step=step,
            processed_reads=processed_reads,
            processed_bases=processed_bases,
            total_reads=discover_progress_total_reads(config, totals),
            total_bases=discover_progress_total_bases(config, totals),
            extra=discover_progress_extra(candidate_reads, totals),
        )
    )


def poll_discover_count_future(
    count_future: Future[SequenceStats] | None,
    totals: DiscoverProgressTotals,
    logger: logging.Logger,
    config: DiscoverConfig,
) -> None:
    if totals.count_done or count_future is None or not count_future.done():
        return
    stats = count_future.result()
    totals.total_reads = stats.record_count
    totals.total_read_bases = stats.total_bases
    totals.count_done = True
    read_files = len(config.reads) if not isinstance(config.reads, Path) else 1
    logger.info(
        "input_summary read_files=%s total_reads=%s total_bases=%s max_read_length=%s",
        read_files,
        stats.record_count,
        stats.total_bases,
        stats.max_read_length,
    )


def maybe_enable_auto_discovery_budget(
    config: DiscoverConfig,
    totals: DiscoverProgressTotals,
    logger: logging.Logger,
    effective_max_reads: int | None,
    effective_max_read_bases: int | None,
) -> tuple[int | None, int | None]:
    if not config.enable_auto_discovery_budget:
        totals.effective_max_reads = effective_max_reads
        totals.effective_max_read_bases = effective_max_read_bases
        return effective_max_reads, effective_max_read_bases

    if (
        effective_max_read_bases is not None
        and effective_max_reads is None
        and config.max_reads is None
        and config.max_read_bases is None
        and totals.total_reads is not None
        and totals.total_read_bases is not None
    ):
        effective_max_reads = max(
            1,
            round(totals.total_reads * effective_max_read_bases / totals.total_read_bases),
        )

    if effective_max_reads is not None or effective_max_read_bases is not None:
        totals.effective_max_reads = effective_max_reads
        totals.effective_max_read_bases = effective_max_read_bases
        return effective_max_reads, effective_max_read_bases

    budget_bases, reason = infer_auto_discovery_budget(config, totals)
    if budget_bases is None:
        totals.effective_max_reads = effective_max_reads
        totals.effective_max_read_bases = effective_max_read_bases
        return effective_max_reads, effective_max_read_bases

    effective_max_read_bases = budget_bases
    if totals.total_reads and totals.total_read_bases:
        estimated_reads = max(1, round(totals.total_reads * budget_bases / totals.total_read_bases))
        effective_max_reads = estimated_reads
    totals.effective_max_reads = effective_max_reads
    totals.effective_max_read_bases = effective_max_read_bases
    logger.info(
        "auto_discovery_budget enabled reason=%s effective_max_read_bases=%s effective_max_reads=%s",
        reason,
        effective_max_read_bases,
        effective_max_reads,
    )
    return effective_max_reads, effective_max_read_bases


def infer_auto_discovery_budget(
    config: DiscoverConfig,
    totals: DiscoverProgressTotals,
) -> tuple[int | None, str | None]:
    if not config.enable_auto_discovery_budget:
        return None, None

    if config.genome_size is not None:
        genome_budget = max(1, round(config.genome_size * config.target_discovery_coverage))
        budget = min(config.auto_discovery_max_bases, genome_budget)
        if totals.total_read_bases is not None:
            budget = min(budget, totals.total_read_bases)
        return budget, "genome_size_x_target_coverage"

    if totals.total_read_bases is None or totals.total_read_bases < config.auto_discovery_trigger_bases:
        return None, None
    return min(config.auto_discovery_max_bases, totals.total_read_bases), "large_input_default_cap"


def auto_discovery_uses_round_robin(
    sequence_paths: Sequence[Path],
    effective_max_read_bases: int | None,
    config: DiscoverConfig,
) -> bool:
    return (
        effective_max_read_bases is not None
        and len(sequence_paths) > 1
        and config.sample_rate == 1.0
    )


def discover_progress_total_reads(
    config: DiscoverConfig,
    totals: DiscoverProgressTotals,
) -> int | None:
    if totals.effective_max_reads is not None:
        return totals.effective_max_reads
    if config.max_reads is not None:
        return config.max_reads
    if totals.total_reads is None:
        return None
    if config.sample_rate < 1.0:
        return max(1, round(totals.total_reads * config.sample_rate))
    return totals.total_reads


def discover_progress_total_bases(
    config: DiscoverConfig,
    totals: DiscoverProgressTotals,
) -> int | None:
    if totals.effective_max_read_bases is not None:
        return totals.effective_max_read_bases
    if config.max_read_bases is not None:
        return config.max_read_bases
    if config.sample_rate < 1.0:
        return None
    return totals.total_read_bases


def discover_progress_extra(candidate_reads: int, totals: DiscoverProgressTotals) -> str:
    parts = [f"candidates={candidate_reads:,}"]
    if not totals.count_done:
        parts.append("counting_inputs")
    return " ".join(parts)


def scan_discover_read(task: ReadScanTask, config: DiscoverConfig) -> ReadScanResult:
    record = task.record
    sequence = record.sequence
    if len(sequence) < config.min_read_length:
        return ReadScanResult(read_bases=len(sequence), skipped_short_reads=1)

    if config.kmer_backend == "rust":
        candidate, overflow_count = find_best_periodic_candidate_rust(
            record,
            min_period=config.min_monomer_len,
            max_period=config.max_monomer_len,
            min_repeat_span=config.min_repeat_span,
            candidate_index=0,
            kmer_size=config.kmer_size,
            top_periods=config.top_periods,
            min_seed_occurrences=config.min_seed_occurrences,
            min_spacing_support=config.min_spacing_support,
            max_pairs_per_kmer=config.max_pairs_per_kmer,
        )
        return ReadScanResult(
            read_bases=len(sequence),
            candidate=candidate,
            seed_overflow_count=overflow_count,
        )

    short_candidate, short_overflow_count = find_best_short_periodic_candidate_with_stats(
        record,
        min_period=config.min_monomer_len,
        max_period=min(config.max_monomer_len, SHORT_PERIOD_SCAN_MAX),
        min_repeat_span=config.min_repeat_span,
        candidate_index=0,
    )
    if short_candidate is not None:
        return ReadScanResult(
            read_bases=len(sequence),
            candidate=short_candidate,
            seed_overflow_count=short_overflow_count,
        )
    if len(sequence) < config.kmer_size:
        return ReadScanResult(read_bases=len(sequence), skipped_short_kmer=1)
    long_min_period = max(config.min_monomer_len, SHORT_PERIOD_SCAN_MAX + 1)
    if long_min_period > config.max_monomer_len:
        return ReadScanResult(read_bases=len(sequence))
    candidate, overflow_count = find_best_periodic_candidate_with_stats(
        record,
        min_period=long_min_period,
        max_period=config.max_monomer_len,
        min_repeat_span=config.min_repeat_span,
        candidate_index=0,
        kmer_size=config.kmer_size,
        top_periods=config.top_periods,
        min_seed_occurrences=config.min_seed_occurrences,
        min_spacing_support=config.min_spacing_support,
        max_pairs_per_kmer=config.max_pairs_per_kmer,
    )
    return ReadScanResult(
        read_bases=len(sequence),
        candidate=candidate,
        seed_overflow_count=overflow_count,
    )


def scan_discover_chunk(
    tasks: Sequence[ReadScanTask],
    config: DiscoverConfig,
) -> list[ReadScanResult]:
    if config.kmer_backend == "rust":
        return scan_discover_chunk_rust(tasks, config)
    return [scan_discover_read(task, config) for task in tasks]


def scan_discover_chunk_rust(
    tasks: Sequence[ReadScanTask],
    config: DiscoverConfig,
) -> list[ReadScanResult]:
    return scan_discover_chunk_rust_batch(tasks, config)


def scan_discover_chunk_rust_batch(
    tasks: Sequence[ReadScanTask],
    config: DiscoverConfig,
) -> list[ReadScanResult]:
    results: list[ReadScanResult | None] = [None] * len(tasks)
    rust_groups: dict[int, list[tuple[int, FastaRecord, str]]] = {}
    min_batch_length = max(config.min_repeat_span, config.min_monomer_len * 2)

    for index, task in enumerate(tasks):
        record = task.record
        sequence = record.sequence
        if len(sequence) < config.min_read_length:
            results[index] = ReadScanResult(read_bases=len(sequence), skipped_short_reads=1)
            continue

        if len(sequence) < config.kmer_size and config.min_monomer_len > SHORT_PERIOD_SCAN_MAX:
            results[index] = ReadScanResult(read_bases=len(sequence), skipped_short_kmer=1)
            continue
        if len(sequence) < min_batch_length:
            results[index] = ReadScanResult(read_bases=len(sequence))
            continue
        bounded_max_period = min(config.max_monomer_len, len(sequence) // 2)
        if bounded_max_period < config.min_monomer_len:
            results[index] = ReadScanResult(read_bases=len(sequence))
            continue
        rust_groups.setdefault(bounded_max_period, []).append((index, record, sequence))

    for bounded_max_period, group in rust_groups.items():
        rust_results = scan_reads_for_periods(
            [sequence for _, _, sequence in group],
            k=config.kmer_size,
            min_period=config.min_monomer_len,
            max_period=bounded_max_period,
            top_periods=config.top_periods,
            min_seed_occurrences=config.min_seed_occurrences,
            min_spacing_support=config.min_spacing_support,
            max_pairs_per_kmer=config.max_pairs_per_kmer,
            min_repeat_span=config.min_repeat_span,
        )
        for (index, record, sequence), rust_result in zip(group, rust_results):
            candidate = build_candidate_from_period(
                record,
                sequence,
                rust_result.best_period,
                rust_result.periodicity_score,
                candidate_index=0,
                read_start=rust_result.repeat_start,
                read_end=rust_result.repeat_end,
            )
            results[index] = ReadScanResult(
                read_bases=len(sequence),
                candidate=candidate,
                seed_overflow_count=rust_result.overflow_count,
            )

    return [result for result in results if result is not None]


def split_read_scan_tasks(
    tasks: Sequence[ReadScanTask],
    worker_count: int,
) -> list[tuple[ReadScanTask, ...]]:
    if not tasks:
        return []
    batch_count = max(1, min(worker_count, len(tasks)))
    batch_size = (len(tasks) + batch_count - 1) // batch_count
    return [
        tuple(tasks[index : index + batch_size])
        for index in range(0, len(tasks), batch_size)
    ]


def validate_discover_config(config: DiscoverConfig) -> None:
    if config.min_monomer_len <= 0:
        raise ValueError("--min-period must be positive")
    if config.max_monomer_len < config.min_monomer_len:
        raise ValueError("--max-period must be greater than or equal to --min-period")
    if config.min_read_length <= 0:
        raise ValueError("--min-read-length must be positive")
    if config.kmer_size <= 0:
        raise ValueError("--kmer-size must be positive")
    if config.top_periods <= 0:
        raise ValueError("--top-periods must be positive")
    if config.min_seed_occurrences < 2:
        raise ValueError("--min-seed-occurrences must be at least 2")
    if config.min_spacing_support <= 0:
        raise ValueError("--min-spacing-support must be positive")
    if config.max_pairs_per_kmer <= 0:
        raise ValueError("--max-pairs-per-kmer must be positive")
    if config.max_reads is not None and config.max_reads <= 0:
        raise ValueError("--max-reads must be positive when provided")
    if config.max_read_bases is not None and config.max_read_bases <= 0:
        raise ValueError("--max-read-bases must be positive when provided")
    if config.genome_size is not None and config.genome_size <= 0:
        raise ValueError("--genome-size must be positive when provided")
    if not 0 < config.sample_rate <= 1:
        raise ValueError("--sample-rate must be greater than 0 and at most 1")
    if config.target_discovery_coverage <= 0:
        raise ValueError("--target-discovery-coverage must be positive")
    if config.auto_discovery_trigger_bases <= 0:
        raise ValueError("--auto-discovery-trigger-bases must be positive")
    if config.auto_discovery_max_bases <= 0:
        raise ValueError("--auto-discovery-max-bases must be positive")
    if config.progress_every <= 0:
        raise ValueError("--progress-every must be positive")
    if config.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive")
    if config.chunk_bases <= 0:
        raise ValueError("--chunk-bases must be positive")
    if config.threads <= 0:
        raise ValueError("--threads must be positive")
    if config.kmer_backend not in {"python", "rust"}:
        raise ValueError("--kmer-backend must be python or rust")
    if config.kmer_backend == "rust" and config.kmer_size > 31:
        raise ValueError("Rust backend requires --kmer-size at most 31")


def read_fasta(path: Path) -> Iterable[FastaRecord]:
    try:
        for record in read_sequence_records(path):
            yield FastaRecord(
                read_id=record.id,
                description=record.description,
                sequence=record.sequence,
            )
    except SequenceFormatError:
        raise


def read_fasta_many(
    paths: Path | Sequence[Path],
    *,
    check_duplicate_ids_across_files: bool = True,
) -> Iterable[FastaRecord]:
    try:
        for record in read_sequence_records_many(
            paths,
            check_duplicate_ids_across_files=check_duplicate_ids_across_files,
        ):
            yield FastaRecord(
                read_id=record.id,
                description=record.description,
                sequence=record.sequence,
            )
    except SequenceFormatError:
        raise


def read_fasta_many_round_robin(paths: Path | Sequence[Path]) -> Iterable[FastaRecord]:
    sequence_paths = normalize_sequence_paths(paths)
    with DuplicateIdTracker() as seen_ids:
        active = deque(
            (path, iter(read_fasta(path)))
            for path in sequence_paths
        )
        while active:
            path, iterator = active.popleft()
            try:
                record = next(iterator)
            except StopIteration:
                continue
            if not seen_ids.add(record.read_id):
                raise SequenceFormatError(
                    f"Duplicate sequence id across input read files: {record.read_id}"
                )
            yield record
            active.append((path, iterator))


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
    kmer_size: int = 11,
    top_periods: int = 5,
    min_seed_occurrences: int = 2,
    min_spacing_support: int = 2,
    max_pairs_per_kmer: int = 100,
) -> CandidateRepeat | None:
    candidate, _ = find_best_periodic_candidate_with_stats(
        record,
        min_period=min_period,
        max_period=max_period,
        min_repeat_span=min_repeat_span,
        candidate_index=candidate_index,
        kmer_size=kmer_size,
        top_periods=top_periods,
        min_seed_occurrences=min_seed_occurrences,
        min_spacing_support=min_spacing_support,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    return candidate


def find_best_short_periodic_candidate_with_stats(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
) -> tuple[CandidateRepeat | None, int]:
    sequence = record.sequence
    bounded_max_period = min(max_period, len(sequence) // 2)
    if min_period > bounded_max_period or len(sequence) < min_repeat_span:
        return None, 0

    best_period = 0
    best_score = 0.0
    best_start = 0
    best_end = 0
    for period in range(min_period, bounded_max_period + 1):
        score, start, end = best_local_periodicity_score(
            sequence,
            period,
            min_repeat_span,
            acceptance_score=SHORT_PERIOD_ACCEPTANCE_SCORE,
        )
        if (score, end - start, -period) > (best_score, best_end - best_start, -best_period):
            best_period = period
            best_score = score
            best_start = start
            best_end = end

    if best_score < SHORT_PERIOD_ACCEPTANCE_SCORE:
        return None, 0
    return build_candidate_from_period(
        record,
        sequence,
        best_period,
        best_score,
        candidate_index,
        read_start=best_start,
        read_end=best_end,
    ), 0


def find_best_periodic_candidate_with_stats(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
    kmer_size: int,
    top_periods: int,
    min_seed_occurrences: int,
    min_spacing_support: int,
    max_pairs_per_kmer: int,
) -> tuple[CandidateRepeat | None, int]:
    sequence = record.sequence
    if len(sequence) < max(min_repeat_span, min_period * 2):
        return None, 0
    bounded_max_period = min(max_period, len(sequence) // 2)
    if bounded_max_period < min_period or len(sequence) < kmer_size:
        return None, 0

    repeated_positions, overflow_count = extract_repeated_kmer_positions(
        sequence,
        kmer_size,
        min_seed_occurrences=min_seed_occurrences,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    histogram = build_spacing_histogram(
        repeated_positions,
        min_period=min_period,
        max_period=bounded_max_period,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    candidate_periods = select_candidate_periods(
        histogram,
        min_period=min_period,
        max_period=bounded_max_period,
        top_periods=top_periods,
        min_spacing_support=min_spacing_support,
    )
    if not candidate_periods:
        return None, overflow_count

    best_period, best_score, repeat_start, repeat_end = refine_candidate_period_with_interval(
        sequence,
        repeated_positions,
        candidate_periods,
        min_period=min_period,
        max_period=bounded_max_period,
        min_repeat_span=min_repeat_span,
    )

    return build_candidate_from_period(
        record,
        sequence,
        best_period,
        best_score,
        candidate_index,
        read_start=repeat_start,
        read_end=repeat_end,
    ), overflow_count


def find_best_periodic_candidate_rust(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
    kmer_size: int,
    top_periods: int,
    min_seed_occurrences: int,
    min_spacing_support: int,
    max_pairs_per_kmer: int,
) -> tuple[CandidateRepeat | None, int]:
    sequence = record.sequence
    if len(sequence) < max(min_repeat_span, min_period * 2):
        return None, 0
    bounded_max_period = min(max_period, len(sequence) // 2)
    if bounded_max_period < min_period or len(sequence) < kmer_size:
        return None, 0
    result = scan_read_for_periods(
        sequence,
        k=kmer_size,
        min_period=min_period,
        max_period=bounded_max_period,
        top_periods=top_periods,
        min_seed_occurrences=min_seed_occurrences,
        min_spacing_support=min_spacing_support,
        max_pairs_per_kmer=max_pairs_per_kmer,
        min_repeat_span=min_repeat_span,
    )
    return build_candidate_from_period(
        record,
        sequence,
        result.best_period,
        result.periodicity_score,
        candidate_index,
        read_start=result.repeat_start,
        read_end=result.repeat_end,
    ), result.overflow_count


def build_candidate_from_period(
    record: FastaRecord,
    sequence: str,
    best_period: int,
    best_score: float,
    candidate_index: int,
    *,
    read_start: int = 0,
    read_end: int | None = None,
) -> CandidateRepeat | None:
    if best_period <= 0 or best_score < 0.75:
        return None
    if read_end is None:
        read_end = len(sequence)
    read_start = max(0, read_start)
    read_end = min(len(sequence), read_end)
    repeat_span = read_end - read_start
    if repeat_span < best_period * 2:
        return None
    monomer = consensus_monomer(sequence[read_start:read_end], best_period)
    if len(monomer) != best_period:
        return None
    warnings = []
    low_complexity = is_low_complexity(monomer)
    if low_complexity:
        warnings.append("low_complexity_candidate")
    if best_period < 20:
        warnings.append("short_period_candidate")
    confidence = "high" if best_score >= 0.9 and not low_complexity else "medium"
    return CandidateRepeat(
        read_id=record.read_id,
        candidate_id=f"TXC{candidate_index:06d}",
        sequence=monomer,
        read_start=read_start,
        read_end=read_end,
        strand=parse_strand(record.description),
        period_bp=best_period,
        repeat_span_bp=repeat_span,
        unit_count=repeat_span / best_period,
        score=best_score,
        low_complexity_flag=low_complexity,
        confidence=confidence,
        warning=";".join(warnings),
    )


def consensus_monomer(repeat_sequence: str, period: int) -> str:
    """Build a deterministic majority consensus from complete units."""
    unit_count = len(repeat_sequence) // period
    if period <= 0 or unit_count < 2:
        return ""
    units = [
        repeat_sequence[index * period : (index + 1) * period]
        for index in range(unit_count)
    ]
    consensus = []
    for offset in range(period):
        counts = Counter(unit[offset] for unit in units if unit[offset] in "ACGT")
        consensus.append(
            max("ACGT", key=lambda base: (counts.get(base, 0), -"ACGT".index(base)))
            if counts
            else "N"
        )
    return orient_monomer("".join(consensus))


def renumber_candidate(candidate: CandidateRepeat, candidate_id: str) -> CandidateRepeat:
    return CandidateRepeat(
        read_id=candidate.read_id,
        candidate_id=candidate_id,
        sequence=candidate.sequence,
        read_start=candidate.read_start,
        read_end=candidate.read_end,
        strand=candidate.strand,
        period_bp=candidate.period_bp,
        repeat_span_bp=candidate.repeat_span_bp,
        unit_count=candidate.unit_count,
        score=candidate.score,
        low_complexity_flag=candidate.low_complexity_flag,
        confidence=candidate.confidence,
        warning=candidate.warning,
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
    if len(sequence) <= 2:
        return True
    counts = {base: sequence.count(base) for base in "ACGT"}
    if max(counts.values()) / len(sequence) >= 0.8:
        return True
    if len(sequence) >= 4 and len(sequence) % 2 == 0:
        first_dinucleotide = sequence[:2]
        if all(sequence[index : index + 2] == first_dinucleotide for index in range(2, len(sequence), 2)):
            return True
    return False


def cluster_candidates(
    candidates: Sequence[CandidateRepeat],
    min_support_reads: int,
    period_tolerance_bp: int = 5,
) -> list[RepeatFamily]:
    clusters: list[CandidateCluster] = []
    for candidate in sorted(candidates, key=lambda item: item.period_bp):
        candidate_kmers = candidate_kmer_sketch(candidate.sequence)
        placed = False
        for cluster in clusters:
            center = cluster.center
            tolerance = 0 if min(candidate.period_bp, center) <= SHORT_PERIOD_SCAN_MAX else period_tolerance_bp
            if (
                abs(candidate.period_bp - center) <= tolerance
                and candidate_sequences_compatible(
                    candidate.sequence,
                    candidate_kmers,
                    cluster.representative.sequence,
                    cluster.representative_kmers,
                )
            ):
                cluster.add(candidate, candidate_kmers)
                placed = True
                break
        if not placed:
            clusters.append(CandidateCluster.from_candidate(candidate, candidate_kmers))

    supported = [cluster for cluster in clusters if cluster.support_read_count >= min_support_reads]
    supported.sort(key=lambda cluster: (-cluster.support_read_count, cluster.center))

    families = []
    for index, cluster in enumerate(supported, start=1):
        dominant_period = Counter(member.period_bp for member in cluster.members).most_common(1)[0][0]
        period_members = [member for member in cluster.members if member.period_bp == dominant_period]
        monomer_sequence = cluster_consensus(period_members)
        mean_identity = cluster.mean_identity
        low_complexity = is_low_complexity(monomer_sequence)
        warning = "low_complexity_family" if low_complexity else ""
        confidence = "high" if cluster.support_read_count >= max(3, min_support_reads) and mean_identity >= 0.9 else "medium"
        families.append(
            RepeatFamily(
                family_id=f"TXF{index:06d}",
                monomer_id=f"TXM{index:06d}",
                monomer_sequence=monomer_sequence,
                monomer_length_bp=dominant_period,
                support_read_count=cluster.support_read_count,
                support_span_bp=cluster.support_span_sum,
                mean_identity=mean_identity,
                low_complexity_flag=low_complexity,
                confidence=confidence,
                warning=warning,
            )
        )
    return families


def orient_monomer(sequence: str) -> str:
    """Canonicalize strand and cyclic phase for tandem-repeat monomers."""
    normalized = sequence.upper()
    if not normalized:
        return normalized
    reverse = reverse_complement(normalized)
    return min(minimal_rotation(normalized), minimal_rotation(reverse))


def minimal_rotation(sequence: str) -> str:
    """Return the lexicographically minimal rotation using Booth's algorithm."""
    if not sequence:
        return sequence
    doubled = sequence + sequence
    length = len(sequence)
    left = 0
    right = 1
    offset = 0
    while left < length and right < length and offset < length:
        a = doubled[left + offset]
        b = doubled[right + offset]
        if a == b:
            offset += 1
            continue
        if a > b:
            left = left + offset + 1
            if left <= right:
                left = right + 1
        else:
            right = right + offset + 1
            if right <= left:
                right = left + 1
        offset = 0
    start = min(left, right)
    return doubled[start : start + length]


def candidate_kmer_sketch(sequence: str, k: int = 11) -> frozenset[int]:
    if not sequence:
        return frozenset()
    effective_k = min(k, len(sequence))
    circular = sequence + sequence[: effective_k - 1]
    return frozenset(
        code
        for position, code in iter_canonical_kmer_codes(circular, effective_k)
        if position < len(sequence)
    )


def candidate_sequences_compatible(
    sequence_a: str,
    kmers_a: frozenset[int],
    sequence_b: str,
    kmers_b: frozenset[int],
) -> bool:
    if not kmers_a or not kmers_b:
        return sequence_a == sequence_b
    shared_fraction = len(kmers_a & kmers_b) / min(len(kmers_a), len(kmers_b))
    threshold = 0.65 if min(len(sequence_a), len(sequence_b)) < 20 else 0.35
    return shared_fraction >= threshold


def cluster_consensus(candidates: Sequence[CandidateRepeat]) -> str:
    if not candidates:
        return ""
    template = max(candidates, key=lambda item: (item.score, item.repeat_span_bp)).sequence
    aligned = [best_cyclic_alignment(template, candidate.sequence) for candidate in candidates]
    consensus = []
    for offset in range(len(template)):
        counts = Counter(sequence[offset] for sequence in aligned if sequence[offset] in "ACGT")
        consensus.append(
            max("ACGT", key=lambda base: (counts.get(base, 0), -"ACGT".index(base)))
            if counts
            else "N"
        )
    return orient_monomer("".join(consensus))


def best_cyclic_alignment(reference: str, sequence: str) -> str:
    if len(reference) != len(sequence) or not reference:
        return sequence[: len(reference)].ljust(len(reference), "N")
    length = len(reference)
    effective_k = min(9, length)
    reference_anchor = reference[:effective_k]
    candidate_shifts: set[tuple[str, int]] = set()
    for oriented in (sequence, reverse_complement(sequence)):
        doubled = oriented + oriented[: effective_k - 1]
        start = doubled.find(reference_anchor)
        while start >= 0 and start < length and len(candidate_shifts) < 128:
            candidate_shifts.add((oriented, start))
            start = doubled.find(reference_anchor, start + 1)
        candidate_shifts.add((oriented, 0))
    if len(candidate_shifts) <= 2 and length <= 256:
        candidate_shifts = {
            (oriented, shift)
            for oriented in (sequence, reverse_complement(sequence))
            for shift in range(length)
        }
    best = reference
    best_matches = -1
    for oriented, shift in candidate_shifts:
        rotated = oriented[shift:] + oriented[:shift]
        matches = sum(left == right for left, right in zip(reference, rotated))
        if matches > best_matches or (matches == best_matches and rotated < best):
            best_matches = matches
            best = rotated
    return best


def compare_families(families: Sequence[RepeatFamily], k: int = 11) -> list[FamilySimilarity]:
    """Compare representative monomers to flag possible family redundancy."""
    similarities: list[FamilySimilarity] = []
    for index, family_a in enumerate(families):
        for family_b in families[index + 1 :]:
            similarities.append(compare_family_pair(family_a, family_b, k=k))
    return similarities


def compare_family_pair(family_a: RepeatFamily, family_b: RepeatFamily, k: int = 11) -> FamilySimilarity:
    kmers_a = canonical_kmer_set(family_a.monomer_sequence, k)
    kmers_b = canonical_kmer_set(family_b.monomer_sequence, k)
    shared = len(kmers_a & kmers_b)
    union = len(kmers_a | kmers_b)
    kmer_jaccard = shared / union if union else 0.0
    shared_fraction = shared / min(len(kmers_a), len(kmers_b)) if kmers_a and kmers_b else 0.0
    identity, overlap, orientation = best_local_identity(
        family_a.monomer_sequence,
        family_b.monomer_sequence,
    )
    shorter = max(1, min(family_a.monomer_length_bp, family_b.monomer_length_bp))
    longer = max(family_a.monomer_length_bp, family_b.monomer_length_bp)
    overlap_fraction_shorter = overlap / shorter
    length_ratio = longer / shorter
    relationship, redundant_candidate, notes = classify_family_relationship(
        kmer_jaccard=kmer_jaccard,
        shared_kmer_fraction=shared_fraction,
        local_identity=identity,
        local_overlap_fraction_shorter=overlap_fraction_shorter,
        length_ratio=length_ratio,
    )
    return FamilySimilarity(
        family_a=family_a.family_id,
        family_b=family_b.family_id,
        length_a_bp=family_a.monomer_length_bp,
        length_b_bp=family_b.monomer_length_bp,
        kmer_jaccard=kmer_jaccard,
        shared_kmer_fraction=shared_fraction,
        local_identity=identity,
        local_overlap_bp=overlap,
        local_overlap_fraction_shorter=overlap_fraction_shorter,
        length_ratio=length_ratio,
        orientation=orientation,
        relationship=relationship,
        redundant_candidate=redundant_candidate,
        notes=notes,
    )


def canonical_kmer_set(sequence: str, k: int) -> set[str]:
    if k <= 0 or len(sequence) < k:
        return set()
    return {canonical_kmer(sequence[index : index + k]) for index in range(len(sequence) - k + 1)}


def best_local_identity(sequence_a: str, sequence_b: str) -> tuple[float, int, str]:
    """Return the best ungapped local identity over both orientations."""
    min_overlap = min(50, len(sequence_a), len(sequence_b))
    best_identity = 0.0
    best_overlap = 0
    best_orientation = "forward"
    for orientation, oriented_b in (
        ("forward", sequence_b),
        ("reverse", reverse_complement(sequence_b)),
    ):
        for offset in range(-len(oriented_b) + 1, len(sequence_a)):
            start_a = max(0, offset)
            start_b = max(0, -offset)
            overlap = min(len(sequence_a) - start_a, len(oriented_b) - start_b)
            if overlap < min_overlap:
                continue
            matches = sum(
                1
                for base_index in range(overlap)
                if sequence_a[start_a + base_index] == oriented_b[start_b + base_index]
            )
            identity = matches / overlap
            if identity > best_identity or (identity == best_identity and overlap > best_overlap):
                best_identity = identity
                best_overlap = overlap
                best_orientation = orientation
    return best_identity, best_overlap, best_orientation


def classify_family_relationship(
    *,
    kmer_jaccard: float,
    shared_kmer_fraction: float,
    local_identity: float,
    local_overlap_fraction_shorter: float,
    length_ratio: float,
) -> tuple[str, bool, str]:
    near_integer_ratio = abs(length_ratio - round(length_ratio)) <= 0.05 and round(length_ratio) >= 2
    if (
        near_integer_ratio
        and local_identity >= 0.75
        and local_overlap_fraction_shorter >= 0.75
        and (kmer_jaccard >= 0.25 or shared_kmer_fraction >= 0.5)
    ):
        return (
            "possible_higher_order_or_partial",
            False,
            "Moderate or high similarity with near-integer length ratio; inspect as possible higher-order unit, dimer, or related family.",
        )
    if (
        local_identity >= 0.9
        and local_overlap_fraction_shorter >= 0.85
        and shared_kmer_fraction >= 0.8
        and length_ratio <= 1.20
    ):
        return (
            "likely_redundant",
            True,
            "High sequence similarity; shorter or lower-support representative may be redundant.",
        )
    if (
        local_identity >= 0.75
        and local_overlap_fraction_shorter >= 0.75
        and (kmer_jaccard >= 0.25 or shared_kmer_fraction >= 0.5)
    ):
        return (
            "possible_higher_order_or_partial",
            False,
            "Moderate local similarity; inspect as possible higher-order unit, partial duplicate, or related family.",
        )
    return "distinct", False, ""


def annotate_family_redundancy(
    families: Sequence[RepeatFamily],
    similarities: Sequence[FamilySimilarity],
) -> list[RepeatFamily]:
    warnings_by_family: dict[str, list[str]] = {family.family_id: [] for family in families}
    for similarity in similarities:
        if similarity.relationship == "distinct":
            continue
        warning = f"{similarity.relationship}:{similarity.family_a}-{similarity.family_b}"
        warnings_by_family[similarity.family_a].append(warning)
        warnings_by_family[similarity.family_b].append(warning)

    annotated: list[RepeatFamily] = []
    for family in families:
        warnings = [family.warning] if family.warning else []
        warnings.extend(warnings_by_family[family.family_id])
        annotated.append(replace(family, warning=";".join(warnings)))
    return annotated


def collapse_redundant_families(
    families: Sequence[RepeatFamily],
    similarities: Sequence[FamilySimilarity],
) -> tuple[list[RepeatFamily], list[FamilyCollapse]]:
    """Collapse only pairs classified as likely redundant."""
    family_by_id = {family.family_id: family for family in families}
    parent = {family.family_id: family.family_id for family in families}

    def find(family_id: str) -> str:
        while parent[family_id] != family_id:
            parent[family_id] = parent[parent[family_id]]
            family_id = parent[family_id]
        return family_id

    def union(left: str, right: str) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for similarity in similarities:
        if similarity.relationship == "likely_redundant":
            union(similarity.family_a, similarity.family_b)

    components: dict[str, list[RepeatFamily]] = {}
    for family in families:
        components.setdefault(find(family.family_id), []).append(family)

    retained_by_family: dict[str, str] = {}
    for component in components.values():
        retained = choose_retained_family(component)
        for family in component:
            retained_by_family[family.family_id] = retained.family_id

    similarity_lookup = {
        frozenset((similarity.family_a, similarity.family_b)): similarity
        for similarity in similarities
        if similarity.relationship == "likely_redundant"
    }
    records: list[FamilyCollapse] = []
    collapsed: list[RepeatFamily] = []
    for family in families:
        retained_id = retained_by_family[family.family_id]
        if family.family_id == retained_id:
            collapsed.append(family)
            records.append(
                FamilyCollapse(
                    original_family_id=family.family_id,
                    retained_family_id=retained_id,
                    action="retained",
                    reason="selected_representative" if len(components[find(family.family_id)]) > 1 else "no_likely_redundant_match",
                    relationship="",
                    similarity_metrics="",
                    notes="Possible higher-order or partial relationships are not collapsed.",
                )
            )
            continue
        similarity = similarity_lookup.get(frozenset((family.family_id, retained_id)))
        records.append(
            FamilyCollapse(
                original_family_id=family.family_id,
                retained_family_id=retained_id,
                action="collapsed",
                reason="likely_redundant_to_retained_family",
                relationship=similarity.relationship if similarity else "likely_redundant",
                similarity_metrics=format_similarity_metrics(similarity) if similarity else "",
                notes="Collapsed only because the relationship was likely_redundant.",
            )
        )
    return collapsed, records


def choose_retained_family(families: Sequence[RepeatFamily]) -> RepeatFamily:
    return max(
        families,
        key=lambda family: (
            family.support_read_count,
            family.support_span_bp,
            family.mean_identity,
            -family.monomer_length_bp,
            family.family_id,
        ),
    )


def format_similarity_metrics(similarity: FamilySimilarity | None) -> str:
    if similarity is None:
        return ""
    return (
        f"kmer_jaccard={similarity.kmer_jaccard:.4f};"
        f"shared_kmer_fraction={similarity.shared_kmer_fraction:.4f};"
        f"local_identity={similarity.local_identity:.4f};"
        f"local_overlap_bp={similarity.local_overlap_bp};"
        f"length_ratio={similarity.length_ratio:.4f}"
    )


def sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def log_discover_progress(
    logger: logging.Logger,
    processed_reads: int,
    processed_bases: int,
    candidate_reads: int,
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
    estimated_remaining = (
        f"{min(remaining_estimates):.1f}" if remaining_estimates else "unknown"
    )
    logger.info(
        "progress processed_reads=%s processed_bases=%s candidate_reads=%s "
        "elapsed_seconds=%.3f reads_per_second=%.3f mb_per_second=%.3f "
        "estimated_remaining_seconds=%s",
        processed_reads,
        processed_bases,
        candidate_reads,
        elapsed,
        reads_per_second,
        mb_per_second,
        estimated_remaining,
    )


def candidate_reads_header() -> str:
    return (
        "read_id\tcandidate_id\tread_start\tread_end\tstrand\tperiod_bp\t"
        "repeat_span_bp\tunit_count\tscore\tlow_complexity_flag\tconfidence\twarning"
    )


def format_candidate_read(candidate: CandidateRepeat) -> str:
    return "\t".join(
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


def write_candidate_reads(path: Path, candidates: Sequence[CandidateRepeat]) -> None:
    lines = [candidate_reads_header()]
    lines.extend(format_candidate_read(candidate) for candidate in candidates)
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


def family_similarity_header() -> str:
    return (
        "family_a\tfamily_b\tlength_a_bp\tlength_b_bp\tkmer_jaccard\t"
        "shared_kmer_fraction\tlocal_identity\tlocal_overlap_bp\t"
        "local_overlap_fraction_shorter\tlength_ratio\torientation\t"
        "relationship\tredundant_candidate\tnotes"
    )


def format_family_similarity(similarity: FamilySimilarity) -> str:
    return "\t".join(
        [
            similarity.family_a,
            similarity.family_b,
            str(similarity.length_a_bp),
            str(similarity.length_b_bp),
            f"{similarity.kmer_jaccard:.4f}",
            f"{similarity.shared_kmer_fraction:.4f}",
            f"{similarity.local_identity:.4f}",
            str(similarity.local_overlap_bp),
            f"{similarity.local_overlap_fraction_shorter:.4f}",
            f"{similarity.length_ratio:.4f}",
            similarity.orientation,
            similarity.relationship,
            str(similarity.redundant_candidate).lower(),
            similarity.notes,
        ]
    )


def write_family_similarity(path: Path, similarities: Sequence[FamilySimilarity]) -> None:
    lines = [family_similarity_header()]
    lines.extend(format_family_similarity(similarity) for similarity in similarities)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def family_collapse_header() -> str:
    return "original_family_id\tretained_family_id\taction\treason\trelationship\tsimilarity_metrics\tnotes"


def format_family_collapse(record: FamilyCollapse) -> str:
    return "\t".join(
        [
            record.original_family_id,
            record.retained_family_id,
            record.action,
            record.reason,
            record.relationship,
            record.similarity_metrics,
            record.notes,
        ]
    )


def write_family_collapse(path: Path, records: Sequence[FamilyCollapse]) -> None:
    lines = [family_collapse_header()]
    lines.extend(format_family_collapse(record) for record in records)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
