"""Optional wrapper around TandemX's compiled read-local discovery backend."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


class RustBackendUnavailable(RuntimeError):
    """Raised when the optional compiled extension cannot be imported."""


RUST_MAX_KMER_SIZE = 31


@dataclass(frozen=True)
class RustScanResult:
    candidate_periods: tuple[int, ...]
    spacing_support: tuple[tuple[int, int], ...]
    best_period: int
    periodicity_score: float
    overflow_count: int
    status: str


@dataclass(frozen=True)
class RustSequenceStats:
    record_count: int
    total_bases: int
    max_read_length: int


def rust_backend_available() -> bool:
    try:
        from tandemx import _rust_core  # noqa: F401
    except ImportError:
        return False
    return True


def scan_read_for_periods(
    sequence: str,
    *,
    k: int,
    min_period: int,
    max_period: int,
    top_periods: int,
    min_seed_occurrences: int,
    min_spacing_support: int,
    max_pairs_per_kmer: int,
) -> RustScanResult:
    try:
        from tandemx import _rust_core
    except ImportError as exc:
        raise RustBackendUnavailable(
            "Rust backend is unavailable. Install from the repository with "
            "`pip install -e .` or run `maturin develop`."
        ) from exc

    result = _rust_core.scan_read_for_periods(
        sequence,
        k,
        min_period,
        max_period,
        top_periods,
        min_seed_occurrences,
        min_spacing_support,
        max_pairs_per_kmer,
    )
    return RustScanResult(
        candidate_periods=tuple(result.candidate_periods),
        spacing_support=tuple(tuple(item) for item in result.spacing_support),
        best_period=result.best_period,
        periodicity_score=result.periodicity_score,
        overflow_count=result.overflow_count,
        status=result.status,
    )


def scan_reads_for_periods(
    sequences: Sequence[str],
    *,
    k: int,
    min_period: int,
    max_period: int,
    top_periods: int,
    min_seed_occurrences: int,
    min_spacing_support: int,
    max_pairs_per_kmer: int,
) -> tuple[RustScanResult, ...]:
    try:
        from tandemx import _rust_core
    except ImportError as exc:
        raise RustBackendUnavailable(
            "Rust backend is unavailable. Install from the repository with "
            "`pip install -e .` or run `maturin develop`."
        ) from exc

    results = _rust_core.scan_reads_for_periods(
        list(sequences),
        k,
        min_period,
        max_period,
        top_periods,
        min_seed_occurrences,
        min_spacing_support,
        max_pairs_per_kmer,
    )
    return tuple(
        RustScanResult(
            candidate_periods=tuple(result.candidate_periods),
            spacing_support=tuple(tuple(item) for item in result.spacing_support),
            best_period=result.best_period,
            periodicity_score=result.periodicity_score,
            overflow_count=result.overflow_count,
            status=result.status,
        )
        for result in results
    )


class RustDiagnosticKmerCounter:
    """Small target-only counter; this is not a global k-mer counting backend."""

    def __init__(self, k: int, targets: set[str]) -> None:
        try:
            from tandemx import _rust_core
        except ImportError as exc:
            raise RustBackendUnavailable(
                "Rust backend is unavailable. Install from the repository with "
                "`pip install -e .` or run `maturin develop`."
            ) from exc
        self._counter = _rust_core.DiagnosticKmerCounter(k, sorted(targets))

    def count_sequence(self, sequence: str) -> None:
        self._counter.count_sequence(sequence)

    def counts(self) -> dict[str, int]:
        return dict(self._counter.counts())


def count_sequence_paths_stats(paths: Sequence[Path], *, threads: int) -> RustSequenceStats:
    try:
        from tandemx import _rust_core
    except ImportError as exc:
        raise RustBackendUnavailable(
            "Rust backend is unavailable. Install from the repository with "
            "`pip install -e .` or run `maturin develop`."
        ) from exc
    path_list = list(paths)
    worker_count = max(1, min(threads, len(path_list)))

    def count_one(path: Path):
        return _rust_core.count_sequence_file_stats(str(path))

    if worker_count == 1:
        stats = [count_one(path) for path in path_list]
    else:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            stats = list(executor.map(count_one, path_list))
    return RustSequenceStats(
        record_count=sum(item.record_count for item in stats),
        total_bases=sum(item.total_bases for item in stats),
        max_read_length=max((item.max_read_length for item in stats), default=0),
    )
