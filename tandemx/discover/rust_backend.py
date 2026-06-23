"""Optional wrapper around TandemX's compiled read-local discovery backend."""

from __future__ import annotations

from dataclasses import dataclass
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
