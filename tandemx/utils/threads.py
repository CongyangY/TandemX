"""Thread-count policy helpers for CPU-heavy TandemX steps."""

from __future__ import annotations

import os


DEFAULT_DISCOVER_THREADS = 8
HARD_MAX_THREADS = 64


def discover_thread_limit(cpu_count: int | None = None) -> int:
    """Return the maximum user-facing discover thread count for this host."""
    logical_cpus = cpu_count if cpu_count is not None else os.cpu_count()
    half_available = max(1, (logical_cpus or 2) // 2)
    return min(HARD_MAX_THREADS, half_available)


def resolve_discover_threads(requested: int | None) -> int:
    """Resolve CLI thread input; default to 8 but cap defaults on small hosts."""
    limit = discover_thread_limit()
    if requested is None:
        return min(DEFAULT_DISCOVER_THREADS, limit)
    if requested <= 0:
        raise ValueError("--threads must be positive")
    if requested > limit:
        raise ValueError(
            f"--threads must be at most {limit} on this host "
            f"(minimum of {HARD_MAX_THREADS} and half of available logical CPUs)"
        )
    return requested


def effective_discover_threads(configured: int) -> int:
    """Clamp programmatic discover configuration to the host policy."""
    return max(1, min(configured, discover_thread_limit()))
