from __future__ import annotations

import pytest

from tandemx.utils.threads import (
    DEFAULT_DISCOVER_THREADS,
    HARD_MAX_THREADS,
    discover_thread_limit,
    resolve_count_threads,
    resolve_discover_threads,
)


def test_discover_thread_limit_uses_half_available_with_hard_cap() -> None:
    assert discover_thread_limit(cpu_count=4) == 2
    assert discover_thread_limit(cpu_count=256) == HARD_MAX_THREADS
    assert discover_thread_limit(cpu_count=None) >= 1


def test_resolve_discover_threads_defaults_to_eight_or_host_cap(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("tandemx.utils.threads.os.cpu_count", lambda: 64)
    assert resolve_discover_threads(None) == DEFAULT_DISCOVER_THREADS

    monkeypatch.setattr("tandemx.utils.threads.os.cpu_count", lambda: 8)
    assert resolve_discover_threads(None) == 4


def test_resolve_discover_threads_rejects_explicit_values_above_host_cap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("tandemx.utils.threads.os.cpu_count", lambda: 8)
    with pytest.raises(ValueError, match="--threads must be at most 4"):
        resolve_discover_threads(5)


def test_resolve_count_threads_uses_at_most_four_threads() -> None:
    assert resolve_count_threads(None, discover_threads=8) == 4
    assert resolve_count_threads(None, discover_threads=2) == 2
    assert resolve_count_threads(1, discover_threads=8) == 1
    with pytest.raises(ValueError, match="--count-threads must be at most 4"):
        resolve_count_threads(5, discover_threads=8)
