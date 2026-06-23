"""Small terminal progress rendering helpers."""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from typing import Callable, TextIO


@dataclass(frozen=True)
class ProgressSnapshot:
    command: str
    step: str
    processed_reads: int = 0
    processed_bases: int = 0
    total_reads: int | None = None
    total_bases: int | None = None
    extra: str = ""


class TerminalProgress:
    """Render read-processing progress without adding a runtime dependency."""

    def __init__(
        self,
        *,
        stream: TextIO | None = None,
        enabled: bool = True,
        width: int = 24,
        clock: Callable[[], float] | None = None,
    ) -> None:
        self.stream = stream if stream is not None else sys.stderr
        self.enabled = enabled
        self.width = width
        self.clock = clock if clock is not None else time.perf_counter
        self.started = self.clock()
        self.last_rendered_length = 0
        self.is_terminal = bool(getattr(self.stream, "isatty", lambda: False)())

    def update(self, snapshot: ProgressSnapshot) -> None:
        if not self.enabled:
            return
        elapsed = max(self.clock() - self.started, 1e-9)
        line = self._format(snapshot, elapsed)
        if self.is_terminal:
            padding = " " * max(0, self.last_rendered_length - len(line))
            self.stream.write("\r" + line + padding)
            self.last_rendered_length = len(line)
        else:
            self.stream.write(line + "\n")
        self.stream.flush()

    def finish(self, command: str, status: str, *, extra: str = "") -> None:
        if not self.enabled:
            return
        elapsed = max(self.clock() - self.started, 0.0)
        line = f"{command} | {status} | elapsed {format_duration(elapsed)}"
        if extra:
            line += f" | {extra}"
        if self.is_terminal:
            padding = " " * max(0, self.last_rendered_length - len(line))
            self.stream.write("\r" + line + padding + "\n")
        else:
            self.stream.write(line + "\n")
        self.stream.flush()

    def _format(self, snapshot: ProgressSnapshot, elapsed: float) -> str:
        fraction = progress_fraction(snapshot)
        bar = progress_bar(fraction, self.width)
        reads_per_minute = snapshot.processed_reads / elapsed * 60
        mb_per_minute = snapshot.processed_bases / 1_000_000 / elapsed * 60
        eta = estimated_remaining_seconds(snapshot, elapsed)
        total_reads = f"/{snapshot.total_reads:,}" if snapshot.total_reads is not None else ""
        total_bases = f"/{format_bases(snapshot.total_bases)}" if snapshot.total_bases is not None else ""
        base = (
            f"{snapshot.command} | {snapshot.step} | {bar} "
            f"{format_percent(fraction)} | {snapshot.processed_reads:,}{total_reads} reads | "
            f"{format_bases(snapshot.processed_bases)}{total_bases} bases | "
            f"{reads_per_minute:,.1f} reads/min | {mb_per_minute:.2f} MB/min | "
            f"ETA {format_eta(eta)}"
        )
        if snapshot.extra:
            base += f" | {snapshot.extra}"
        return base


def progress_fraction(snapshot: ProgressSnapshot) -> float | None:
    fractions: list[float] = []
    if snapshot.total_reads is not None and snapshot.total_reads > 0:
        fractions.append(snapshot.processed_reads / snapshot.total_reads)
    if snapshot.total_bases is not None and snapshot.total_bases > 0:
        fractions.append(snapshot.processed_bases / snapshot.total_bases)
    if not fractions:
        return None
    return min(1.0, max(0.0, max(fractions)))


def estimated_remaining_seconds(snapshot: ProgressSnapshot, elapsed: float) -> float | None:
    estimates: list[float] = []
    if snapshot.total_reads is not None and snapshot.processed_reads > 0:
        reads_per_second = snapshot.processed_reads / elapsed
        estimates.append(max(0, snapshot.total_reads - snapshot.processed_reads) / reads_per_second)
    if snapshot.total_bases is not None and snapshot.processed_bases > 0:
        bases_per_second = snapshot.processed_bases / elapsed
        estimates.append(max(0, snapshot.total_bases - snapshot.processed_bases) / bases_per_second)
    if not estimates:
        return None
    return min(estimates)


def progress_bar(fraction: float | None, width: int) -> str:
    if fraction is None:
        return "[" + "-" * width + "]"
    filled = round(width * fraction)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def format_percent(fraction: float | None) -> str:
    if fraction is None:
        return "unknown"
    return f"{fraction * 100:5.1f}%"


def format_eta(seconds: float | None) -> str:
    if seconds is None:
        return "unknown"
    return format_duration(seconds)


def format_duration(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    if seconds < 60:
        return f"{seconds}s"
    minutes, remaining_seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}m {remaining_seconds:02d}s"
    hours, remaining_minutes = divmod(minutes, 60)
    return f"{hours}h {remaining_minutes:02d}m"


def format_bases(bases: int | None) -> str:
    if bases is None:
        return "unknown"
    if bases < 1_000:
        return f"{bases} bp"
    if bases < 1_000_000:
        return f"{bases / 1_000:.1f} kb"
    if bases < 1_000_000_000:
        return f"{bases / 1_000_000:.1f} Mb"
    return f"{bases / 1_000_000_000:.1f} Gb"
