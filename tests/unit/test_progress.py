from __future__ import annotations

from io import StringIO

from tandemx.utils.progress import ProgressSnapshot, TerminalProgress


def test_terminal_progress_formats_read_rate_eta_and_step() -> None:
    ticks = iter([0.0, 60.0])
    stream = StringIO()
    progress = TerminalProgress(stream=stream, clock=lambda: next(ticks))

    progress.update(
        ProgressSnapshot(
            command="discover",
            step="scan_reads",
            processed_reads=50,
            processed_bases=1_000_000,
            total_reads=100,
            extra="candidates=7",
        )
    )

    rendered = stream.getvalue()
    assert "discover | scan_reads" in rendered
    assert "50.0%" in rendered
    assert "50/100 reads" in rendered
    assert "1.0 Mb bases" in rendered
    assert "elapsed 1m 00s" in rendered
    assert "total est 2m 00s" in rendered
    assert "remaining 1m 00s" in rendered
    assert "50.0 reads/min" in rendered
    assert "candidates=7" in rendered


def test_terminal_progress_handles_unbounded_total() -> None:
    ticks = iter([0.0, 30.0])
    stream = StringIO()
    progress = TerminalProgress(stream=stream, clock=lambda: next(ticks))

    progress.update(
        ProgressSnapshot(
            command="quantify",
            step="scan_reads",
            processed_reads=10,
            processed_bases=500_000,
        )
    )

    rendered = stream.getvalue()
    assert "quantify | scan_reads" in rendered
    assert "--" in rendered
    assert "total est --" in rendered
    assert "remaining --" in rendered
    assert "20.0 reads/min" in rendered
