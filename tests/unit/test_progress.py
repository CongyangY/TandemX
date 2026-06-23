from __future__ import annotations

from io import StringIO

from tandemx.utils.progress import ProgressSnapshot, TerminalProgress


def test_terminal_progress_formats_read_rate_eta_and_step() -> None:
    ticks = iter([0.0, 60.0])
    stream = StringIO()
    progress = TerminalProgress(stream=stream, clock=lambda: next(ticks), dynamic=False)

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
    assert "discover scan_reads" in rendered
    assert "50.0%" in rendered
    assert "50/100 reads" in rendered
    assert "1.0Mb" in rendered
    assert "elapsed 1m 00s" in rendered
    assert "est 2m 00s" in rendered
    assert "rem 1m 00s" in rendered
    assert "50 r/min" in rendered
    assert "cand=7" in rendered


def test_terminal_progress_handles_unbounded_total() -> None:
    ticks = iter([0.0, 30.0])
    stream = StringIO()
    progress = TerminalProgress(stream=stream, clock=lambda: next(ticks), dynamic=False)

    progress.update(
        ProgressSnapshot(
            command="quantify",
            step="scan_reads",
            processed_reads=10,
            processed_bases=500_000,
        )
    )

    rendered = stream.getvalue()
    assert "quantify scan_reads" in rendered
    assert "--" in rendered
    assert "est --" in rendered
    assert "rem --" in rendered
    assert "20 r/min" in rendered


def test_dynamic_progress_rewrites_one_bounded_terminal_line() -> None:
    ticks = iter([0.0, 10.0, 20.0])
    stream = StringIO()
    progress = TerminalProgress(
        stream=stream,
        clock=lambda: next(ticks),
        max_columns=80,
    )

    for processed_reads in (1_000, 2_000):
        progress.update(
            ProgressSnapshot(
                command="discover",
                step="scan_reads",
                processed_reads=processed_reads,
                processed_bases=processed_reads * 14_000,
                total_reads=984_870,
                total_bases=13_900_000_000,
                extra="candidates=14",
            )
        )

    rendered = stream.getvalue()
    assert "\n" not in rendered
    refreshes = [part for part in rendered.split("\r\x1b[2K") if part]
    assert len(refreshes) == 2
    assert all(len(part) <= 79 for part in refreshes)
