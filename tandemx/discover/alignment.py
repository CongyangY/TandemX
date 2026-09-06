"""Read-local banded self alignment; Python reference for the native kernel.

Coordinates are zero-based, half-open. The band follows an approximate repeat
period; it is not a whole-read all-versus-all alignment. A fixed score drawdown
terminates paths through unrelated sequence. This is a heuristic, not an exact
Smith-Waterman optimum, and does not infer unobserved array sequence.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from statistics import median_low


MAX_TRACE_CELLS = 32_000_000


@dataclass(frozen=True)
class AlignmentHit:
    start: int
    end: int
    period: int
    matches: int
    columns: int
    gaps: int
    alignment_score: int
    pairs: tuple[tuple[int, int], ...]

    @property
    def identity(self) -> float:
        return self.matches / self.columns if self.columns else 0.0


def banded_self_align(
    sequence: str,
    period: int,
    min_span: int,
    *,
    band: int | None = None,
    x_drop: int = 40,
    backend: str = "python",
) -> list[AlignmentHit]:
    """Return nonoverlapping local paths near one forward period offset.

    Scores: match +2, substitution -3, each gap base -4. Both Ns and other
    ambiguous bases break paths. Identity includes gap columns in its denominator.
    The period is the lower median of aligned coordinate offsets, not a claim
    that every individual repeat unit has that length.
    """
    if period <= 0 or min_span <= 0 or x_drop <= 0:
        raise ValueError("period, min_span and x_drop must be positive")
    if band is None:
        band = min(period - 1, max(3, ceil(period * 0.08)))
    if band < 0 or band >= period:
        raise ValueError("band must be nonnegative and smaller than period")
    sequence = sequence.upper()
    if not sequence.isascii():
        raise ValueError("sequence must be ASCII")
    if (len(sequence) + 1) * (2 * band + 1) > MAX_TRACE_CELLS:
        raise ValueError("Self-alignment trace exceeds 32 million cells; split this read or narrow the period band")
    if backend == "rust":
        from tandemx.discover.rust_backend import RustBackendUnavailable

        try:
            from tandemx import _rust_core
        except ImportError as exc:
            raise RustBackendUnavailable("Rebuild the TandemX Rust extension for elastic discovery") from exc
        if not hasattr(_rust_core, "banded_self_align"):
            raise RustBackendUnavailable("Rust extension lacks elastic alignment; reinstall TandemX")
        return [AlignmentHit(a, b, p, m, c, g, s, tuple(pairs))
                for a, b, p, m, c, g, s, pairs in
                _rust_core.banded_self_align(sequence, period, min_span, band, x_drop)]
    if backend != "python":
        raise ValueError("alignment backend must be python or rust")
    return _self_align_python(sequence, period, min_span, band, x_drop)


def _self_align_python(sequence: str, period: int, min_span: int,
                       band: int, x_drop: int) -> list[AlignmentHit]:
    n = len(sequence)
    width = 2 * band + 1
    lowest = period - band
    directions = bytearray((n + 1) * width)
    previous = [0] * width
    previous_peak = [0] * width
    endpoints: list[tuple[int, int, int]] = []
    minimum_score = max(min(40, min_span), ceil(0.7 * max(period, min_span - period)))
    for i in range(1, n - lowest + 1):
        current = [0] * width
        current_peak = [0] * width
        row_best = (0, 0)
        for b in range(width):
            j = i + lowest + b
            if j > n:
                break
            if sequence[i - 1] not in "ACGT" or sequence[j - 1] not in "ACGT":
                continue
            score = previous[b] + (2 if sequence[i - 1] == sequence[j - 1] else -3)
            direction = 1
            peak = previous_peak[b]
            if b + 1 < width and previous[b + 1] - 4 > score:
                score, direction, peak = previous[b + 1] - 4, 2, previous_peak[b + 1]
            if b > 0 and current[b - 1] - 4 > score:
                score, direction, peak = current[b - 1] - 4, 3, current_peak[b - 1]
            if score <= 0 or peak - score > x_drop:
                continue
            current[b] = score
            current_peak[b] = max(peak, score)
            directions[i * width + b] = direction
            if score > row_best[0]:
                row_best = (score, b)
        if row_best[0] >= minimum_score:
            endpoints.append((row_best[0], i, row_best[1]))
        previous, previous_peak = current, current_peak
    endpoints.sort(key=lambda item: (-item[0], item[1], item[2]))
    hits: list[AlignmentHit] = []
    # Mark comparison rows after accepting a path; later endpoints inside that
    # same evidence cannot create nested duplicate hits for this band.
    claimed = bytearray(n + 1)
    for score, endpoint, end_band in endpoints:
        if claimed[endpoint]:
            continue
        i, b = endpoint, end_band
        end = i + lowest + b
        pairs: list[tuple[int, int]] = []
        matches = columns = gaps = 0
        while i > 0 and 0 <= b < width:
            direction = directions[i * width + b]
            if not direction:
                break
            j = i + lowest + b
            columns += 1
            if direction == 1:
                pairs.append((i - 1, j - 1))
                matches += sequence[i - 1] == sequence[j - 1]
                i -= 1
            elif direction == 2:
                gaps += 1
                i -= 1
                b += 1
            else:
                gaps += 1
                b -= 1
        if not pairs or matches / columns < 0.75:
            continue
        pairs.reverse()
        measured_period = median_low([right - left for left, right in pairs])
        if end - i < min_span or endpoint - i < 0.8 * measured_period:
            continue
        if any(max(0, min(end, h.end) - max(i, h.start)) * 2 >= min(end - i, h.end - h.start)
               for h in hits):
            continue
        hits.append(AlignmentHit(i, end, measured_period, matches, columns, gaps, score, tuple(pairs)))
        claimed[i: end + 1] = b"\1" * (end + 1 - i)
    return sorted(hits, key=lambda hit: (hit.start, hit.end, hit.period))


def global_align_ops(reference: str, query: str, *, backend: str = "python") -> str:
    """Banded global alignment operations: M paired, D reference-only, I query-only."""
    n, m = len(reference), len(query)
    band = abs(n - m) + max(4, ceil(max(n, m) * 0.12))
    if (n + 1) * (2 * band + 1) > MAX_TRACE_CELLS:
        raise ValueError("Unit alignment exceeds 32 million trace cells")
    if not reference.isascii() or not query.isascii():
        raise ValueError("sequence must be ASCII")
    if backend == "rust":
        from tandemx.discover.rust_backend import RustBackendUnavailable
        try:
            from tandemx import _rust_core
        except ImportError as exc:
            raise RustBackendUnavailable("Rebuild the Rust extension for elastic consensus") from exc
        if not hasattr(_rust_core, "global_align_ops"):
            raise RustBackendUnavailable("Rust extension lacks elastic consensus; reinstall TandemX")
        return _rust_core.global_align_ops(reference, query, band)
    if backend != "python":
        raise ValueError("alignment backend must be python or rust")
    width = 2 * band + 1
    directions = bytearray((n + 1) * width)
    impossible = -1_000_000_000
    previous = [impossible] * width
    for j in range(min(m, band) + 1):
        previous[j + band] = -4 * j
        directions[j + band] = 3 if j else 0
    for i in range(1, n + 1):
        current = [impossible] * width
        for j in range(max(0, i - band), min(m, i + band) + 1):
            b = j - i + band
            if j == 0:
                current[b] = -4 * i
                directions[i * width + b] = 2
                continue
            score = previous[b] + (2 if reference[i - 1] == query[j - 1] else -3)
            direction = 1
            if b + 1 < width and previous[b + 1] - 4 > score:
                score, direction = previous[b + 1] - 4, 2
            if b > 0 and current[b - 1] - 4 > score:
                score, direction = current[b - 1] - 4, 3
            current[b] = score
            directions[i * width + b] = direction
        previous = current
    i, j = n, m
    ops = []
    while i or j:
        direction = directions[i * width + j - i + band]
        if direction == 1:
            ops.append("M"); i -= 1; j -= 1
        elif direction == 2:
            ops.append("D"); i -= 1
        elif direction == 3:
            ops.append("I"); j -= 1
        else:
            raise ValueError("Unit alignment traceback failed")
    return "".join(reversed(ops))
