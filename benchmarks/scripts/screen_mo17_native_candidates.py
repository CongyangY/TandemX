"""Reproduce the frozen six-region Mo17 assembly-only periodicity screen."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PERIOD_RANGES = {
    "Mo17_S9_R5": range(107, 118),    # sat112
    "Mo17_S8_R7": range(151, 162),    # CentC
    "Mo17_S9_R7": range(263, 274),    # sat268
    "Mo17_S7_R9": range(175, 186),    # knob180
    "Mo17_S8_R18": range(151, 162),   # CentC
    "Mo17_S6_R9": range(340, 361),    # TR-1
}


def read_regions(path: Path) -> dict[str, str]:
    regions: dict[str, str] = {}
    key: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split("|", 1)[0]
            if key in regions:
                raise ValueError(f"duplicate candidate {key}")
            regions[key] = ""
        elif key is None:
            raise ValueError("sequence before header")
        else:
            regions[key] += line.upper()
    if set(regions) != set(PERIOD_RANGES) or any(set(x) - set("ACGT") for x in regions.values()):
        raise ValueError("unexpected candidate set or non-ACGT sequence")
    return regions


def shift_identity(sequence: str, period: int) -> float:
    if not 0 < period < len(sequence):
        raise ValueError("period out of range")
    return sum(a == b for a, b in zip(sequence[:-period], sequence[period:])) / (len(sequence) - period)


def screen(regions: dict[str, str]) -> list[dict[str, str | int | float]]:
    rows = []
    for region_id, sequence in regions.items():
        best: tuple[float, int, int] | None = None
        # Each extracted record has exactly 3 kb of assembly margin on both sides.
        # Windows are wholly inside the published region.
        for start in range(3000, len(sequence) - 6000 + 1, 250):
            window = sequence[start:start + 3000]
            for period in PERIOD_RANGES[region_id]:
                score = shift_identity(window, period)
                if best is None or score > best[0]:
                    best = (score, start, period)
        if best is None:
            raise ValueError(f"no 3 kb window in {region_id}")
        rows.append({"published_region_id": region_id,
                     "window_offset_in_extracted_region_bp": best[1],
                     "best_period_bp": best[2],
                     "period_shift_identity": round(best[0], 9)})
    return sorted(rows, key=lambda row: -row["period_shift_identity"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=["published_region_id",
                                                "window_offset_in_extracted_region_bp",
                                                "best_period_bp", "period_shift_identity"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(screen(read_regions(args.regions)))


if __name__ == "__main__":
    main()
