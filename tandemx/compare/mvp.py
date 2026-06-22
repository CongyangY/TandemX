"""Toy-scale read-vs-assembly abundance comparison MVP."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Sequence


DEFAULT_COLLAPSE_THRESHOLD = 0.6
DEFAULT_OVEREXPANSION_THRESHOLD = 1.5


class ArrayLike(Protocol):
    start: int
    end: int
    family_id: str


@dataclass(frozen=True)
class CompareConfig:
    copy_number: Path
    arrays: Path
    outdir: Path
    collapse_threshold: float = DEFAULT_COLLAPSE_THRESHOLD
    overexpansion_threshold: float = DEFAULT_OVEREXPANSION_THRESHOLD


@dataclass(frozen=True)
class AssemblyArray:
    chrom: str
    start: int
    end: int
    family_id: str
    score: int
    strand: str
    confidence: str
    warning: str


@dataclass(frozen=True)
class AssemblyReadComparison:
    family_id: str
    read_estimated_bp: float
    assembly_estimated_bp: float
    assembly_read_ratio: float
    status: str
    confidence: str
    warning: str


def compare_toy_abundance(config: CompareConfig) -> list[AssemblyReadComparison]:
    validate_compare_config(config)
    arrays = read_arrays_bed(config.arrays)
    comparisons = compare_assembly_to_reads(
        arrays,
        config.copy_number,
        collapse_threshold=config.collapse_threshold,
        overexpansion_threshold=config.overexpansion_threshold,
    )
    config.outdir.mkdir(parents=True, exist_ok=True)
    write_comparisons(config.outdir / "assembly_vs_read_cn.tsv", comparisons)
    return comparisons


def validate_compare_config(config: CompareConfig) -> None:
    if config.collapse_threshold <= 0:
        raise ValueError("--collapse-threshold must be positive")
    if config.overexpansion_threshold <= 0:
        raise ValueError("--overexpansion-threshold must be positive")
    if config.collapse_threshold >= config.overexpansion_threshold:
        raise ValueError("--collapse-threshold must be less than --overexpansion-threshold")


def read_copy_number(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"copy_number.tsv is empty: {path}")
    header = lines[0].split("\t")
    if "family_id" not in header:
        raise ValueError(f"{path} is missing required field: family_id")
    family_index = header.index("family_id")
    if "estimated_bp" in header:
        bp_index = header.index("estimated_bp")
    elif "estimated_repeat_bp" in header:
        bp_index = header.index("estimated_repeat_bp")
    else:
        raise ValueError(f"{path} is missing required field: estimated_bp")
    values: dict[str, float] = defaultdict(float)
    for line_number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != len(header):
            raise ValueError(f"{path} line {line_number} has {len(parts)} fields, expected {len(header)}")
        family_id = parts[family_index]
        if not family_id:
            raise ValueError(f"{path} line {line_number} has empty family_id")
        try:
            values[family_id] += float(parts[bp_index])
        except ValueError as exc:
            raise ValueError(f"{path} line {line_number} estimated_bp is not numeric: {parts[bp_index]}") from exc
    return dict(values)


def read_arrays_bed(path: Path) -> list[AssemblyArray]:
    arrays: list[AssemblyArray] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"{path} line {line_number} must have at least 4 BED fields including family_id")
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"{path} line {line_number} start/end must be integers") from exc
        if start < 0 or end <= start:
            raise ValueError(f"{path} line {line_number} must use 0-based half-open coordinates with end > start")
        family_id = fields[3]
        if not family_id:
            raise ValueError(f"{path} line {line_number} has empty family_id")
        try:
            score = int(fields[4]) if len(fields) > 4 and fields[4] else 0
        except ValueError as exc:
            raise ValueError(f"{path} line {line_number} score must be an integer") from exc
        strand = fields[5] if len(fields) > 5 and fields[5] else "."
        confidence = fields[6] if len(fields) > 6 and fields[6] else "medium"
        warning = fields[7] if len(fields) > 7 else ""
        arrays.append(
            AssemblyArray(
                chrom=fields[0],
                start=start,
                end=end,
                family_id=family_id,
                score=score,
                strand=strand,
                confidence=confidence,
                warning=warning,
            )
        )
    return arrays


def compare_assembly_to_reads(
    arrays: Sequence[ArrayLike],
    copy_number_path: Path | None,
    *,
    collapse_threshold: float = DEFAULT_COLLAPSE_THRESHOLD,
    overexpansion_threshold: float = DEFAULT_OVEREXPANSION_THRESHOLD,
) -> list[AssemblyReadComparison]:
    read_bp = read_copy_number(copy_number_path)
    assembly_bp: dict[str, float] = defaultdict(float)
    for array in arrays:
        assembly_bp[array.family_id] += array.end - array.start
    family_ids = sorted(set(read_bp) | set(assembly_bp))
    comparisons = []
    for family_id in family_ids:
        read_value = read_bp.get(family_id, 0.0)
        assembly_value = assembly_bp.get(family_id, 0.0)
        ratio = assembly_value / read_value if read_value > 0 else 0.0
        status, confidence, warning = classify_assembly_read_ratio(
            read_value,
            assembly_value,
            collapse_threshold=collapse_threshold,
            overexpansion_threshold=overexpansion_threshold,
        )
        comparisons.append(
            AssemblyReadComparison(
                family_id=family_id,
                read_estimated_bp=read_value,
                assembly_estimated_bp=assembly_value,
                assembly_read_ratio=ratio,
                status=status,
                confidence=confidence,
                warning=warning,
            )
        )
    return comparisons


def classify_assembly_read_ratio(
    read_bp: float,
    assembly_bp: float,
    *,
    collapse_threshold: float = DEFAULT_COLLAPSE_THRESHOLD,
    overexpansion_threshold: float = DEFAULT_OVEREXPANSION_THRESHOLD,
) -> tuple[str, str, str]:
    if read_bp <= 0 and assembly_bp <= 0:
        return "low_confidence", "low", "missing_read_and_assembly_estimates"
    if read_bp <= 0:
        return "assembly_only", "low", "missing_read_estimate"
    if assembly_bp <= 0:
        return "reads_only", "medium", "missing_assembly_array"
    ratio = assembly_bp / read_bp
    if ratio < collapse_threshold:
        return "possible_collapse", "medium", ""
    if ratio > overexpansion_threshold:
        return "possible_overexpansion", "medium", ""
    return "consistent", "medium", ""


def write_comparisons(path: Path, comparisons: Sequence[AssemblyReadComparison]) -> None:
    lines = [
        "family_id\tread_estimated_bp\tassembly_estimated_bp\tassembly_read_ratio\tstatus\tconfidence\twarning"
    ]
    for comparison in comparisons:
        lines.append(
            "\t".join(
                [
                    comparison.family_id,
                    f"{comparison.read_estimated_bp:.4f}",
                    f"{comparison.assembly_estimated_bp:.4f}",
                    f"{comparison.assembly_read_ratio:.4f}",
                    comparison.status,
                    comparison.confidence,
                    comparison.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
