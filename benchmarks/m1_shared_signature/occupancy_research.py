"""Research-only streaming competitive repeat-window occupancy.

This deliberately has its own schema. It does not implement the frozen
diagnostic-k-mer estimator, expose a public CLI, or claim calibrated physical
copy number.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Iterable, Sequence

try:  # Development dependency only; the research prototype has a Python fallback.
    import edlib
except ImportError:  # pragma: no cover - exercised in clean-environment smoke separately
    edlib = None

from tandemx.io.sequences import SequenceRecord, normalize_sequence_paths, read_sequence_records_many
from tandemx.quantify.mvp import read_monomer_fasta


@dataclass(frozen=True)
class CompetitiveConfig:
    reads: Path | Sequence[Path]
    monomers: Path
    genome_size: int
    outdir: Path
    haploid_depth: float | None = None
    max_reads: int | None = None
    max_read_bases: int | None = None
    window_bp: int = 80
    min_tail_bp: int = 40
    max_edit_fraction: float = 0.225
    min_margin_edits: int = 1
    max_catalogue_families: int = 8
    max_monomer_bp: int = 1024
    max_read_bp: int = 1_000_000


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def periodic_edit_distance(query: str, unit: str) -> int:
    """Minimum edit distance to any phase of an infinitely repeated unit.

    A finite template of query length plus one period minus one covers every
    starting phase. Free reference ends permit an arbitrary phase and offset;
    only query gaps/substitutions/insertions cost edits. Two DP rows bound RAM.
    """
    if not query or not unit or set(query + unit) - set("ACGT"):
        raise ValueError("Periodic edit distance requires nonempty A/C/G/T sequences")
    needed = len(query) + len(unit) - 1
    reference = (unit * math.ceil(needed / len(unit)))[:needed]
    if edlib is not None:
        return int(edlib.align(query, reference, mode="HW", task="distance")["editDistance"])
    previous = [0] * (len(reference) + 1)
    for query_index, base in enumerate(query, start=1):
        current = [query_index] + [0] * len(reference)
        for column, target in enumerate(reference, start=1):
            current[column] = min(
                previous[column - 1] + (base != target),
                previous[column] + 1,
                current[column - 1] + 1,
            )
        previous = current
    return min(previous)


@lru_cache(maxsize=16)
def near_family_pairs(items: tuple[tuple[str, str], ...]) -> frozenset[tuple[str, str]]:
    """Input-only ambiguity for effectively indistinguishable catalogue units."""
    pairs = set()
    for i, (first_name, first_unit) in enumerate(items):
        for second_name, second_unit in items[i + 1 :]:
            if min(periodic_edit_distance(first_unit, second_unit),
                   periodic_edit_distance(first_unit, reverse_complement(second_unit))) <= 1:
                pairs.add(tuple(sorted((first_name, second_name))))
    return frozenset(pairs)


def load_catalogue(path: Path, config: CompetitiveConfig) -> dict[str, str]:
    catalogue: dict[str, str] = {}
    for record in read_monomer_fasta(path):
        if record.family_id in catalogue:
            raise ValueError(f"Duplicate family id: {record.family_id}")
        if len(catalogue) >= config.max_catalogue_families:
            raise ValueError("Experimental competitive catalogue exceeds family cap")
        sequence = record.sequence.upper()
        if not sequence or len(sequence) > config.max_monomer_bp or set(sequence) - set("ACGT"):
            raise ValueError("Experimental competitive monomers must be bounded A/C/G/T sequences")
        catalogue[record.family_id] = sequence
    if not catalogue:
        raise ValueError("No monomers found for competitive occupancy")
    return catalogue


def classify_window(window: str, catalogue: dict[str, str], config: CompetitiveConfig
                    ) -> tuple[str, str | None]:
    """Return assigned, ambiguous, or unknown plus a family/group identifier."""
    if len(window) < config.min_tail_bp or set(window) - set("ACGT"):
        return "unknown", None
    distances = {
        family: min(periodic_edit_distance(window, unit),
                    periodic_edit_distance(window, reverse_complement(unit)))
        for family, unit in catalogue.items()
    }
    best = min(distances.values())
    if best > math.floor(len(window) * config.max_edit_fraction):
        return "unknown", None
    contenders = sorted(family for family, distance in distances.items()
                        if distance - best < config.min_margin_edits)
    pairs = near_family_pairs(tuple(sorted(catalogue.items())))
    contenders = sorted(set(contenders).union(
        family for pair in pairs for family in pair if any(name in pair for name in contenders)
    ))
    if len(contenders) > 1:
        return "ambiguous", "+".join(contenders)
    return "assigned", contenders[0]


def validate_config(config: CompetitiveConfig) -> None:
    paths = normalize_sequence_paths(config.reads)
    if len(set(paths)) != len(paths):
        raise ValueError("--reads must not contain duplicate paths")
    if config.genome_size <= 0 or (config.haploid_depth is not None and (
        not math.isfinite(config.haploid_depth) or config.haploid_depth <= 0
    )):
        raise ValueError("Genome size and optional haploid depth must be positive")
    if (config.max_reads is not None and config.max_reads <= 0) or (
        config.max_read_bases is not None and config.max_read_bases <= 0
    ):
        raise ValueError("Read limits must be positive")
    if not (1 <= config.min_tail_bp <= config.window_bp <= 512):
        raise ValueError("Require 1 <= min tail <= window <= 512 bp")
    if not 0 <= config.max_edit_fraction < 0.5 or config.min_margin_edits < 1:
        raise ValueError("Invalid edit or margin gate")
    if min(config.max_catalogue_families, config.max_monomer_bp, config.max_read_bp) < 1:
        raise ValueError("Resource caps must be positive")


def summarize_records(records: Iterable[SequenceRecord], catalogue: dict[str, str],
                      config: CompetitiveConfig) -> dict:
    """Consume one read at a time; state size is independent of read count."""
    read_count = total_bp = assigned_bp = ambiguous_bp = unknown_bp = windows = 0
    assigned_by_family = Counter({family: 0 for family in catalogue})
    ambiguity_groups: Counter[str] = Counter()
    for record in records:
        if config.max_reads is not None and read_count >= config.max_reads:
            break
        if config.max_read_bases is not None and total_bp >= config.max_read_bases:
            break
        sequence = record.sequence.upper()
        if not sequence or len(sequence) > config.max_read_bp:
            raise ValueError("Read is empty or exceeds the experimental read-length cap")
        read_count += 1
        total_bp += len(sequence)
        for offset in range(0, len(sequence), config.window_bp):
            window = sequence[offset : offset + config.window_bp]
            windows += 1
            outcome, identity = classify_window(window, catalogue, config)
            if outcome == "assigned":
                assert identity is not None
                assigned_bp += len(window)
                assigned_by_family[identity] += len(window)
            elif outcome == "ambiguous":
                assert identity is not None
                ambiguous_bp += len(window)
                ambiguity_groups[identity] += len(window)
            else:
                unknown_bp += len(window)
    if read_count == 0:
        raise ValueError("No reads found for competitive occupancy")
    assert assigned_bp + ambiguous_bp + unknown_bp == total_bp
    depth = config.haploid_depth if config.haploid_depth is not None else total_bp / config.genome_size
    return dict(read_count=read_count, total_read_bp=total_bp, window_count=windows,
                assigned_read_bp=assigned_bp, ambiguous_read_bp=ambiguous_bp,
                unknown_read_bp=unknown_bp, assigned_by_family=dict(assigned_by_family),
                ambiguity_groups=dict(sorted(ambiguity_groups.items())),
                haploid_depth=depth,
                normalization_method="explicit_haploid_depth" if config.haploid_depth is not None
                else "total_read_bp_divided_by_genome_size")


def _atomic_text(path: Path, content: str) -> None:
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                         prefix=path.name + ".", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_outputs(outdir: Path, catalogue: dict[str, str], summary: dict,
                  config: CompetitiveConfig) -> list[dict]:
    rows = []
    header = ("family_id\tmonomer_length_bp\tassigned_read_bp\t"
              "estimated_genomic_bp\testimated_copy_number\thaploid_depth\t"
              "status\twarning")
    lines = [header]
    for family, unit in catalogue.items():
        assigned = summary["assigned_by_family"][family]
        genomic = assigned / summary["haploid_depth"]
        row = dict(family_id=family, monomer_length_bp=len(unit), assigned_read_bp=assigned,
                   estimated_genomic_bp=genomic, estimated_copy_number=genomic / len(unit),
                   haploid_depth=summary["haploid_depth"],
                   status="assigned" if assigned else "no_unique_assignment",
                   warning="experimental_window_occupancy_not_calibrated_physical_copy_number")
        rows.append(row)
        lines.append("\t".join((family, str(len(unit)), str(assigned), f"{genomic:.4f}",
                                f"{row['estimated_copy_number']:.4f}",
                                f"{summary['haploid_depth']:.8f}", row["status"], row["warning"])))
    outdir.mkdir(parents=True, exist_ok=True)
    table = outdir / "competitive_occupancy.tsv"
    receipt = outdir / "competitive_occupancy_summary.json"
    try:
        _atomic_text(table, "\n".join(lines) + "\n")
        payload = dict(complete=True, estimator="experimental_competitive_occupancy_v1",
                       window_bp=config.window_bp, min_tail_bp=config.min_tail_bp,
                       max_edit_fraction=config.max_edit_fraction,
                       min_margin_edits=config.min_margin_edits,
                       genome_size=config.genome_size, **summary)
        _atomic_text(receipt, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    except BaseException:
        table.unlink(missing_ok=True)
        receipt.unlink(missing_ok=True)
        raise
    return rows


def run_competitive_occupancy(config: CompetitiveConfig) -> list[dict]:
    validate_config(config)
    catalogue = load_catalogue(config.monomers, config)
    sequence_paths = normalize_sequence_paths(config.reads)
    summary = summarize_records(read_sequence_records_many(sequence_paths), catalogue, config)
    return write_outputs(config.outdir, catalogue, summary, config)
