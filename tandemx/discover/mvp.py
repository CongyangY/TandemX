"""Toy-scale tandem repeat discovery MVP."""

from __future__ import annotations

import hashlib
import logging
import random
import time
from dataclasses import dataclass, replace
from pathlib import Path
from statistics import mean
from typing import Iterable, Sequence

from tandemx.io.sequences import SequenceFormatError, read_sequence_records
from tandemx.discover.spacing import (
    bounded_periodicity_score,
    build_spacing_histogram,
    canonical_kmer,
    extract_repeated_kmer_positions,
    is_low_complexity_kmer,
    modulo_periodicity_score,
    refine_candidate_period,
    select_candidate_periods,
)
from tandemx.discover.rust_backend import RustBackendUnavailable, scan_read_for_periods
from tandemx.simulate.toy import reverse_complement, wrap_sequence


@dataclass(frozen=True)
class FastaRecord:
    read_id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class CandidateRepeat:
    read_id: str
    candidate_id: str
    sequence: str
    read_start: int
    read_end: int
    strand: str
    period_bp: int
    repeat_span_bp: int
    unit_count: float
    score: float
    low_complexity_flag: bool
    confidence: str
    warning: str


@dataclass(frozen=True)
class RepeatFamily:
    family_id: str
    monomer_id: str
    monomer_sequence: str
    monomer_length_bp: int
    support_read_count: int
    support_span_bp: int
    mean_identity: float
    low_complexity_flag: bool
    confidence: str
    warning: str


@dataclass(frozen=True)
class FamilySimilarity:
    family_a: str
    family_b: str
    length_a_bp: int
    length_b_bp: int
    kmer_jaccard: float
    shared_kmer_fraction: float
    local_identity: float
    local_overlap_bp: int
    local_overlap_fraction_shorter: float
    length_ratio: float
    orientation: str
    relationship: str
    redundant_candidate: bool
    notes: str


@dataclass(frozen=True)
class FamilyCollapse:
    original_family_id: str
    retained_family_id: str
    action: str
    reason: str
    relationship: str
    similarity_metrics: str
    notes: str


@dataclass(frozen=True)
class DiscoverConfig:
    reads: Path
    outdir: Path
    min_monomer_len: int
    max_monomer_len: int
    min_support_reads: int
    min_repeat_span: int
    min_read_length: int = 1
    kmer_size: int = 11
    top_periods: int = 5
    min_seed_occurrences: int = 2
    min_spacing_support: int = 2
    max_pairs_per_kmer: int = 100
    max_reads: int | None = None
    max_read_bases: int | None = None
    sample_rate: float = 1.0
    seed: int = 1
    progress_every: int = 1000
    chunk_size: int = 1000
    kmer_backend: str = "python"
    collapse_redundant_families: bool = False


def discover_toy_repeats(
    config: DiscoverConfig,
    logger: logging.Logger | None = None,
) -> tuple[list[CandidateRepeat], list[RepeatFamily]]:
    """Discover tandem repeat families with a bounded k-mer spacing prefilter."""
    validate_discover_config(config)
    config.outdir.mkdir(parents=True, exist_ok=True)
    logger = logger or logging.getLogger("tandemx.discover")
    candidate_path = config.outdir / "candidate_reads.tsv"
    for stale_path in (
        config.outdir / "monomers.fa",
        config.outdir / "families.tsv",
        config.outdir / "family_similarity.tsv",
        config.outdir / "collapsed_families.tsv",
        config.outdir / "collapsed_monomers.fa",
        config.outdir / "family_collapse.tsv",
    ):
        stale_path.unlink(missing_ok=True)

    candidates: list[CandidateRepeat] = []
    processed_reads = 0
    processed_bases = 0
    skipped_short_reads = 0
    skipped_short_kmer = 0
    skipped_low_complexity = 0
    seed_overflow_count = 0
    started = time.perf_counter()
    rng = random.Random(config.seed)

    with candidate_path.open("wt", encoding="utf-8") as handle:
        handle.write(candidate_reads_header() + "\n")
        handle.flush()
        for record in read_fasta(config.reads):
            if config.max_reads is not None and processed_reads >= config.max_reads:
                break
            if config.sample_rate < 1.0 and rng.random() > config.sample_rate:
                continue
            if (
                config.max_read_bases is not None
                and processed_bases + len(record.sequence) > config.max_read_bases
            ):
                logger.info(
                    "limit_reached=max_read_bases configured_bases=%s next_read_bases=%s",
                    config.max_read_bases,
                    len(record.sequence),
                )
                break

            processed_reads += 1
            processed_bases += len(record.sequence)
            candidate: CandidateRepeat | None = None
            if len(record.sequence) < config.min_read_length:
                skipped_short_reads += 1
            elif len(record.sequence) < config.kmer_size:
                skipped_short_kmer += 1
            elif is_low_complexity(record.sequence):
                skipped_low_complexity += 1
            else:
                finder = (
                    find_best_periodic_candidate_rust
                    if config.kmer_backend == "rust"
                    else find_best_periodic_candidate_with_stats
                )
                try:
                    candidate, overflow_count = finder(
                        record,
                        min_period=config.min_monomer_len,
                        max_period=config.max_monomer_len,
                        min_repeat_span=config.min_repeat_span,
                        candidate_index=len(candidates) + 1,
                        kmer_size=config.kmer_size,
                        top_periods=config.top_periods,
                        min_seed_occurrences=config.min_seed_occurrences,
                        min_spacing_support=config.min_spacing_support,
                        max_pairs_per_kmer=config.max_pairs_per_kmer,
                    )
                except RustBackendUnavailable as exc:
                    logger.error("rust_backend_unavailable=%s", exc)
                    raise ValueError(str(exc)) from exc
                seed_overflow_count += overflow_count

            if candidate is not None:
                candidates.append(candidate)
                handle.write(format_candidate_read(candidate) + "\n")
                handle.flush()

            if processed_reads % config.progress_every == 0:
                log_discover_progress(
                    logger,
                    processed_reads,
                    processed_bases,
                    len(candidates),
                    started,
                    config.max_reads,
                    config.max_read_bases,
                )

        log_discover_progress(
            logger,
            processed_reads,
            processed_bases,
            len(candidates),
            started,
            config.max_reads,
            config.max_read_bases,
        )
        logger.info(
            "filter_summary skipped_short_reads=%s skipped_short_kmer=%s "
            "skipped_low_complexity=%s seed_overflow_count=%s",
            skipped_short_reads,
            skipped_short_kmer,
            skipped_low_complexity,
            seed_overflow_count,
        )

    if not candidates:
        if processed_reads and skipped_short_kmer == processed_reads:
            raise ValueError(
                f"No reads were long enough for --kmer-size {config.kmer_size}; "
                f"processed {processed_reads} reads"
            )
        if processed_reads and skipped_low_complexity == processed_reads:
            raise ValueError(
                f"Only low-complexity reads were found; processed {processed_reads} reads"
            )
        raise ValueError(
            "No tandem repeat candidates found after k-mer spacing prefilter; "
            f"processed {processed_reads} reads"
        )
    if all(candidate.low_complexity_flag for candidate in candidates):
        raise ValueError("Only low-complexity tandem candidates found in reads")

    families = cluster_candidates(candidates, config.min_support_reads)
    if not families:
        raise ValueError(
            "No repeat families passed the minimum support threshold; "
            "try lowering --min-support-reads for toy data or check the reads"
        )
    similarities = compare_families(families, k=config.kmer_size)
    families = annotate_family_redundancy(families, similarities)
    write_monomers(config.outdir / "monomers.fa", families)
    write_families(config.outdir / "families.tsv", families)
    write_family_similarity(config.outdir / "family_similarity.tsv", similarities)
    if config.collapse_redundant_families:
        collapsed_families, collapse_records = collapse_redundant_families(families, similarities)
        write_monomers(config.outdir / "collapsed_monomers.fa", collapsed_families)
        write_families(config.outdir / "collapsed_families.tsv", collapsed_families)
        write_family_collapse(config.outdir / "family_collapse.tsv", collapse_records)
    return candidates, families


def validate_discover_config(config: DiscoverConfig) -> None:
    if config.min_monomer_len <= 0:
        raise ValueError("--min-period must be positive")
    if config.max_monomer_len < config.min_monomer_len:
        raise ValueError("--max-period must be greater than or equal to --min-period")
    if config.min_read_length <= 0:
        raise ValueError("--min-read-length must be positive")
    if config.kmer_size <= 0:
        raise ValueError("--kmer-size must be positive")
    if config.top_periods <= 0:
        raise ValueError("--top-periods must be positive")
    if config.min_seed_occurrences < 2:
        raise ValueError("--min-seed-occurrences must be at least 2")
    if config.min_spacing_support <= 0:
        raise ValueError("--min-spacing-support must be positive")
    if config.max_pairs_per_kmer <= 0:
        raise ValueError("--max-pairs-per-kmer must be positive")
    if config.max_reads is not None and config.max_reads <= 0:
        raise ValueError("--max-reads must be positive when provided")
    if config.max_read_bases is not None and config.max_read_bases <= 0:
        raise ValueError("--max-read-bases must be positive when provided")
    if not 0 < config.sample_rate <= 1:
        raise ValueError("--sample-rate must be greater than 0 and at most 1")
    if config.progress_every <= 0:
        raise ValueError("--progress-every must be positive")
    if config.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive")
    if config.kmer_backend not in {"python", "rust"}:
        raise ValueError("--kmer-backend must be python or rust")
    if config.kmer_backend == "rust" and config.kmer_size > 31:
        raise ValueError("Rust backend requires --kmer-size at most 31")


def read_fasta(path: Path) -> Iterable[FastaRecord]:
    try:
        for record in read_sequence_records(path):
            yield FastaRecord(
                read_id=record.id,
                description=record.description,
                sequence=record.sequence,
            )
    except SequenceFormatError:
        raise


def make_fasta_record(header: str, sequence_parts: Sequence[str]) -> FastaRecord:
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"Invalid FASTA: empty sequence for record {header}")
    read_id = header.split()[0].split(";")[0]
    return FastaRecord(read_id=read_id, description=header, sequence=sequence)


def find_best_periodic_candidate(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
    kmer_size: int = 11,
    top_periods: int = 5,
    min_seed_occurrences: int = 2,
    min_spacing_support: int = 2,
    max_pairs_per_kmer: int = 100,
) -> CandidateRepeat | None:
    candidate, _ = find_best_periodic_candidate_with_stats(
        record,
        min_period=min_period,
        max_period=max_period,
        min_repeat_span=min_repeat_span,
        candidate_index=candidate_index,
        kmer_size=kmer_size,
        top_periods=top_periods,
        min_seed_occurrences=min_seed_occurrences,
        min_spacing_support=min_spacing_support,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    return candidate


def find_best_periodic_candidate_with_stats(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
    kmer_size: int,
    top_periods: int,
    min_seed_occurrences: int,
    min_spacing_support: int,
    max_pairs_per_kmer: int,
) -> tuple[CandidateRepeat | None, int]:
    sequence = record.sequence.upper()
    if len(sequence) < max(min_repeat_span, min_period * 2):
        return None, 0
    bounded_max_period = min(max_period, len(sequence) // 2)
    if bounded_max_period < min_period or len(sequence) < kmer_size:
        return None, 0

    repeated_positions, overflow_count = extract_repeated_kmer_positions(
        sequence,
        kmer_size,
        min_seed_occurrences=min_seed_occurrences,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    histogram = build_spacing_histogram(
        repeated_positions,
        min_period=min_period,
        max_period=bounded_max_period,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    candidate_periods = select_candidate_periods(
        histogram,
        min_period=min_period,
        max_period=bounded_max_period,
        top_periods=top_periods,
        min_spacing_support=min_spacing_support,
    )
    if not candidate_periods:
        return None, overflow_count

    best_period, best_score = refine_candidate_period(
        sequence,
        repeated_positions,
        candidate_periods,
        min_period=min_period,
        max_period=bounded_max_period,
    )

    return build_candidate_from_period(
        record,
        sequence,
        best_period,
        best_score,
        candidate_index,
    ), overflow_count


def find_best_periodic_candidate_rust(
    record: FastaRecord,
    min_period: int,
    max_period: int,
    min_repeat_span: int,
    candidate_index: int,
    kmer_size: int,
    top_periods: int,
    min_seed_occurrences: int,
    min_spacing_support: int,
    max_pairs_per_kmer: int,
) -> tuple[CandidateRepeat | None, int]:
    sequence = record.sequence.upper()
    if len(sequence) < max(min_repeat_span, min_period * 2):
        return None, 0
    bounded_max_period = min(max_period, len(sequence) // 2)
    if bounded_max_period < min_period or len(sequence) < kmer_size:
        return None, 0
    result = scan_read_for_periods(
        sequence,
        k=kmer_size,
        min_period=min_period,
        max_period=bounded_max_period,
        top_periods=top_periods,
        min_seed_occurrences=min_seed_occurrences,
        min_spacing_support=min_spacing_support,
        max_pairs_per_kmer=max_pairs_per_kmer,
    )
    return build_candidate_from_period(
        record,
        sequence,
        result.best_period,
        result.periodicity_score,
        candidate_index,
    ), result.overflow_count


def build_candidate_from_period(
    record: FastaRecord,
    sequence: str,
    best_period: int,
    best_score: float,
    candidate_index: int,
) -> CandidateRepeat | None:
    if best_period <= 0 or best_score < 0.75:
        return None
    repeat_span = len(sequence)
    warning = ""
    low_complexity = is_low_complexity(sequence[:best_period])
    if low_complexity:
        warning = "low_complexity_candidate"
    confidence = "high" if best_score >= 0.9 else "medium"
    return CandidateRepeat(
        read_id=record.read_id,
        candidate_id=f"TXC{candidate_index:06d}",
        sequence=sequence[:best_period],
        read_start=0,
        read_end=len(sequence),
        strand=parse_strand(record.description),
        period_bp=best_period,
        repeat_span_bp=repeat_span,
        unit_count=repeat_span / best_period,
        score=best_score,
        low_complexity_flag=low_complexity,
        confidence=confidence,
        warning=warning,
    )


def periodicity_score(sequence: str, period: int) -> float:
    compared = len(sequence) - period
    if compared <= 0:
        return 0.0
    matches = 0
    valid = 0
    for index in range(compared):
        left = sequence[index]
        right = sequence[index + period]
        if left == "N" or right == "N":
            continue
        valid += 1
        if left == right:
            matches += 1
    if valid == 0:
        return 0.0
    return matches / valid


def parse_strand(description: str) -> str:
    for part in description.split(";"):
        if part.startswith("strand="):
            value = part.split("=", 1)[1]
            if value in {"+", "-"}:
                return value
    return "."


def is_low_complexity(sequence: str) -> bool:
    if not sequence:
        return True
    if len(sequence) <= 2:
        return True
    counts = {base: sequence.count(base) for base in "ACGT"}
    if max(counts.values()) / len(sequence) >= 0.8:
        return True
    first_dinucleotide = sequence[:2]
    for index in range(2, len(sequence) - 1, 2):
        if sequence[index : index + 2] != first_dinucleotide:
            return False
    return True


def cluster_candidates(
    candidates: Sequence[CandidateRepeat],
    min_support_reads: int,
    period_tolerance_bp: int = 5,
) -> list[RepeatFamily]:
    clusters: list[list[CandidateRepeat]] = []
    for candidate in sorted(candidates, key=lambda item: item.period_bp):
        placed = False
        for cluster in clusters:
            center = round(mean(item.period_bp for item in cluster))
            if abs(candidate.period_bp - center) <= period_tolerance_bp:
                cluster.append(candidate)
                placed = True
                break
        if not placed:
            clusters.append([candidate])

    supported = [cluster for cluster in clusters if len({item.read_id for item in cluster}) >= min_support_reads]
    supported.sort(key=lambda cluster: (-len({item.read_id for item in cluster}), round(mean(item.period_bp for item in cluster))))

    families = []
    for index, cluster in enumerate(supported, start=1):
        representative = max(cluster, key=lambda item: (item.score, item.repeat_span_bp))
        monomer_sequence = orient_monomer(representative.sequence)
        mean_identity = mean(item.score for item in cluster)
        low_complexity = any(item.low_complexity_flag for item in cluster)
        warning = "low_complexity_family" if low_complexity else ""
        confidence = "high" if len({item.read_id for item in cluster}) >= max(3, min_support_reads) and mean_identity >= 0.9 else "medium"
        families.append(
            RepeatFamily(
                family_id=f"TXF{index:06d}",
                monomer_id=f"TXM{index:06d}",
                monomer_sequence=monomer_sequence,
                monomer_length_bp=representative.period_bp,
                support_read_count=len({item.read_id for item in cluster}),
                support_span_bp=sum(item.repeat_span_bp for item in cluster),
                mean_identity=mean_identity,
                low_complexity_flag=low_complexity,
                confidence=confidence,
                warning=warning,
            )
        )
    return families


def orient_monomer(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def compare_families(families: Sequence[RepeatFamily], k: int = 11) -> list[FamilySimilarity]:
    """Compare representative monomers to flag possible family redundancy."""
    similarities: list[FamilySimilarity] = []
    for index, family_a in enumerate(families):
        for family_b in families[index + 1 :]:
            similarities.append(compare_family_pair(family_a, family_b, k=k))
    return similarities


def compare_family_pair(family_a: RepeatFamily, family_b: RepeatFamily, k: int = 11) -> FamilySimilarity:
    kmers_a = canonical_kmer_set(family_a.monomer_sequence, k)
    kmers_b = canonical_kmer_set(family_b.monomer_sequence, k)
    shared = len(kmers_a & kmers_b)
    union = len(kmers_a | kmers_b)
    kmer_jaccard = shared / union if union else 0.0
    shared_fraction = shared / min(len(kmers_a), len(kmers_b)) if kmers_a and kmers_b else 0.0
    identity, overlap, orientation = best_local_identity(
        family_a.monomer_sequence,
        family_b.monomer_sequence,
    )
    shorter = max(1, min(family_a.monomer_length_bp, family_b.monomer_length_bp))
    longer = max(family_a.monomer_length_bp, family_b.monomer_length_bp)
    overlap_fraction_shorter = overlap / shorter
    length_ratio = longer / shorter
    relationship, redundant_candidate, notes = classify_family_relationship(
        kmer_jaccard=kmer_jaccard,
        shared_kmer_fraction=shared_fraction,
        local_identity=identity,
        local_overlap_fraction_shorter=overlap_fraction_shorter,
        length_ratio=length_ratio,
    )
    return FamilySimilarity(
        family_a=family_a.family_id,
        family_b=family_b.family_id,
        length_a_bp=family_a.monomer_length_bp,
        length_b_bp=family_b.monomer_length_bp,
        kmer_jaccard=kmer_jaccard,
        shared_kmer_fraction=shared_fraction,
        local_identity=identity,
        local_overlap_bp=overlap,
        local_overlap_fraction_shorter=overlap_fraction_shorter,
        length_ratio=length_ratio,
        orientation=orientation,
        relationship=relationship,
        redundant_candidate=redundant_candidate,
        notes=notes,
    )


def canonical_kmer_set(sequence: str, k: int) -> set[str]:
    if k <= 0 or len(sequence) < k:
        return set()
    return {canonical_kmer(sequence[index : index + k]) for index in range(len(sequence) - k + 1)}


def best_local_identity(sequence_a: str, sequence_b: str) -> tuple[float, int, str]:
    """Return the best ungapped local identity over both orientations."""
    min_overlap = min(50, len(sequence_a), len(sequence_b))
    best_identity = 0.0
    best_overlap = 0
    best_orientation = "forward"
    for orientation, oriented_b in (
        ("forward", sequence_b),
        ("reverse", reverse_complement(sequence_b)),
    ):
        for offset in range(-len(oriented_b) + 1, len(sequence_a)):
            start_a = max(0, offset)
            start_b = max(0, -offset)
            overlap = min(len(sequence_a) - start_a, len(oriented_b) - start_b)
            if overlap < min_overlap:
                continue
            matches = sum(
                1
                for base_index in range(overlap)
                if sequence_a[start_a + base_index] == oriented_b[start_b + base_index]
            )
            identity = matches / overlap
            if identity > best_identity or (identity == best_identity and overlap > best_overlap):
                best_identity = identity
                best_overlap = overlap
                best_orientation = orientation
    return best_identity, best_overlap, best_orientation


def classify_family_relationship(
    *,
    kmer_jaccard: float,
    shared_kmer_fraction: float,
    local_identity: float,
    local_overlap_fraction_shorter: float,
    length_ratio: float,
) -> tuple[str, bool, str]:
    near_integer_ratio = abs(length_ratio - round(length_ratio)) <= 0.05 and round(length_ratio) >= 2
    if (
        near_integer_ratio
        and local_identity >= 0.75
        and local_overlap_fraction_shorter >= 0.75
        and (kmer_jaccard >= 0.25 or shared_kmer_fraction >= 0.5)
    ):
        return (
            "possible_higher_order_or_partial",
            False,
            "Moderate or high similarity with near-integer length ratio; inspect as possible higher-order unit, dimer, or related family.",
        )
    if (
        local_identity >= 0.9
        and local_overlap_fraction_shorter >= 0.85
        and shared_kmer_fraction >= 0.8
        and length_ratio <= 1.20
    ):
        return (
            "likely_redundant",
            True,
            "High sequence similarity; shorter or lower-support representative may be redundant.",
        )
    if (
        local_identity >= 0.75
        and local_overlap_fraction_shorter >= 0.75
        and (kmer_jaccard >= 0.25 or shared_kmer_fraction >= 0.5)
    ):
        return (
            "possible_higher_order_or_partial",
            False,
            "Moderate local similarity; inspect as possible higher-order unit, partial duplicate, or related family.",
        )
    return "distinct", False, ""


def annotate_family_redundancy(
    families: Sequence[RepeatFamily],
    similarities: Sequence[FamilySimilarity],
) -> list[RepeatFamily]:
    warnings_by_family: dict[str, list[str]] = {family.family_id: [] for family in families}
    for similarity in similarities:
        if similarity.relationship == "distinct":
            continue
        warning = f"{similarity.relationship}:{similarity.family_a}-{similarity.family_b}"
        warnings_by_family[similarity.family_a].append(warning)
        warnings_by_family[similarity.family_b].append(warning)

    annotated: list[RepeatFamily] = []
    for family in families:
        warnings = [family.warning] if family.warning else []
        warnings.extend(warnings_by_family[family.family_id])
        annotated.append(replace(family, warning=";".join(warnings)))
    return annotated


def collapse_redundant_families(
    families: Sequence[RepeatFamily],
    similarities: Sequence[FamilySimilarity],
) -> tuple[list[RepeatFamily], list[FamilyCollapse]]:
    """Collapse only pairs classified as likely redundant."""
    family_by_id = {family.family_id: family for family in families}
    parent = {family.family_id: family.family_id for family in families}

    def find(family_id: str) -> str:
        while parent[family_id] != family_id:
            parent[family_id] = parent[parent[family_id]]
            family_id = parent[family_id]
        return family_id

    def union(left: str, right: str) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for similarity in similarities:
        if similarity.relationship == "likely_redundant":
            union(similarity.family_a, similarity.family_b)

    components: dict[str, list[RepeatFamily]] = {}
    for family in families:
        components.setdefault(find(family.family_id), []).append(family)

    retained_by_family: dict[str, str] = {}
    for component in components.values():
        retained = choose_retained_family(component)
        for family in component:
            retained_by_family[family.family_id] = retained.family_id

    similarity_lookup = {
        frozenset((similarity.family_a, similarity.family_b)): similarity
        for similarity in similarities
        if similarity.relationship == "likely_redundant"
    }
    records: list[FamilyCollapse] = []
    collapsed: list[RepeatFamily] = []
    for family in families:
        retained_id = retained_by_family[family.family_id]
        if family.family_id == retained_id:
            collapsed.append(family)
            records.append(
                FamilyCollapse(
                    original_family_id=family.family_id,
                    retained_family_id=retained_id,
                    action="retained",
                    reason="selected_representative" if len(components[find(family.family_id)]) > 1 else "no_likely_redundant_match",
                    relationship="",
                    similarity_metrics="",
                    notes="Possible higher-order or partial relationships are not collapsed.",
                )
            )
            continue
        similarity = similarity_lookup.get(frozenset((family.family_id, retained_id)))
        records.append(
            FamilyCollapse(
                original_family_id=family.family_id,
                retained_family_id=retained_id,
                action="collapsed",
                reason="likely_redundant_to_retained_family",
                relationship=similarity.relationship if similarity else "likely_redundant",
                similarity_metrics=format_similarity_metrics(similarity) if similarity else "",
                notes="Collapsed only because the relationship was likely_redundant.",
            )
        )
    return collapsed, records


def choose_retained_family(families: Sequence[RepeatFamily]) -> RepeatFamily:
    return max(
        families,
        key=lambda family: (
            family.support_read_count,
            family.support_span_bp,
            family.mean_identity,
            -family.monomer_length_bp,
            family.family_id,
        ),
    )


def format_similarity_metrics(similarity: FamilySimilarity | None) -> str:
    if similarity is None:
        return ""
    return (
        f"kmer_jaccard={similarity.kmer_jaccard:.4f};"
        f"shared_kmer_fraction={similarity.shared_kmer_fraction:.4f};"
        f"local_identity={similarity.local_identity:.4f};"
        f"local_overlap_bp={similarity.local_overlap_bp};"
        f"length_ratio={similarity.length_ratio:.4f}"
    )


def sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def log_discover_progress(
    logger: logging.Logger,
    processed_reads: int,
    processed_bases: int,
    candidate_reads: int,
    started: float,
    max_reads: int | None,
    max_read_bases: int | None,
) -> None:
    elapsed = max(time.perf_counter() - started, 1e-9)
    reads_per_second = processed_reads / elapsed
    mb_per_second = (processed_bases / 1_000_000) / elapsed
    remaining_estimates: list[float] = []
    if max_reads is not None and reads_per_second > 0:
        remaining_estimates.append(max(0, max_reads - processed_reads) / reads_per_second)
    if max_read_bases is not None and mb_per_second > 0:
        remaining_mb = max(0, max_read_bases - processed_bases) / 1_000_000
        remaining_estimates.append(remaining_mb / mb_per_second)
    estimated_remaining = (
        f"{min(remaining_estimates):.1f}" if remaining_estimates else "unknown"
    )
    logger.info(
        "progress processed_reads=%s processed_bases=%s candidate_reads=%s "
        "elapsed_seconds=%.3f reads_per_second=%.3f mb_per_second=%.3f "
        "estimated_remaining_seconds=%s",
        processed_reads,
        processed_bases,
        candidate_reads,
        elapsed,
        reads_per_second,
        mb_per_second,
        estimated_remaining,
    )


def candidate_reads_header() -> str:
    return (
        "read_id\tcandidate_id\tread_start\tread_end\tstrand\tperiod_bp\t"
        "repeat_span_bp\tunit_count\tscore\tlow_complexity_flag\tconfidence\twarning"
    )


def format_candidate_read(candidate: CandidateRepeat) -> str:
    return "\t".join(
        [
            candidate.read_id,
            candidate.candidate_id,
            str(candidate.read_start),
            str(candidate.read_end),
            candidate.strand,
            str(candidate.period_bp),
            str(candidate.repeat_span_bp),
            f"{candidate.unit_count:.3f}",
            f"{candidate.score:.4f}",
            str(candidate.low_complexity_flag).lower(),
            candidate.confidence,
            candidate.warning,
        ]
    )


def write_candidate_reads(path: Path, candidates: Sequence[CandidateRepeat]) -> None:
    lines = [candidate_reads_header()]
    lines.extend(format_candidate_read(candidate) for candidate in candidates)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_monomers(path: Path, families: Sequence[RepeatFamily]) -> None:
    lines = []
    for family in families:
        lines.append(
            (
                f">family_id={family.family_id};monomer_id={family.monomer_id};"
                f"length_bp={family.monomer_length_bp};confidence={family.confidence}"
            )
        )
        lines.append(wrap_sequence(family.monomer_sequence))
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_families(path: Path, families: Sequence[RepeatFamily]) -> None:
    lines = [
        (
            "family_id\tmonomer_id\tmonomer_length_bp\tconsensus_md5\tgc_fraction\t"
            "support_read_count\tsupport_span_bp\tmean_identity\tlow_complexity_flag\t"
            "confidence\twarning"
        )
    ]
    for family in families:
        gc_fraction = (family.monomer_sequence.count("G") + family.monomer_sequence.count("C")) / len(family.monomer_sequence)
        lines.append(
            "\t".join(
                [
                    family.family_id,
                    family.monomer_id,
                    str(family.monomer_length_bp),
                    sequence_md5(family.monomer_sequence),
                    f"{gc_fraction:.4f}",
                    str(family.support_read_count),
                    str(family.support_span_bp),
                    f"{family.mean_identity:.4f}",
                    str(family.low_complexity_flag).lower(),
                    family.confidence,
                    family.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def family_similarity_header() -> str:
    return (
        "family_a\tfamily_b\tlength_a_bp\tlength_b_bp\tkmer_jaccard\t"
        "shared_kmer_fraction\tlocal_identity\tlocal_overlap_bp\t"
        "local_overlap_fraction_shorter\tlength_ratio\torientation\t"
        "relationship\tredundant_candidate\tnotes"
    )


def format_family_similarity(similarity: FamilySimilarity) -> str:
    return "\t".join(
        [
            similarity.family_a,
            similarity.family_b,
            str(similarity.length_a_bp),
            str(similarity.length_b_bp),
            f"{similarity.kmer_jaccard:.4f}",
            f"{similarity.shared_kmer_fraction:.4f}",
            f"{similarity.local_identity:.4f}",
            str(similarity.local_overlap_bp),
            f"{similarity.local_overlap_fraction_shorter:.4f}",
            f"{similarity.length_ratio:.4f}",
            similarity.orientation,
            similarity.relationship,
            str(similarity.redundant_candidate).lower(),
            similarity.notes,
        ]
    )


def write_family_similarity(path: Path, similarities: Sequence[FamilySimilarity]) -> None:
    lines = [family_similarity_header()]
    lines.extend(format_family_similarity(similarity) for similarity in similarities)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def family_collapse_header() -> str:
    return "original_family_id\tretained_family_id\taction\treason\trelationship\tsimilarity_metrics\tnotes"


def format_family_collapse(record: FamilyCollapse) -> str:
    return "\t".join(
        [
            record.original_family_id,
            record.retained_family_id,
            record.action,
            record.reason,
            record.relationship,
            record.similarity_metrics,
            record.notes,
        ]
    )


def write_family_collapse(path: Path, records: Sequence[FamilyCollapse]) -> None:
    lines = [family_collapse_header()]
    lines.extend(format_family_collapse(record) for record in records)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
