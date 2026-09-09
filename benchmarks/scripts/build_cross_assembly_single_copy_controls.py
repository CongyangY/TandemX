"""Build deterministic real-genome controls unique in two enrolled assemblies.

This validation-only helper does not alter TandemX public behavior. It samples
bounded candidate k-mers across the newer assembly, excludes every catalogue
monomer k-mer, and uses the native target-only counter to require exactly one
canonical occurrence in both old and new assemblies.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file
from tandemx.discover.rust_backend import RustDiagnosticKmerCounter
from tandemx.io.sequences import read_sequence_records
from tandemx.quantify.mvp import read_monomer_fasta
from tandemx.utils.kmers import (
    canonical_kmer,
    is_low_complexity_kmer,
    iter_circular_canonical_kmers,
)


def selection_score(word: str, seed: int) -> int:
    """Return a stable whole-word priority independent of Python hash state."""
    payload = f"{seed}\0{word}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest(), "big")


def catalogue_words(catalogue: Path, k: int) -> set[str]:
    """Return all canonical tandem-context k-mers in the frozen catalogue."""
    return {
        word
        for monomer in read_monomer_fasta(catalogue)
        for word in iter_circular_canonical_kmers(monomer.sequence, k)
    }


def candidate_words(
    assembly: Path,
    catalogue: Path,
    *,
    k: int,
    stride: int,
    candidate_limit: int,
    seed: int,
) -> tuple[list[str], dict[str, int]]:
    """Keep the lowest stable-hash background candidates with bounded memory."""
    if not 1 <= k <= 31:
        raise ValueError("k must be in 1..31 for the native target counter")
    if stride < 1 or candidate_limit < 1:
        raise ValueError("stride and candidate_limit must be positive")
    excluded = catalogue_words(catalogue, k)
    heap: list[tuple[int, str]] = []
    retained: set[str] = set()
    sampled_positions = invalid_words = low_complexity_words = catalogue_hits = 0

    for record in read_sequence_records(assembly):
        sequence = record.sequence.upper()
        for start in range(0, len(sequence) - k + 1, stride):
            sampled_positions += 1
            raw = sequence[start : start + k]
            if set(raw) - set("ACGT"):
                invalid_words += 1
                continue
            word = canonical_kmer(raw)
            if is_low_complexity_kmer(word):
                low_complexity_words += 1
                continue
            if word in excluded:
                catalogue_hits += 1
                continue
            if word in retained:
                continue
            score = selection_score(word, seed)
            entry = (-score, word)
            if len(heap) < candidate_limit:
                heapq.heappush(heap, entry)
                retained.add(word)
            elif score < -heap[0][0]:
                _, removed = heapq.heapreplace(heap, entry)
                retained.remove(removed)
                retained.add(word)

    ordered = sorted(retained, key=lambda word: (selection_score(word, seed), word))
    audit = {
        "sampled_positions": sampled_positions,
        "invalid_sampled_words": invalid_words,
        "low_complexity_sampled_words": low_complexity_words,
        "catalogue_overlap_sampled_words": catalogue_hits,
        "retained_candidates": len(ordered),
        "catalogue_kmers_excluded": len(excluded),
    }
    return ordered, audit


def count_targets(assembly: Path, k: int, targets: set[str]) -> dict[str, int]:
    """Count selected canonical targets while holding at most one contig in memory."""
    counter = RustDiagnosticKmerCounter(k, targets)
    for record in read_sequence_records(assembly):
        counter.count_sequence(record.sequence)
    observed = counter.counts()
    return {word: int(observed.get(word, 0)) for word in targets}


def build_controls(
    old_assembly: Path,
    new_assembly: Path,
    catalogue: Path,
    output_tsv: Path,
    receipt_json: Path,
    *,
    k: int,
    desired: int,
    stride: int,
    candidate_limit: int,
    seed: int,
) -> dict[str, object]:
    """Select, cross-count, write, and describe exact single-copy controls."""
    if desired < 1 or candidate_limit < desired:
        raise ValueError("candidate_limit must be at least desired and both must be positive")
    for path in (old_assembly, new_assembly, catalogue):
        if not path.is_file():
            raise ValueError(f"Required input is not a file: {path}")
    for path in (output_tsv, receipt_json):
        if path.exists():
            raise ValueError(f"Refusing to overwrite existing output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)

    candidates, candidate_audit = candidate_words(
        new_assembly,
        catalogue,
        k=k,
        stride=stride,
        candidate_limit=candidate_limit,
        seed=seed,
    )
    targets = set(candidates)
    old_counts = count_targets(old_assembly, k, targets)
    new_counts = count_targets(new_assembly, k, targets)
    eligible = [
        word for word in candidates if old_counts[word] == 1 and new_counts[word] == 1
    ]
    if len(eligible) < desired:
        raise ValueError(
            f"Only {len(eligible)} candidates are exactly single-copy in both assemblies; "
            f"requested {desired}"
        )
    selected = eligible[:desired]
    with output_tsv.open("w", encoding="utf-8") as handle:
        handle.write("kmer\texpected_copy_number\told_assembly_count\tnew_assembly_count\tselection_rank\n")
        for rank, word in enumerate(selected, start=1):
            handle.write(f"{word}\t1\t1\t1\t{rank}\n")

    receipt: dict[str, object] = {
        "schema_version": 1,
        "status": "complete",
        "k": k,
        "seed": seed,
        "selection_stride_bp_per_contig": stride,
        "candidate_limit": candidate_limit,
        "requested_controls": desired,
        "eligible_cross_assembly_single_copy_candidates": len(eligible),
        "selected_controls": len(selected),
        "candidate_audit": candidate_audit,
        "criterion": "canonical_kmer_exactly_once_in_each_old_and_new_assembly;absent_from_all_catalogue_monomer_tandem_contexts;not_low_complexity",
        "warning": "assembly-derived_controls_normalize_read_depth_but_do_not_make_either_assembly_copy_number_truth",
        "inputs": {
            "old_assembly": str(old_assembly.resolve()),
            "old_assembly_sha256": digest_file(old_assembly),
            "new_assembly": str(new_assembly.resolve()),
            "new_assembly_sha256": digest_file(new_assembly),
            "catalogue": str(catalogue.resolve()),
            "catalogue_sha256": digest_file(catalogue),
        },
        "output_tsv": str(output_tsv.resolve()),
    }
    receipt_json.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-assembly", required=True, type=Path)
    parser.add_argument("--new-assembly", required=True, type=Path)
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--receipt-json", required=True, type=Path)
    parser.add_argument("--k", required=True, type=int)
    parser.add_argument("--desired", type=int, default=20_000)
    parser.add_argument("--stride", type=int, default=251)
    parser.add_argument("--candidate-limit", type=int, default=200_000)
    parser.add_argument("--seed", type=int, default=6101)
    args = parser.parse_args()
    build_controls(
        args.old_assembly,
        args.new_assembly,
        args.catalogue,
        args.output_tsv,
        args.receipt_json,
        k=args.k,
        desired=args.desired,
        stride=args.stride,
        candidate_limit=args.candidate_limit,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
