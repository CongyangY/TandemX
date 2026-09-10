#!/usr/bin/env python3
"""Score a discovered catalogue against a fixed known-repeat clone library.

This is a post hoc exclusion analysis.  It never supplies known sequences to
discovery, and a no-match result applies only to the named library.  For each
pair, the complete shorter sequence is aligned to the doubled longer sequence
in both orientations.  Doubling permits a circular boundary match without
enumerating rotations; edlib ``HW`` keeps the shorter sequence end-to-end.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import edlib

from benchmarks.scripts.integrate_tr_candidate_evidence import read_fasta
from tandemx.annotation import RepeatRecord, read_discovered_catalog, read_known_repeats
from tandemx.simulate.toy import reverse_complement


PAIR_FIELDS = (
    "family_id",
    "family_length_bp",
    "known_id",
    "known_length_bp",
    "shorter_length_bp",
    "best_orientation",
    "edit_distance",
    "glocal_edit_identity",
    "length_ratio",
    "integer_multiple",
    "multiple_relative_error",
    "exclusion_state",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_masked_known_fasta(path: Path) -> tuple[list[RepeatRecord], int]:
    """Read a fixed library and mask non-ACGT IUPAC characters as N."""
    records: list[RepeatRecord] = []
    masked_count = 0
    for identifier, sequence in read_fasta(path).items():
        masked = []
        for base in sequence.upper():
            if base in "ACGT":
                masked.append(base)
            else:
                masked.append("N")
                masked_count += 1
        records.append(RepeatRecord(identifier=identifier, sequence="".join(masked)))
    return records, masked_count


def circular_glocal_identity(first: str, second: str) -> tuple[float, int, str]:
    """Align the complete shorter sequence across a circular longer sequence."""
    shorter, longer = (first, second) if len(first) <= len(second) else (second, first)
    query = shorter.upper().replace("N", "X")
    target = (longer.upper() + longer.upper()).replace("N", "Y")
    best: tuple[int, str] | None = None
    for orientation, oriented in (
        ("forward", query),
        ("reverse", reverse_complement(query)),
    ):
        distance = edlib.align(oriented, target, mode="HW", task="distance")["editDistance"]
        if distance < 0:
            raise RuntimeError("edlib did not return a glocal edit distance")
        if best is None or distance < best[0]:
            best = (distance, orientation)
    assert best is not None
    return 1.0 - best[0] / len(shorter), best[0], best[1]


def classify(identity: float, shorter_length: int, possible_threshold: float, strong_threshold: float,
             minimum_aligned_bp: int) -> str:
    if shorter_length < minimum_aligned_bp:
        return "insufficient_length_for_exclusion"
    if identity >= strong_threshold:
        return "strong_known_clone_match"
    if identity >= possible_threshold:
        return "possible_known_family_match"
    if identity >= 0.70:
        return "weak_known_similarity"
    return "no_match_in_limited_library"


def score_pair(first: str, second: str, *, possible_threshold: float, strong_threshold: float,
               minimum_aligned_bp: int) -> dict[str, object]:
    identity, distance, orientation = circular_glocal_identity(first, second)
    shorter = min(len(first), len(second))
    longer = max(len(first), len(second))
    multiple = max(1, round(longer / shorter))
    relative_error = abs(longer - multiple * shorter) / longer
    return {
        "shorter_length_bp": shorter,
        "best_orientation": orientation,
        "edit_distance": distance,
        "glocal_edit_identity": identity,
        "length_ratio": longer / shorter,
        "integer_multiple": multiple if multiple >= 2 and relative_error <= 0.10 else "NA",
        "multiple_relative_error": relative_error,
        "exclusion_state": classify(
            identity,
            shorter,
            possible_threshold,
            strong_threshold,
            minimum_aligned_bp,
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--known", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--possible-threshold", type=float, default=0.80)
    parser.add_argument("--strong-threshold", type=float, default=0.90)
    parser.add_argument("--minimum-aligned-bp", type=int, default=50)
    parser.add_argument(
        "--mask-ambiguous-known-bases",
        action="store_true",
        help="Mask non-ACGT characters in the known FASTA as N and record the count.",
    )
    args = parser.parse_args()
    if not 0 < args.possible_threshold <= args.strong_threshold <= 1:
        parser.error("require 0 < possible-threshold <= strong-threshold <= 1")
    if args.minimum_aligned_bp < 1:
        parser.error("minimum-aligned-bp must be positive")

    catalog = read_discovered_catalog(args.catalog)
    if args.mask_ambiguous_known_bases:
        known, masked_known_base_count = read_masked_known_fasta(args.known)
    else:
        known = read_known_repeats(args.known)
        masked_known_base_count = 0
    if not catalog or not known:
        parser.error("catalog and known library must both contain sequences")
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for family in catalog:
        for reference in known:
            row = {
                "family_id": family.identifier,
                "family_length_bp": len(family.sequence),
                "known_id": reference.identifier,
                "known_length_bp": len(reference.sequence),
            }
            row.update(
                score_pair(
                    family.sequence,
                    reference.sequence,
                    possible_threshold=args.possible_threshold,
                    strong_threshold=args.strong_threshold,
                    minimum_aligned_bp=args.minimum_aligned_bp,
                )
            )
            rows.append(row)

    pair_path = args.outdir / "known_repeat_pairwise.tsv"
    with pair_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            rendered = dict(row)
            rendered["glocal_edit_identity"] = f"{float(row['glocal_edit_identity']):.6f}"
            rendered["length_ratio"] = f"{float(row['length_ratio']):.6f}"
            rendered["multiple_relative_error"] = f"{float(row['multiple_relative_error']):.6f}"
            writer.writerow(rendered)

    state_rank = {
        "strong_known_clone_match": 4,
        "possible_known_family_match": 3,
        "weak_known_similarity": 2,
        "insufficient_length_for_exclusion": 1,
        "no_match_in_limited_library": 0,
    }
    summary_path = args.outdir / "known_repeat_family_summary.tsv"
    summary_fields = (
        "family_id",
        "family_length_bp",
        "best_known_id",
        "best_known_length_bp",
        "shorter_length_bp",
        "best_orientation",
        "edit_distance",
        "glocal_edit_identity",
        "length_ratio",
        "integer_multiple",
        "multiple_relative_error",
        "exclusion_state",
        "interpretation_boundary",
    )
    best_rows: list[dict[str, object]] = []
    for family in catalog:
        candidates = [row for row in rows if row["family_id"] == family.identifier]
        best = max(
            candidates,
            key=lambda row: (
                state_rank[str(row["exclusion_state"])],
                float(row["glocal_edit_identity"]),
                int(row["shorter_length_bp"]),
                str(row["known_id"]),
            ),
        )
        best_rows.append(best)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in best_rows:
            writer.writerow(
                {
                    "family_id": row["family_id"],
                    "family_length_bp": row["family_length_bp"],
                    "best_known_id": row["known_id"],
                    "best_known_length_bp": row["known_length_bp"],
                    "shorter_length_bp": row["shorter_length_bp"],
                    "best_orientation": row["best_orientation"],
                    "edit_distance": row["edit_distance"],
                    "glocal_edit_identity": f"{float(row['glocal_edit_identity']):.6f}",
                    "length_ratio": f"{float(row['length_ratio']):.6f}",
                    "integer_multiple": row["integer_multiple"],
                    "multiple_relative_error": f"{float(row['multiple_relative_error']):.6f}",
                    "exclusion_state": row["exclusion_state"],
                    "interpretation_boundary": "no-match is limited to the fixed input library; not proof of novelty",
                }
            )

    counts: dict[str, int] = {}
    for row in best_rows:
        state = str(row["exclusion_state"])
        counts[state] = counts.get(state, 0) + 1
    receipt = {
        "schema_version": 1,
        "catalog": str(args.catalog.resolve()),
        "catalog_sha256": sha256(args.catalog),
        "known": str(args.known.resolve()),
        "known_sha256": sha256(args.known),
        "catalog_count": len(catalog),
        "known_count": len(known),
        "pair_count": len(rows),
        "possible_threshold": args.possible_threshold,
        "strong_threshold": args.strong_threshold,
        "minimum_aligned_bp": args.minimum_aligned_bp,
        "mask_ambiguous_known_bases": args.mask_ambiguous_known_bases,
        "masked_known_base_count": masked_known_base_count,
        "alignment": "complete_shorter_vs_doubled_longer_edlib_HW_both_orientations",
        "edlib_version": version("edlib"),
        "family_state_counts": counts,
        "pairwise_sha256": sha256(pair_path),
        "summary_sha256": sha256(summary_path),
        "interpretation_boundary": "post hoc exclusion against this library only; no novelty claim",
    }
    (args.outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
