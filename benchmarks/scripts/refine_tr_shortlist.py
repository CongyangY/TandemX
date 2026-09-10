#!/usr/bin/env python3
"""Apply a fixed additional known-sequence exclusion to a TR shortlist."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.integrate_tr_candidate_evidence import read_fasta


ADDED_FIELDS = (
    "additional_library_name",
    "additional_best_known_id",
    "additional_best_known_identity",
    "additional_exclusion_state",
    "refined_candidate_state",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_additional_exclusion(state: str) -> str:
    if state == "no_match_in_limited_library":
        return "retained_preliminary_candidate"
    if state == "insufficient_length_for_exclusion":
        return "excluded_unresolved_additional_library"
    return "excluded_known_or_related_sequence"


def write_table(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--candidate-fasta", required=True, type=Path)
    parser.add_argument("--exclusion-summary", required=True, type=Path)
    parser.add_argument("--library-name", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    for path in (args.candidate_table, args.candidate_fasta, args.exclusion_summary):
        if not path.is_file():
            parser.error(f"input does not exist: {path}")

    with args.candidate_table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original_fields = list(reader.fieldnames or [])
        candidates = list(reader)
    with args.exclusion_summary.open(encoding="utf-8", newline="") as handle:
        exclusions = {row["family_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    sequences = read_fasta(args.candidate_fasta)
    candidate_ids = {row["family_id"] for row in candidates}
    if candidate_ids != set(sequences) or candidate_ids.difference(exclusions):
        parser.error("candidate table, FASTA, and exclusion summary do not contain matching candidate IDs")

    all_rows: list[dict[str, str]] = []
    retained: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    for candidate in candidates:
        family_id = candidate["family_id"]
        exclusion = exclusions[family_id]
        refined_state = classify_additional_exclusion(exclusion["exclusion_state"])
        row = {
            **candidate,
            "additional_library_name": args.library_name,
            "additional_best_known_id": exclusion["best_known_id"],
            "additional_best_known_identity": exclusion["glocal_edit_identity"],
            "additional_exclusion_state": exclusion["exclusion_state"],
            "refined_candidate_state": refined_state,
        }
        all_rows.append(row)
        (retained if refined_state == "retained_preliminary_candidate" else excluded).append(row)

    args.outdir.mkdir(parents=True, exist_ok=True)
    fields = original_fields + list(ADDED_FIELDS)
    all_path = args.outdir / "refined_shortlist_all.tsv"
    retained_path = args.outdir / "refined_shortlist_retained.tsv"
    excluded_path = args.outdir / "refined_shortlist_excluded.tsv"
    fasta_path = args.outdir / "refined_shortlist_retained.fa"
    write_table(all_path, fields, all_rows)
    write_table(retained_path, fields, retained)
    write_table(excluded_path, fields, excluded)
    with fasta_path.open("w", encoding="utf-8") as handle:
        for row in retained:
            family_id = row["family_id"]
            handle.write(f">{family_id} retained_preliminary_candidate_not_novel\n")
            sequence = sequences[family_id]
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")

    output_paths = (all_path, retained_path, excluded_path, fasta_path)
    receipt = {
        "schema_version": 1,
        "library_name": args.library_name,
        "candidate_count": len(candidates),
        "retained_count": len(retained),
        "excluded_count": len(excluded),
        "retained_family_ids": [row["family_id"] for row in retained],
        "excluded_family_ids": [row["family_id"] for row in excluded],
        "input_sha256": {
            args.candidate_table.name: sha256(args.candidate_table),
            args.candidate_fasta.name: sha256(args.candidate_fasta),
            args.exclusion_summary.name: sha256(args.exclusion_summary),
        },
        "output_sha256": {path.name: sha256(path) for path in output_paths},
        "interpretation_boundary": "additional exclusion only; retained rows remain preliminary and are not validated novel repeats",
    }
    (args.outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
