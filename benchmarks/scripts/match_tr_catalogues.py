#!/usr/bin/env python3
"""Match fixed query TR consensuses against a second operational catalogue.

The complete shorter sequence is aligned to a doubled longer sequence with
edlib HW in both orientations.  The result measures cross-depth catalogue
recurrence; it is not an independent biological replication or novelty test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

from benchmarks.scripts.integrate_tr_candidate_evidence import read_fasta
from benchmarks.scripts.score_known_repeat_exclusion import circular_glocal_identity


FIELDS = (
    "query_id",
    "query_length_bp",
    "target_id",
    "target_length_bp",
    "best_orientation",
    "edit_distance",
    "glocal_edit_identity",
    "length_ratio",
    "nearest_integer_multiple",
    "multiple_relative_error",
    "match_state",
    "interpretation_boundary",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_match(
    identity: float,
    length_ratio: float,
    direct_identity: float,
    direct_minimum_ratio: float,
    direct_maximum_ratio: float,
    related_identity: float,
    maximum_multiple_error: float,
) -> tuple[str, int, float]:
    nearest_multiple = max(1, round(length_ratio))
    multiple_error = abs(length_ratio - nearest_multiple) / nearest_multiple
    if direct_identity <= identity and direct_minimum_ratio <= length_ratio <= direct_maximum_ratio:
        return "direct_catalogue_recurrence", nearest_multiple, multiple_error
    if identity >= related_identity and multiple_error <= maximum_multiple_error:
        return "related_period_multiple", nearest_multiple, multiple_error
    return "no_qualifying_match", nearest_multiple, multiple_error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--direct-identity", type=float, default=0.95)
    parser.add_argument("--direct-minimum-length-ratio", type=float, default=0.90)
    parser.add_argument("--direct-maximum-length-ratio", type=float, default=1.10)
    parser.add_argument("--related-identity", type=float, default=0.90)
    parser.add_argument("--maximum-multiple-error", type=float, default=0.05)
    args = parser.parse_args()
    for path in (args.query, args.target):
        if not path.is_file():
            parser.error(f"input does not exist: {path}")
    if not (0 <= args.related_identity <= args.direct_identity <= 1):
        parser.error("identity thresholds must satisfy 0 <= related <= direct <= 1")
    if not (0 < args.direct_minimum_length_ratio <= args.direct_maximum_length_ratio):
        parser.error("direct length-ratio bounds are invalid")

    query = read_fasta(args.query)
    target = read_fasta(args.target)
    if not query or not target:
        parser.error("query and target FASTA files must both contain records")
    rows: list[dict[str, object]] = []
    best_rows: list[dict[str, object]] = []
    priority = {"direct_catalogue_recurrence": 2, "related_period_multiple": 1, "no_qualifying_match": 0}
    for query_id, query_sequence in query.items():
        query_rows: list[dict[str, object]] = []
        for target_id, target_sequence in target.items():
            identity, edit_distance, orientation = circular_glocal_identity(query_sequence, target_sequence)
            length_ratio = max(len(query_sequence), len(target_sequence)) / min(len(query_sequence), len(target_sequence))
            state, nearest_multiple, multiple_error = classify_match(
                identity,
                length_ratio,
                args.direct_identity,
                args.direct_minimum_length_ratio,
                args.direct_maximum_length_ratio,
                args.related_identity,
                args.maximum_multiple_error,
            )
            row = {
                "query_id": query_id,
                "query_length_bp": len(query_sequence),
                "target_id": target_id,
                "target_length_bp": len(target_sequence),
                "best_orientation": orientation,
                "edit_distance": edit_distance,
                "glocal_edit_identity": f"{identity:.6f}",
                "length_ratio": f"{length_ratio:.6f}",
                "nearest_integer_multiple": nearest_multiple,
                "multiple_relative_error": f"{multiple_error:.6f}",
                "match_state": state,
                "interpretation_boundary": "same-run nested-depth catalogue comparison; not independent replication or proof of novelty",
            }
            rows.append(row)
            query_rows.append(row)
        best_rows.append(
            max(
                query_rows,
                key=lambda row: (
                    priority[str(row["match_state"])],
                    float(row["glocal_edit_identity"]),
                    -abs(float(row["length_ratio"]) - 1.0),
                    str(row["target_id"]),
                ),
            )
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    all_path = args.outdir / "catalogue_matches.tsv"
    best_path = args.outdir / "best_match_per_query.tsv"
    for path, table in ((all_path, rows), (best_path, best_rows)):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(table)
    counts: dict[str, int] = {}
    for row in best_rows:
        state = str(row["match_state"])
        counts[state] = counts.get(state, 0) + 1
    receipt = {
        "schema_version": 1,
        "query": str(args.query.resolve()),
        "query_sha256": sha256(args.query),
        "target": str(args.target.resolve()),
        "target_sha256": sha256(args.target),
        "query_count": len(query),
        "target_count": len(target),
        "pair_count": len(rows),
        "best_match_state_counts": counts,
        "thresholds": {
            "direct_identity": args.direct_identity,
            "direct_minimum_length_ratio": args.direct_minimum_length_ratio,
            "direct_maximum_length_ratio": args.direct_maximum_length_ratio,
            "related_identity": args.related_identity,
            "maximum_multiple_error": args.maximum_multiple_error,
        },
        "alignment": "complete_shorter_vs_doubled_longer_edlib_HW_both_orientations",
        "edlib_version": version("edlib"),
        "output_sha256": {all_path.name: sha256(all_path), best_path.name: sha256(best_path)},
        "interpretation_boundary": "nested-depth catalogue recurrence only; not independent sample/platform support or novelty proof",
    }
    (args.outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
