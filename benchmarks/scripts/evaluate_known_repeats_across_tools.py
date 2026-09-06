"""Score source-backed repeat queries against actual tool consensus outputs.

This is a post hoc sequence-recovery endpoint. It does not infer that a query is
present in the tested donor, and it does not turn unmatched predictions into
false positives.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import canonical_monomer
from benchmarks.challenge.schema import digest_file, iter_table, write_table
from benchmarks.challenge.sequence_metrics import (
    cyclic_edit_similarity,
    cyclic_reaches_threshold,
    score_threshold_recovery,
)


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("prediction must be TOOL=PATH")
    tool, raw_path = value.split("=", 1)
    if not tool or not tool.replace("-", "").replace("_", "").isalnum() or not raw_path:
        raise argparse.ArgumentTypeError("prediction must use a simple nonempty TOOL and PATH")
    return tool, Path(raw_path)


def table_sequences(path: Path) -> tuple[list[tuple[str, str]], int]:
    sequences: list[tuple[str, str]] = []
    rows = 0
    for row in iter_table(path, {"sequence"}):
        rows += 1
        sequence = row["sequence"].upper()
        if sequence and set(sequence) - set("ACGTN"):
            raise ValueError(f"Non-ACGTN consensus in {path}")
        if sequence:
            if row.get("read_id"):
                coordinates = (
                    f":{row['start']}-{row['end']}"
                    if row.get("start") not in (None, "") and row.get("end") not in (None, "")
                    else ""
                )
                label = f"{row['read_id']}{coordinates}"
            else:
                label = f"row_{rows}"
            sequences.append((label, sequence))
    return sequences, rows


def compatible_with_any_query(sequence: str, queries: dict[str, str], threshold: float) -> bool:
    return any(
        abs(len(sequence) - len(query))
        <= math.floor((1 - threshold + 1e-12) * max(len(sequence), len(query)))
        for query in queries.values()
    )


def evaluate(
    known_fasta: Path,
    query_ids: list[str],
    catalogs: dict[str, Path],
    consensus_tables: dict[str, Path],
    outdir: Path,
    material: str,
    evidence_boundary: str,
    threshold: float = 0.9,
) -> dict:
    if outdir.exists():
        raise ValueError(f"Output directory already exists: {outdir}")
    if not material.strip() or not evidence_boundary.strip() or not 0 < threshold <= 1:
        raise ValueError("Require material, evidence boundary and threshold in (0,1]")
    if set(catalogs) & set(consensus_tables) or not catalogs and not consensus_tables:
        raise ValueError("Prediction tool names must be unique and nonempty")
    known = read_fasta(known_fasta)
    if len(query_ids) != len(set(query_ids)) or not query_ids:
        raise ValueError("Require unique selected query IDs")
    missing = set(query_ids) - set(known)
    if missing:
        raise ValueError(f"Unknown selected query IDs: {sorted(missing)}")
    queries = {query_id: known[query_id] for query_id in query_ids}

    inputs: dict[str, dict] = {}
    loaded: dict[str, tuple[str, list[tuple[str, str]], int, Path]] = {}
    for tool, path in sorted(catalogs.items()):
        if not path.is_file():
            raise ValueError(f"Missing catalog: {path}")
        records = read_fasta(path)
        loaded[tool] = ("family_catalog", list(records.items()), len(records), path)
    for tool, path in sorted(consensus_tables.items()):
        if not path.is_file():
            raise ValueError(f"Missing consensus table: {path}")
        sequences, rows = table_sequences(path)
        loaded[tool] = ("array_consensus", sequences, rows, path)

    summary_rows: list[dict] = []
    detail_rows: list[dict] = []
    for tool, (unit, raw_sequences, source_rows, path) in sorted(loaded.items()):
        sources: dict[str, dict[str, object]] = {}
        for label, sequence in raw_sequences:
            canonical = canonical_monomer(sequence)
            entry = sources.setdefault(canonical, {"count": 0, "examples": []})
            entry["count"] = int(entry["count"]) + 1
            examples = entry["examples"]
            if isinstance(examples, list) and len(examples) < 5:
                examples.append(label)
        distinct = sorted(sources)
        compatible = [sequence for sequence in distinct if compatible_with_any_query(sequence, queries, threshold)]
        metrics, details = score_threshold_recovery(compatible, queries, threshold)
        summary_rows.append(
            dict(
                material=material,
                tool=tool,
                prediction_unit=unit,
                source_rows=source_rows,
                nonempty_consensus_rows=len(raw_sequences),
                distinct_canonical_consensus_count=len(distinct),
                length_compatible_consensus_count=len(compatible),
                selected_query_count=len(queries),
                recovered_query_count=metrics["recovered_family_count"],
                selected_query_recall=metrics["cyclic_monomer_recall"],
                homologous_consensus_fraction=metrics["homologous_consensus_fraction"],
                threshold=threshold,
                warning=(
                    "posthoc_source_query_recovery_not_donor_presence_or_complete_family_recall;"
                    "homologous_consensus_fraction_not_precision;prediction_units_differ"
                ),
            )
        )
        for detail in details:
            query = queries[detail["truth_id"]]
            supported = [
                sequence for sequence in compatible
                if cyclic_reaches_threshold(query, sequence, threshold)
            ]
            best_score, best = max(
                ((cyclic_edit_similarity(query, sequence), sequence) for sequence in supported),
                default=(None, None),
                key=lambda item: (item[0], item[1]) if item[0] is not None else (-1.0, ""),
            )
            assigned = detail["assigned_sequence_index"]
            matched = compatible[assigned] if isinstance(assigned, int) else None
            detail_rows.append(
                dict(
                    material=material,
                    tool=tool,
                    prediction_unit=unit,
                    query_id=detail["truth_id"],
                    query_length=len(queries[detail["truth_id"]]),
                    recovered=detail["recovered"],
                    matched_consensus_length=len(matched) if matched else None,
                    matched_cyclic_edit_similarity=(
                        cyclic_edit_similarity(queries[detail["truth_id"]], matched) if matched else None
                    ),
                    matched_consensus_sha256=(
                        hashlib.sha256(matched.encode()).hexdigest() if matched else None
                    ),
                    matched_source_occurrence_count=(int(sources[matched]["count"]) if matched else None),
                    matched_source_examples=(
                        "|".join(str(value) for value in sources[matched]["examples"]) if matched else None
                    ),
                    supported_distinct_consensus_count=len(supported),
                    best_supported_cyclic_edit_similarity=best_score,
                    best_supported_consensus_length=len(best) if best else None,
                    best_supported_consensus_sha256=(
                        hashlib.sha256(best.encode()).hexdigest() if best else None
                    ),
                    best_supported_source_occurrence_count=(int(sources[best]["count"]) if best else None),
                    best_supported_source_examples=(
                        "|".join(str(value) for value in sources[best]["examples"]) if best else None
                    ),
                    threshold=threshold,
                    evidence_boundary=evidence_boundary,
                )
            )
        inputs[tool] = {
            "path": str(path.resolve()),
            "sha256": digest_file(path),
            "prediction_unit": unit,
            "source_rows": source_rows,
        }

    outdir.mkdir(parents=True)
    summary_path = outdir / "known_query_summary.tsv"
    detail_path = outdir / "known_query_details.tsv"
    write_table(summary_path, summary_rows, list(summary_rows[0]))
    write_table(detail_path, detail_rows, list(detail_rows[0]))
    receipt = {
        "complete": True,
        "material": material,
        "selected_query_ids": query_ids,
        "threshold": threshold,
        "known_fasta": str(known_fasta.resolve()),
        "known_fasta_sha256": digest_file(known_fasta),
        "prediction_inputs": inputs,
        "script_sha256": digest_file(Path(__file__)),
        "outputs": {
            summary_path.name: digest_file(summary_path),
            detail_path.name: digest_file(detail_path),
        },
        "evidence_boundary": evidence_boundary,
        "interpretation": (
            "selected source-query sequence recovery only; not verified donor presence, "
            "genome-wide family recall, family precision or array-coordinate accuracy"
        ),
    }
    (outdir / "evaluation_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--known", required=True, type=Path)
    parser.add_argument("--query-id", required=True, action="append")
    parser.add_argument("--catalog", action="append", default=[], type=parse_named_path,
                        help="Family catalog as TOOL=FASTA; repeat for tools.")
    parser.add_argument("--consensus-table", action="append", default=[], type=parse_named_path,
                        help="Normalized array consensuses as TOOL=TSV; repeat for tools.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--material", required=True)
    parser.add_argument("--evidence-boundary", required=True)
    parser.add_argument("--threshold", type=float, default=0.9)
    args = parser.parse_args()
    catalogs, tables = dict(args.catalog), dict(args.consensus_table)
    if len(catalogs) != len(args.catalog) or len(tables) != len(args.consensus_table):
        parser.error("Duplicate tool names are not allowed")
    try:
        evaluate(args.known, args.query_id, catalogs, tables, args.outdir,
                 args.material, args.evidence_boundary, args.threshold)
    except ValueError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
