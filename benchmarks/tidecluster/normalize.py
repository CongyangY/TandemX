"""Normalize TideCluster GFF3 and score it against explicit assembly truth."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Iterable

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import score_arrays, score_families
from benchmarks.challenge.schema import ArrayRecord, digest_file, write_table


@dataclass(frozen=True)
class TideClusterRecord:
    sequence_id: str
    start: int
    end: int
    family_id: str
    period: int
    consensus_sequence: str
    copy_number: float

    def array(self) -> ArrayRecord:
        return ArrayRecord(
            self.sequence_id, self.start, self.end, self.period,
            self.consensus_sequence, self.family_id,
        )


def parse_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in value.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise ValueError(f"Malformed GFF3 attribute: {field}")
        key, item = field.split("=", 1)
        if not key or key in result:
            raise ValueError(f"Duplicate or empty GFF3 attribute: {key}")
        result[key] = item
    return result


def read_gff(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9 GFF3 fields: {path}:{line_number}")
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"Invalid 1-based GFF3 interval: {path}:{line_number}")
            rows.append({
                "sequence_id": fields[0], "start": start - 1, "end": end,
                "source": fields[1], "type": fields[2],
                "attributes": parse_attributes(fields[8]),
            })
    return rows


def normalize_tidecluster(tidehunter_gff: Path, clustering_gff: Path) -> list[TideClusterRecord]:
    tidehunter: dict[tuple[str, int, int], dict[str, str]] = {}
    for row in read_gff(tidehunter_gff):
        key = (str(row["sequence_id"]), int(row["start"]), int(row["end"]))
        attributes = row["attributes"]
        if key in tidehunter:
            raise ValueError(f"Duplicate TideHunter interval: {key}")
        if not {"consensus_sequence", "consensus_length", "copy_number"} <= set(attributes):
            raise ValueError(f"Missing TideHunter attributes: {key}")
        tidehunter[key] = attributes
    normalized: list[TideClusterRecord] = []
    for row in read_gff(clustering_gff):
        key = (str(row["sequence_id"]), int(row["start"]), int(row["end"]))
        if key not in tidehunter:
            raise ValueError(f"TideCluster interval lacks TideHunter provenance: {key}")
        cluster_attributes = row["attributes"]
        if "Name" not in cluster_attributes:
            raise ValueError(f"TideCluster interval lacks Name: {key}")
        source = tidehunter[key]
        sequence = source["consensus_sequence"].upper()
        period = int(source["consensus_length"])
        if period != len(sequence):
            raise ValueError(f"TideHunter consensus length differs from sequence: {key}")
        normalized.append(TideClusterRecord(
            sequence_id=key[0], start=key[1], end=key[2],
            family_id=cluster_attributes["Name"], period=period,
            consensus_sequence=sequence, copy_number=float(source["copy_number"]),
        ))
    if not normalized:
        raise ValueError("No TideCluster tandem-repeat records found")
    return sorted(normalized, key=lambda x: (x.sequence_id, x.start, x.end, x.family_id))


def read_truth(path: Path, sequences: dict[str, str]) -> list[ArrayRecord]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"chrom", "family_id", "start", "end", "period"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Truth TSV is empty or lacks required columns")
    result = []
    for row in rows:
        family = row["family_id"]
        if family not in sequences:
            raise ValueError(f"Truth family absent from catalogue: {family}")
        result.append(ArrayRecord(
            row["chrom"], int(row["start"]), int(row["end"]), int(row["period"]),
            sequences[family], family,
        ))
    return result


def fasta_lengths(path: Path) -> dict[str, int]:
    return {name: len(sequence) for name, sequence in read_fasta(path).items()}


def json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def evaluate(
    tidehunter_gff: Path,
    clustering_gff: Path,
    assembly: Path,
    truth_tsv: Path,
    catalogue: Path,
    outdir: Path,
) -> dict[str, object]:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f"Choose a new empty output directory: {outdir}")
    records = normalize_tidecluster(tidehunter_gff, clustering_gff)
    catalogue_sequences = read_fasta(catalogue)
    truth = read_truth(truth_tsv, catalogue_sequences)
    metrics, matches = score_arrays(
        [record.array() for record in records], truth, fasta_lengths(assembly), 0.5
    )
    family_metrics, family_rows = score_families(
        [record.consensus_sequence for record in records], catalogue_sequences
    )
    metrics.update(family_metrics)
    metrics.update({
        "normalized_coordinate_system": "0-based_half-open",
        "source_coordinate_system": "GFF3_1-based_inclusive",
        "warning": "known_planted_assembly_truth;smoke_test_not_publication_scale",
    })
    outdir.mkdir(parents=True, exist_ok=True)
    normalized_rows = [{
        "sequence_id": record.sequence_id, "start": record.start, "end": record.end,
        "family_id": record.family_id, "period": record.period,
        "consensus_sequence": record.consensus_sequence, "copy_number": record.copy_number,
        "source": "TideCluster_clustering_joined_to_TideHunter",
        "warning": "period_and_consensus_from_exact_interval_TideHunter_record",
    } for record in records]
    write_table(outdir / "normalized_arrays.tsv", normalized_rows, list(normalized_rows[0]))
    write_table(outdir / "matches.tsv", matches, list(matches[0]))
    write_table(outdir / "family_recovery.tsv", family_rows, list(family_rows[0]))
    (outdir / "metrics.json").write_text(json.dumps(json_safe(metrics), indent=2) + "\n")
    receipt = {
        "complete": True,
        "inputs": {str(path): {"sha256": digest_file(path), "bytes": path.stat().st_size}
                   for path in (tidehunter_gff, clustering_gff, assembly, truth_tsv, catalogue)},
        "outputs": {name: digest_file(outdir / name) for name in
                    ("normalized_arrays.tsv", "matches.tsv", "family_recovery.tsv", "metrics.json")},
    }
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tidehunter-gff", required=True, type=Path)
    parser.add_argument("--clustering-gff", required=True, type=Path)
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        evaluate(args.tidehunter_gff, args.clustering_gff, args.assembly, args.truth, args.catalogue, args.outdir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
