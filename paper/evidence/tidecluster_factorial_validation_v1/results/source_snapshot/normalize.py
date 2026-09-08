"""Normalize TideCluster GFF3 and score it against explicit assembly truth."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
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


@dataclass(frozen=True)
class ResolvedTideClusterRecord:
    sequence_id: str
    start: int
    end: int
    family_id: str
    period: int
    consensus_sequence: str
    copy_number: float | None
    representative_tidehunter_id: str
    representative_selection_source: str
    supporting_tidehunter_count: int
    copy_number_source: str

    def array(self) -> ArrayRecord:
        return ArrayRecord(
            self.sequence_id,
            self.start,
            self.end,
            self.period,
            self.consensus_sequence,
            self.family_id,
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


def normalize_resolved_tidecluster(
    tidehunter_gff: Path,
    intermediate_clustering_gff: Path,
    clustering_gff: Path,
    family_consensus_fasta: Path | None = None,
    *,
    allow_empty: bool = False,
) -> list[ResolvedTideClusterRecord]:
    """Join merged TideCluster intervals through its representative-ID map.

    TideCluster may merge adjacent or overlapping TideHunter hits before writing
    its final GFF.  In that case no single TideHunter interval has the final
    coordinates.  The intermediate clustering GFF retains the representative
    TideHunter ID for each final interval and is therefore the authoritative
    bridge to period and consensus provenance.
    """
    tidehunter_by_id: dict[str, dict[str, str]] = {}
    tidehunter_by_interval: dict[tuple[str, int, int], list[dict[str, str]]] = {}
    for row in read_gff(tidehunter_gff):
        attributes = row["attributes"]
        required = {"ID", "consensus_sequence", "consensus_length", "copy_number"}
        if not required <= set(attributes):
            raise ValueError("TideHunter row lacks resolved-normalization attributes")
        identifier = attributes["ID"]
        if identifier in tidehunter_by_id:
            raise ValueError(f"Duplicate TideHunter ID: {identifier}")
        tidehunter_by_id[identifier] = attributes
        key = (str(row["sequence_id"]), int(row["start"]), int(row["end"]))
        tidehunter_by_interval.setdefault(key, []).append(attributes)

    family_by_representative: dict[str, str] = {}
    if family_consensus_fasta is not None:
        for line_number, line in enumerate(family_consensus_fasta.open(encoding="utf-8"), 1):
            if not line.startswith(">"):
                continue
            match = re.fullmatch(r">(TRC_\d+)_(rep\S+)", line.rstrip())
            if match is None:
                raise ValueError(
                    f"Malformed TideCluster consensus header: "
                    f"{family_consensus_fasta}:{line_number}"
                )
            family_id, representative_id = match.groups()
            if representative_id in family_by_representative:
                raise ValueError(
                    f"Duplicate TideCluster consensus representative: {representative_id}"
                )
            family_by_representative[representative_id] = family_id

    representative_by_interval: dict[tuple[str, int, int], str] = {}
    intermediate_by_sequence: dict[str, list[tuple[int, int, str]]] = {}
    for row in read_gff(intermediate_clustering_gff):
        key = (str(row["sequence_id"]), int(row["start"]), int(row["end"]))
        attributes = row["attributes"]
        if "Name" not in attributes or key in representative_by_interval:
            raise ValueError(f"Invalid intermediate TideCluster interval: {key}")
        representative_by_interval[key] = attributes["Name"]
        intermediate_by_sequence.setdefault(key[0], []).append(
            (key[1], key[2], attributes["Name"])
        )

    missing_representatives = sorted(
        set(representative_by_interval.values()) - set(tidehunter_by_id)
    )
    if missing_representatives:
        raise ValueError(
            "TideCluster representatives are absent from TideHunter output: "
            + ",".join(missing_representatives[:5])
        )
    if family_consensus_fasta is not None:
        missing_membership = sorted(
            set(representative_by_interval.values()) - set(family_by_representative)
        )
        if missing_membership:
            raise ValueError(
                "TideCluster representatives lack family membership: "
                + ",".join(missing_membership[:5])
            )

    normalized: list[ResolvedTideClusterRecord] = []
    for row in read_gff(clustering_gff):
        key = (str(row["sequence_id"]), int(row["start"]), int(row["end"]))
        attributes = row["attributes"]
        if "Name" not in attributes:
            raise ValueError(f"Final TideCluster interval lacks family name: {key}")
        family_id = attributes["Name"]
        if key in representative_by_interval:
            representative_id = representative_by_interval[key]
            representative_selection_source = "exact_intermediate_interval"
            supporting_tidehunter_count = 1
            if (
                family_consensus_fasta is not None
                and family_by_representative[representative_id] != family_id
            ):
                raise ValueError(
                    f"Exact TideCluster representative-family mismatch: {key}"
                )
        else:
            if family_consensus_fasta is None:
                raise ValueError(
                    f"Final TideCluster interval lacks exact intermediate provenance: {key}; "
                    "family consensus membership is required to resolve clipping or merging"
                )
            candidates = []
            for start, end, candidate_id in intermediate_by_sequence.get(key[0], []):
                overlap = min(key[2], end) - max(key[1], start)
                if overlap > 0 and family_by_representative[candidate_id] == family_id:
                    candidates.append((overlap, start, end, candidate_id))
            if not candidates:
                raise ValueError(
                    f"Final TideCluster interval lacks family-consistent overlap provenance: {key}"
                )
            clipped = sorted(
                (max(key[1], start), min(key[2], end))
                for _overlap, start, end, _candidate_id in candidates
            )
            cursor = key[1]
            for start, end in clipped:
                if start > cursor:
                    break
                cursor = max(cursor, end)
            if cursor < key[2]:
                raise ValueError(
                    f"Family-consistent TideCluster provenance does not cover final interval: {key}"
                )
            selected = min(
                candidates,
                key=lambda item: (-item[0], item[1], item[2], item[3]),
            )
            representative_id = selected[3]
            supporting_tidehunter_count = len(candidates)
            representative_selection_source = (
                "unique_family_consistent_overlap_after_clipping"
                if len(candidates) == 1
                else "longest_family_consistent_overlap_after_merge"
            )
        representative = tidehunter_by_id[representative_id]
        sequence = representative["consensus_sequence"].upper()
        period = int(representative["consensus_length"])
        if period != len(sequence):
            raise ValueError(
                f"TideHunter representative length differs from sequence: {representative_id}"
            )
        exact = tidehunter_by_interval.get(key, [])
        if len(exact) == 1:
            copy_number = float(exact[0]["copy_number"])
            copy_number_source = "exact_tidehunter_interval"
        else:
            copy_number = None
            copy_number_source = "unavailable_after_interval_merge_or_resolution"
        normalized.append(
            ResolvedTideClusterRecord(
                sequence_id=key[0],
                start=key[1],
                end=key[2],
                family_id=family_id,
                period=period,
                consensus_sequence=sequence,
                copy_number=copy_number,
                representative_tidehunter_id=representative_id,
                representative_selection_source=representative_selection_source,
                supporting_tidehunter_count=supporting_tidehunter_count,
                copy_number_source=copy_number_source,
            )
        )
    if not normalized and not allow_empty:
        raise ValueError("No resolved TideCluster tandem-repeat records found")
    return sorted(
        normalized, key=lambda item: (item.sequence_id, item.start, item.end, item.family_id)
    )


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
