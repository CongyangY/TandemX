"""Cross-sample pan-repeat catalogue and abundance/representation matrices."""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from tandemx.discover.clustering import cluster_monomers
from tandemx.discover.mvp import CandidateRepeat, RepeatFamily
from tandemx.quantify.mvp import read_monomer_fasta


MANIFEST_FIELDS = {"sample_id", "monomers", "copy_number", "comparison"}
CONFIDENCE_ORDER = {"low": 0, "medium": 1, "high": 2}


@dataclass(frozen=True)
class CohortSample:
    sample_id: str
    monomers: Path
    copy_number: Path
    comparison: Path | None


@dataclass(frozen=True)
class LocalFamily:
    sample_id: str
    family_id: str
    sequence: str
    candidate_id: str


def _resolve_input(base: Path, value: str, *, optional: bool = False) -> Path | None:
    if optional and value in {"", "NA", "--"}:
        return None
    path = Path(value)
    path = path if path.is_absolute() else base / path
    path = path.resolve()
    if not path.is_file():
        raise ValueError(f"cohort input is not a file: {path}")
    return path


def read_manifest(path: Path) -> list[CohortSample]:
    base = path.resolve().parent
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not MANIFEST_FIELDS.issubset(reader.fieldnames):
            missing = sorted(MANIFEST_FIELDS - set(reader.fieldnames or ()))
            raise ValueError(f"cohort manifest is missing fields: {','.join(missing)}")
        samples = []
        seen = set()
        for line_number, row in enumerate(reader, 2):
            sample_id = row["sample_id"].strip()
            if not sample_id:
                raise ValueError(f"cohort manifest line {line_number} has empty sample_id")
            if sample_id in seen:
                raise ValueError(f"duplicate cohort sample_id: {sample_id}")
            seen.add(sample_id)
            samples.append(CohortSample(
                sample_id,
                _resolve_input(base, row["monomers"].strip()),
                _resolve_input(base, row["copy_number"].strip()),
                _resolve_input(base, row["comparison"].strip(), optional=True),
            ))
    if len(samples) < 2:
        raise ValueError("cohort requires at least two samples")
    return samples


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or ()))
            raise ValueError(f"{path} is missing fields: {','.join(missing)}")
        return list(reader)


def load_local_families(samples: list[CohortSample]) -> tuple[list[LocalFamily], list[CandidateRepeat]]:
    local: list[LocalFamily] = []
    candidates: list[CandidateRepeat] = []
    seen: set[tuple[str, str]] = set()
    for sample in samples:
        for monomer in read_monomer_fasta(sample.monomers):
            key = (sample.sample_id, monomer.family_id)
            if key in seen:
                raise ValueError(f"duplicate local family in monomer catalogue: {sample.sample_id}/{monomer.family_id}")
            seen.add(key)
            candidate_id = f"TXC{len(local) + 1:09d}"
            local.append(LocalFamily(sample.sample_id, monomer.family_id, monomer.sequence, candidate_id))
            candidates.append(CandidateRepeat(
                read_id=sample.sample_id,
                candidate_id=candidate_id,
                sequence=monomer.sequence,
                read_start=0,
                read_end=len(monomer.sequence),
                strand=".",
                period_bp=len(monomer.sequence),
                repeat_span_bp=len(monomer.sequence),
                unit_count=1.0,
                score=1.0,
                low_complexity_flag=False,
                confidence="medium",
                warning="cohort_catalogue_member",
            ))
    if not local:
        raise ValueError("cohort monomer catalogues contain no families")
    return local, candidates


def _float(row: dict[str, str], field: str, path: Path) -> float:
    try:
        return float(row[field])
    except ValueError as error:
        raise ValueError(f"{path} field {field} is not numeric for family {row.get('family_id', 'NA')}") from error


def _index_family_rows(
    rows: list[dict[str, str]],
    path: Path,
    expected: set[str],
    *,
    require_complete: bool,
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        family_id = row["family_id"].strip()
        if not family_id:
            raise ValueError(f"{path} contains an empty family_id")
        if family_id in indexed:
            raise ValueError(f"{path} contains duplicate family_id: {family_id}")
        indexed[family_id] = row
    unknown = sorted(set(indexed) - expected)
    if unknown:
        raise ValueError(
            f"{path} contains families absent from its monomer catalogue: {','.join(unknown)}"
        )
    missing = sorted(expected - set(indexed))
    if require_complete and missing:
        raise ValueError(
            f"{path} lacks quantified catalogue families: {','.join(missing)}"
        )
    return indexed


def _least_confidence(rows: Iterable[dict[str, str]]) -> str:
    observed = [row.get("confidence", "").strip() for row in rows]
    if not observed or any(value not in CONFIDENCE_ORDER for value in observed):
        return "low"
    return min(observed, key=CONFIDENCE_ORDER.__getitem__)


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _pan_mapping(
    local: list[LocalFamily], families: list[RepeatFamily], membership: list[dict]
) -> tuple[dict[tuple[str, str], str], list[dict[str, Any]], dict[str, RepeatFamily]]:
    candidate_to_local = {record.candidate_id: record for record in local}
    family_order = {family.family_id: index for index, family in enumerate(families, 1)}
    source_to_pan = {source: f"TXP{index:06d}" for source, index in family_order.items()}
    mapping: dict[tuple[str, str], str] = {}
    rows = []
    for row in membership:
        record = candidate_to_local[row["candidate_id"]]
        pan_id = source_to_pan.get(row["family_id"], "NA")
        if pan_id == "NA":
            raise ValueError(f"cohort family unexpectedly fell below support: {record.sample_id}/{record.family_id}")
        mapping[(record.sample_id, record.family_id)] = pan_id
        rows.append({
            "sample_id": record.sample_id,
            "local_family_id": record.family_id,
            "pan_family_id": pan_id,
            "monomer_length_bp": len(record.sequence),
            "edit_distance_upper_bound": row["edit_distance_upper_bound"],
            "similarity_lower_bound": row["similarity_lower_bound"],
            "compatible_pan_family_count": row["compatible_cluster_count"],
            "status": row["status"],
            "warning": row["warning"],
        })
    by_pan = {source_to_pan[family.family_id]: family for family in families}
    return mapping, sorted(rows, key=lambda row: (row["sample_id"], row["local_family_id"])), by_pan


def _matrix_rows(
    pan_ids: list[str], samples: list[CohortSample], values: dict[tuple[str, str], str]
) -> list[dict[str, str]]:
    return [
        {"pan_family_id": pan_id, **{
            sample.sample_id: values.get((sample.sample_id, pan_id), "NA") for sample in samples
        }}
        for pan_id in pan_ids
    ]


def build_cohort(manifest: Path, outdir: Path, *, cluster_identity: float = 0.95,
                 backend: str = "python", top_families: int = 30) -> dict[str, Any]:
    if not 0 < cluster_identity <= 1:
        raise ValueError("--cluster-identity must be in (0,1]")
    if backend not in {"python", "rust"}:
        raise ValueError("--backend must be python or rust")
    if top_families < 1:
        raise ValueError("--top-families must be positive")
    samples = read_manifest(manifest)
    local, candidates = load_local_families(samples)
    families, membership = cluster_monomers(candidates, 1, cluster_identity, backend)
    mapping, membership_rows, by_pan = _pan_mapping(local, families, membership)
    pan_ids = sorted(by_pan)
    outdir.mkdir(parents=True, exist_ok=True)

    sample_counts: dict[str, set[str]] = {pan_id: set() for pan_id in pan_ids}
    member_counts = {pan_id: 0 for pan_id in pan_ids}
    for row in membership_rows:
        sample_counts[row["pan_family_id"]].add(row["sample_id"])
        member_counts[row["pan_family_id"]] += 1
    pan_rows = []
    fasta = []
    for pan_id in pan_ids:
        family = by_pan[pan_id]
        sequence_hash = hashlib.sha256(family.monomer_sequence.encode()).hexdigest()
        pan_rows.append({
            "pan_family_id": pan_id,
            "monomer_length_bp": family.monomer_length_bp,
            "sample_count": len(sample_counts[pan_id]),
            "member_family_count": member_counts[pan_id],
            "representative_sha256": sequence_hash,
            "mean_assignment_identity": f"{family.mean_identity:.6f}",
            "confidence": family.confidence,
            "warning": "operational_cross_sample_cluster;uncalibrated_cohort_family",
        })
        fasta.append(
            f">pan_family_id={pan_id};length_bp={len(family.monomer_sequence)};samples={len(sample_counts[pan_id])}\n"
            f"{family.monomer_sequence}\n"
        )
    write_tsv(outdir / "pan_families.tsv", list(pan_rows[0]), pan_rows)
    (outdir / "pan_monomers.fa").write_text("".join(fasta), encoding="utf-8")
    write_tsv(outdir / "family_membership.tsv", list(membership_rows[0]), membership_rows)

    local_lengths = {(record.sample_id, record.family_id): len(record.sequence) for record in local}
    local_ids = {
        sample.sample_id: {
            record.family_id for record in local if record.sample_id == sample.sample_id
        }
        for sample in samples
    }
    abundance_rows = []
    abundance_matrix: dict[tuple[str, str], str] = {}
    abundance_interval_low_matrix: dict[tuple[str, str], str] = {}
    abundance_interval_high_matrix: dict[tuple[str, str], str] = {}
    representation_rows = []
    representation_matrix: dict[tuple[str, str], str] = {}
    input_qc_rows = []
    for sample in samples:
        copy_rows = read_tsv(sample.copy_number, {
            "family_id", "estimated_copy_number", "estimated_bp",
            "copy_number_interval_low", "copy_number_interval_high", "confidence", "warning",
        })
        copy_index = _index_family_rows(
            copy_rows,
            sample.copy_number,
            local_ids[sample.sample_id],
            require_complete=True,
        )
        copy_by_pan: dict[str, list[dict[str, str]]] = {}
        for row in copy_index.values():
            pan_id = mapping.get((sample.sample_id, row["family_id"]))
            if pan_id is None:
                raise ValueError(
                    f"catalogue family lacks a pan-family assignment: "
                    f"{sample.sample_id}/{row['family_id']}"
                )
            copy_by_pan.setdefault(pan_id, []).append(row)
        comparison_by_pan: dict[str, list[dict[str, str]]] = {}
        comparison_index: dict[str, dict[str, str]] = {}
        if sample.comparison is not None:
            comparison_rows = read_tsv(sample.comparison, {
                "family_id", "read_estimated_bp", "assembly_estimated_bp",
                "assembly_read_ratio", "status", "confidence", "warning",
            })
            comparison_index = _index_family_rows(
                comparison_rows,
                sample.comparison,
                local_ids[sample.sample_id],
                require_complete=False,
            )
            for row in comparison_index.values():
                pan_id = mapping.get((sample.sample_id, row["family_id"]))
                if pan_id is None:
                    raise ValueError(
                        f"comparison family lacks a pan-family assignment: "
                        f"{sample.sample_id}/{row['family_id']}"
                    )
                comparison_by_pan.setdefault(pan_id, []).append(row)
        comparison_missing = sorted(local_ids[sample.sample_id] - set(comparison_index))
        input_qc_rows.append({
            "sample_id": sample.sample_id,
            "monomer_family_count": len(local_ids[sample.sample_id]),
            "copy_number_family_count": len(copy_index),
            "copy_number_matched_count": len(copy_index),
            "comparison_provided": str(sample.comparison is not None).lower(),
            "comparison_family_count": len(comparison_index),
            "comparison_matched_count": len(comparison_index),
            "monomers_sha256": hashlib.sha256(sample.monomers.read_bytes()).hexdigest(),
            "copy_number_sha256": hashlib.sha256(sample.copy_number.read_bytes()).hexdigest(),
            "comparison_sha256": (
                "NA" if sample.comparison is None
                else hashlib.sha256(sample.comparison.read_bytes()).hexdigest()
            ),
            "status": "complete",
            "warning": (
                "comparison_not_provided"
                if sample.comparison is None
                else (
                    "comparison_missing_catalogue_families:" + ",".join(comparison_missing)
                    if comparison_missing else ""
                )
            ),
        })
        for pan_id in pan_ids:
            rows = copy_by_pan.get(pan_id, [])
            if rows:
                estimated_bp = sum(_float(row, "estimated_bp", sample.copy_number) for row in rows)
                interval_low = sum(
                    _float(row, "copy_number_interval_low", sample.copy_number)
                    * local_lengths[(sample.sample_id, row["family_id"])] for row in rows
                )
                interval_high = sum(
                    _float(row, "copy_number_interval_high", sample.copy_number)
                    * local_lengths[(sample.sample_id, row["family_id"])] for row in rows
                )
                warnings = sorted({item for row in rows for item in row["warning"].split(";") if item})
                warnings.append("summed_marginal_interval_endpoints_not_joint_ci")
                confidence = _least_confidence(rows)
                status = "quantified"
                value = f"{estimated_bp:.4f}"
            else:
                estimated_bp = interval_low = interval_high = None
                warnings = ["pan_family_not_observed_or_not_quantified"]
                confidence = "low"
                status = "not_observed"
                value = "NA"
            abundance_rows.append({
                "sample_id": sample.sample_id,
                "pan_family_id": pan_id,
                "local_family_count": len(rows),
                "estimated_bp": value,
                "estimated_bp_interval_low": "NA" if interval_low is None else f"{interval_low:.4f}",
                "estimated_bp_interval_high": "NA" if interval_high is None else f"{interval_high:.4f}",
                "status": status,
                "confidence": confidence,
                "warning": ";".join(warnings),
            })
            abundance_matrix[(sample.sample_id, pan_id)] = value
            abundance_interval_low_matrix[(sample.sample_id, pan_id)] = (
                "NA" if interval_low is None else f"{interval_low:.4f}"
            )
            abundance_interval_high_matrix[(sample.sample_id, pan_id)] = (
                "NA" if interval_high is None else f"{interval_high:.4f}"
            )

            comparisons = comparison_by_pan.get(pan_id, [])
            if comparisons:
                read_bp = sum(_float(row, "read_estimated_bp", sample.comparison) for row in comparisons)
                assembly_bp = sum(_float(row, "assembly_estimated_bp", sample.comparison) for row in comparisons)
                ratio = assembly_bp / read_bp if read_bp > 0 else None
                statuses = sorted({row["status"] for row in comparisons})
                rep_status = statuses[0] if len(statuses) == 1 else "mixed_local_status"
                rep_warnings = sorted({item for row in comparisons for item in row["warning"].split(";") if item})
                rep_confidence = _least_confidence(comparisons)
            else:
                read_bp = assembly_bp = ratio = None
                rep_status = "not_evaluated" if sample.comparison is None else "not_observed"
                rep_warnings = ["assembly_read_comparison_unavailable"]
                rep_confidence = "low"
            ratio_value = "NA" if ratio is None else f"{ratio:.6f}"
            representation_rows.append({
                "sample_id": sample.sample_id,
                "pan_family_id": pan_id,
                "read_estimated_bp": "NA" if read_bp is None else f"{read_bp:.4f}",
                "assembly_estimated_bp": "NA" if assembly_bp is None else f"{assembly_bp:.4f}",
                "assembly_read_ratio": ratio_value,
                "status": rep_status,
                "confidence": rep_confidence,
                "warning": ";".join(rep_warnings),
            })
            representation_matrix[(sample.sample_id, pan_id)] = ratio_value

    write_tsv(outdir / "cohort_input_qc.tsv", list(input_qc_rows[0]), input_qc_rows)
    write_tsv(outdir / "sample_family_abundance.tsv", list(abundance_rows[0]), abundance_rows)
    write_tsv(outdir / "abundance_matrix.tsv", ["pan_family_id", *(sample.sample_id for sample in samples)],
              _matrix_rows(pan_ids, samples, abundance_matrix))
    write_tsv(
        outdir / "abundance_interval_low_matrix.tsv",
        ["pan_family_id", *(sample.sample_id for sample in samples)],
        _matrix_rows(pan_ids, samples, abundance_interval_low_matrix),
    )
    write_tsv(
        outdir / "abundance_interval_high_matrix.tsv",
        ["pan_family_id", *(sample.sample_id for sample in samples)],
        _matrix_rows(pan_ids, samples, abundance_interval_high_matrix),
    )
    write_tsv(outdir / "sample_family_representation.tsv", list(representation_rows[0]), representation_rows)
    write_tsv(outdir / "representation_matrix.tsv", ["pan_family_id", *(sample.sample_id for sample in samples)],
              _matrix_rows(pan_ids, samples, representation_matrix))

    from tandemx.visualize.cohort import render_cohort_overview

    figure_outputs = render_cohort_overview(outdir, top_families)

    outputs = [
        "pan_families.tsv", "pan_monomers.fa", "family_membership.tsv", "cohort_input_qc.tsv",
        "sample_family_abundance.tsv", "abundance_matrix.tsv",
        "abundance_interval_low_matrix.tsv", "abundance_interval_high_matrix.tsv",
        "sample_family_representation.tsv", "representation_matrix.tsv",
        *(path.name for path in figure_outputs),
    ]
    receipt = {
        "schema_version": 1,
        "complete": True,
        "sample_count": len(samples),
        "local_family_count": len(local),
        "pan_family_count": len(families),
        "cluster_identity": cluster_identity,
        "backend": backend,
        "top_families": top_families,
        "output_sha256": {
            name: hashlib.sha256((outdir / name).read_bytes()).hexdigest() for name in outputs
        },
        "warning": "pan families are operational sequence clusters; NA is not zero abundance",
    }
    (outdir / "cohort_summary.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt
