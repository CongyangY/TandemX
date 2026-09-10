"""Build conservative, portable family summaries for a TandemX run.

This module is deliberately a reporting adapter: it never re-clusters monomers,
re-estimates abundance, or upgrades architecture evidence.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence
from xml.sax.saxutils import escape as xml_escape


SUMMARY_FIELDS = (
    "family_id", "representative_monomer_id", "monomer_length_bp", "gc_fraction",
    "support_read_count", "candidate_array_count", "estimated_abundance_bp",
    "assembly_representation_bp", "assembly_read_ratio", "abundance_deficit_bp",
    "discovery_confidence", "abundance_confidence", "comparison_confidence", "confidence_scope",
    "confidence", "warning",
)


def _record_value(record: object, name: str, default: object = "") -> object:
    return record.get(name, default) if isinstance(record, Mapping) else getattr(record, name, default)


def _allowed_steps(records: Sequence[object]) -> set[str] | None:
    """Return validated successful stages, or ``None`` for standalone reports."""
    if not records:
        return None
    allowed = set()
    for record in records:
        step = str(_record_value(record, "step", ""))
        notes = str(_record_value(record, "notes", ""))
        if (_record_value(record, "exit_status") == 0
                and bool(_record_value(record, "output_validated", False))
                and not notes.startswith("skipped_missing_")):
            allowed.add(step)
    return allowed


def _read_tsv(path: Path, allowed: bool) -> list[dict[str, str]]:
    if not allowed or not path.is_file():
        return []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_arrays(path: Path, allowed: bool) -> list[dict[str, str]]:
    """Read the documented headerless eight-column candidate-array BED."""
    if not allowed or not path.is_file():
        return []
    names = ("chrom", "start", "end", "family_id", "score", "strand", "confidence", "warning")
    rows = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        values = raw.split("\t")
        if not raw:
            continue
        if len(values) != len(names):
            raise ValueError(f"{path} line {line_number} has {len(values)} fields; expected 8")
        try:
            start, end, score = int(values[1]), int(values[2]), int(values[4])
        except ValueError as exc:
            raise ValueError(f"{path} line {line_number} has non-integer BED coordinates or score") from exc
        if start < 0 or end <= start or not 0 <= score <= 1000:
            raise ValueError(f"{path} line {line_number} has invalid BED interval or score")
        rows.append(dict(zip(names, values)))
    return rows


def _number(value: str | None, field: str = "numeric value") -> float | None:
    if value in (None, "", "NA", "null"):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"Invalid non-finite {field}: {value!r}")
    return number


def _integer(value: str | None, field: str = "integer value") -> int | None:
    number = _number(value, field)
    return int(number) if number is not None else None


def _clean(value: object) -> str:
    return str(value).strip() if value not in (None, "NA") else ""


def _fasta_records(path: Path, allowed: bool) -> dict[str, tuple[str, str]]:
    """Read actual representative records keyed by family ID."""
    if not allowed or not path.is_file():
        return {}
    records: dict[str, tuple[str, str]] = {}
    header = ""
    sequence: list[str] = []
    def save() -> None:
        if not header:
            return
        fields = dict(part.split("=", 1) for part in header.split(";") if "=" in part)
        family_id = fields.get("family_id")
        if family_id:
            if family_id in records:
                raise ValueError(f"Duplicate family_id in {path}: {family_id}")
            records[family_id] = (header, "".join(sequence))
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            save(); header = line[1:]; sequence = []
        else:
            sequence.append(line.strip())
    save()
    return records


def _warnings(*values: str) -> str:
    return ";".join(sorted({part for value in values for part in _clean(value).split(";") if part}))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _index_unique(rows: Sequence[Mapping[str, str]], key: str, path: Path) -> dict[str, Mapping[str, str]]:
    indexed: dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row.get(key, "")
        if value in indexed:
            raise ValueError(f"Duplicate {key} in {path}: {value}")
        indexed[value] = row
    return indexed


def _least_confidence(values: Sequence[str | None]) -> str | None:
    observed = [value for value in values if value]
    if not observed:
        return None
    ranks = {"high": 3, "medium": 2, "low": 1, "unresolved": 0}
    return min(observed, key=lambda value: ranks.get(str(value).lower(), -1))


def _run_metadata(outdir: Path, records: Sequence[object], *, quantification_available: bool) -> tuple[list[dict[str, str]], list[dict[str, object]], list[dict[str, object]]]:
    """Read local run metadata without inferring values from it."""
    warnings: list[dict[str, str]] = []
    files: list[dict[str, object]] = []
    defaults_path = outdir / "automatic_defaults.json"
    defaults: dict[str, object] = {}
    if defaults_path.is_file():
        try:
            defaults = json.loads(defaults_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            warnings.append({"code": "automatic_defaults_unreadable", "message": "Automatic-default metadata could not be read."})
    if defaults_path.is_file():
        source = defaults.get("genome_size_source")
        if source == "assembly_total_length_provisional":
            warnings.append({"code": "assembly_length_proxy", "message": "Genome size used an assembly-length proxy; abundance estimates are provisional."})
        elif source == "required" or defaults.get("genome_size_bp") is None:
            warnings.append({"code": "missing_genome_size", "message": "No genome-size value was available; read-abundance estimates may be absent."})
    elif not quantification_available:
        warnings.append({"code": "quantification_not_available", "message": "No validated read-abundance table was available for this report."})
    for relative in ("automatic_defaults.json", "run_config.yaml", "quantify/run_config.yaml",
                     "source_reuse_manifest.json", "pipeline_summary.tsv", "pipeline_summary.json"):
        path = outdir / relative
        if path.is_file():
            files.append({"path": relative, "sha256": _sha256(path), "bytes": path.stat().st_size})
    commands = [{"step": str(_record_value(record, "step", "")), "command": str(_record_value(record, "command", "")),
                 "exit_status": _record_value(record, "exit_status"), "output_validated": _record_value(record, "output_validated")}
                for record in records]
    return warnings, files, commands


def build_report_data(outdir: Path, records: Sequence[object] = ()) -> dict[str, Any]:
    """Build a JSON-serializable evidence view without interpreting stale outputs."""
    outdir = Path(outdir)
    allowed = _allowed_steps(records)
    use = lambda step: allowed is None or step in allowed
    families_path = outdir / "discover" / "families.tsv"
    family_rows = _read_tsv(families_path, use("discover"))
    family_ids = [row.get("family_id", "") for row in family_rows]
    duplicates = sorted(family_id for family_id, count in Counter(family_ids).items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate family_id values in {families_path}: {', '.join(duplicates)}")
    candidates = _read_tsv(outdir / "discover" / "candidate_reads.tsv", use("discover"))
    membership = _read_tsv(outdir / "discover" / "monomer_membership.tsv", use("discover"))
    hierarchy = _read_tsv(outdir / "discover" / "family_hierarchy.tsv", use("discover"))
    copy_rows = _read_tsv(outdir / "quantify" / "copy_number.tsv", use("quantify"))
    arrays = _read_arrays(outdir / "locate" / "arrays.bed", use("locate"))
    comparison_path = outdir / "compare" / "assembly_vs_read_cn.tsv"
    comparisons = _read_tsv(comparison_path, use("compare"))
    canonical_comparison = bool(comparisons)
    comparison_source_path = comparison_path
    if not comparisons:
        comparison_source_path = outdir / "locate" / "assembly_vs_read_cn.tsv"
        comparisons = _read_tsv(comparison_source_path, use("locate"))
    fasta = _fasta_records(outdir / "discover" / "monomers.fa", use("discover"))
    recovery = _read_tsv(outdir / "recovery" / "recovery_candidates.tsv", use("recovery"))

    copy_by_family = _index_unique(copy_rows, "family_id", outdir / "quantify" / "copy_number.tsv")
    comparison_by_family = _index_unique(comparisons, "family_id", comparison_source_path)
    arrays_measured = use("locate") and (outdir / "locate" / "arrays.bed").is_file()
    arrays_by_family: dict[str, int] = defaultdict(int)
    for row in arrays:
        if row.get("family_id"):
            arrays_by_family[row["family_id"]] += 1
    assigned_candidate_families = {
        row.get("candidate_id", ""): row.get("family_id", "")
        for row in membership if row.get("status") == "assigned" and row.get("family_id") not in ("", "NA")
    }
    candidates_by_family: dict[str, int] = defaultdict(int)
    for row in candidates:
        family_id = assigned_candidate_families.get(row.get("candidate_id", ""))
        if family_id:
            candidates_by_family[family_id] += 1

    results: list[dict[str, Any]] = []
    for row in family_rows:
        family_id = row.get("family_id", "")
        cn, comp = copy_by_family.get(family_id, {}), comparison_by_family.get(family_id, {})
        read_bp = _number(cn.get("estimated_bp"), f"estimated_bp for {family_id}")
        assembly_bp = _number(comp.get("assembly_estimated_bp"), f"assembly_estimated_bp for {family_id}")
        deficit = max(read_bp - assembly_bp, 0.0) if read_bp is not None and assembly_bp is not None else None
        discovery_confidence = row.get("confidence") or None
        abundance_confidence = cn.get("confidence") or None
        comparison_confidence = comp.get("confidence") or None
        comparison_status = comp.get("status") or None
        comparison_warning = comp.get("warning", "")
        reported_ratio = _number(comp.get("assembly_read_ratio"), f"assembly_read_ratio for {family_id}")
        reported_read: float | None = None
        if comp and read_bp is not None:
            reported_read = _number(comp.get("read_estimated_bp"), f"read_estimated_bp for {family_id}")
            if canonical_comparison:
                if reported_read is None or abs(reported_read - read_bp) > 1e-3:
                    raise ValueError(f"Comparison read_estimated_bp conflicts with copy_number estimated_bp for {family_id}")
                if read_bp > 0:
                    expected_ratio = assembly_bp / read_bp if assembly_bp is not None else None
                    if (reported_ratio is None or expected_ratio is None
                            or abs(reported_ratio - expected_ratio) > 5e-4):
                        raise ValueError(f"Comparison assembly_read_ratio conflicts with read/assembly values for {family_id}")
                elif reported_ratio == 0:
                    reported_ratio = None
            elif reported_read is None or abs(reported_read - read_bp) > 1e-3:
                comparison_status = "not_comparable"
                comparison_confidence = None
                reported_ratio = None
                comparison_warning = _warnings(comparison_warning, "missing_validated_comparison", "locate_compatibility_read_estimate_mismatch")
            else:
                comparison_warning = _warnings(comparison_warning, "missing_validated_comparison")
        elif not canonical_comparison and arrays_measured:
            comparison_warning = _warnings(comparison_warning, "missing_validated_comparison")
        elif comp:
            reported_read = _number(comp.get("read_estimated_bp"), f"read_estimated_bp for {family_id}")
        if reported_read == 0 and reported_ratio == 0:
            reported_ratio = None
        confidence_values = [discovery_confidence, abundance_confidence, comparison_confidence]
        confidence_scope = "+".join(name for name, value in zip(
            ("discovery", "abundance", "comparison"), confidence_values) if value) or None
        results.append({
            "family_id": family_id,
            "representative_monomer_id": row.get("monomer_id") or None,
            "monomer_length_bp": _integer(row.get("monomer_length_bp"), f"monomer_length_bp for {family_id}"),
            "gc_fraction": _number(row.get("gc_fraction"), f"gc_fraction for {family_id}"),
            "support_read_count": _integer(row.get("support_read_count"), f"support_read_count for {family_id}"),
            "support_span_bp": _number(row.get("support_span_bp"), f"support_span_bp for {family_id}"),
            "candidate_read_count": candidates_by_family.get(family_id) if membership else None,
            "candidate_array_count": arrays_by_family.get(family_id, 0) if arrays_measured else None,
            "estimated_abundance_bp": read_bp,
            "assembly_representation_bp": assembly_bp,
            "assembly_read_ratio": reported_ratio,
            "abundance_deficit_bp": deficit,
            "comparison_status": comparison_status,
            "discovery_confidence": discovery_confidence,
            "abundance_confidence": abundance_confidence,
            "comparison_confidence": comparison_confidence,
            "confidence_scope": confidence_scope,
            "confidence": _least_confidence(confidence_values),
            "warning": _warnings(row.get("warning", ""), cn.get("warning", ""), comparison_warning),
            "monomer_sequence": fasta.get(family_id, ("", ""))[1] or None,
        })
    results.sort(key=lambda row: row["family_id"])
    high_confidence = sum(row["confidence"] == "high" for row in results)
    possible_underrepresented = [row for row in results if row["comparison_status"] == "possible_collapse"]
    unresolved_count = sum(str(row["confidence"] or "").lower() == "unresolved" for row in results)
    low_confidence_count = sum(str(row["confidence"] or "").lower() == "low" for row in results)
    recommended = sorted(possible_underrepresented, key=lambda row: row["abundance_deficit_bp"] or float("-inf"), reverse=True)[:10]
    recovery_status_counts: dict[str, int] = defaultdict(int)
    for row in recovery:
        if row.get("recovery_status"):
            recovery_status_counts[row["recovery_status"]] += 1
    sources = []
    consumed = (
        ("discover", "discover/families.tsv"), ("discover", "discover/monomers.fa"),
        ("discover", "discover/candidate_reads.tsv"), ("discover", "discover/monomer_membership.tsv"),
        ("discover", "discover/family_hierarchy.tsv"), ("quantify", "quantify/copy_number.tsv"),
        ("locate", "locate/arrays.bed"),
        (("compare" if canonical_comparison else "locate"), str(comparison_source_path.relative_to(outdir))),
    )
    for step, relative in consumed:
        path = outdir / relative
        used = (allowed is None or step in allowed) and path.is_file()
        sources.append({"step": step, "path": relative, "used": used,
                        "sha256": _sha256(path) if used else None, "bytes": path.stat().st_size if used else None})
    recovery_path = outdir / "recovery" / "recovery_candidates.tsv"
    recovery_used = use("recovery") and recovery_path.is_file()
    sources.append({"step": "recovery", "path": "recovery/recovery_candidates.tsv", "used": recovery_used,
                    "sha256": _sha256(recovery_path) if recovery_used else None,
                    "bytes": recovery_path.stat().st_size if recovery_used else None})
    run_warnings, metadata_files, commands = _run_metadata(outdir, records, quantification_available=bool(copy_rows))
    return {"schema_version": 1, "run": {"outdir": str(outdir), "validated_steps": sorted(allowed) if allowed is not None else None},
            "summary": {"family_count": len(results), "high_confidence_family_count": high_confidence,
                        "candidate_read_count": len(candidates), "candidate_array_count": len(arrays),
                        "architecture_edge_count": len(hierarchy), "possible_underrepresented_count": len(possible_underrepresented),
                        "unresolved_family_count": unresolved_count, "low_confidence_family_count": low_confidence_count,
                        "recommended_review_family_ids": [row["family_id"] for row in recommended],
                        "recovery_candidate_count": len(recovery),
                        "recovery_status_counts": dict(sorted(recovery_status_counts.items()))},
            "families": results, "architecture_edges": hierarchy, "sources": sources,
            "membership": membership, "recovery_candidates": recovery, "recommended_families": recommended,
            "data_availability_warnings": run_warnings, "provenance": {"files": metadata_files, "commands": commands}}


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row.get(key) for key in fields})


def _graphml(edges: Sequence[Mapping[str, str]], family_ids: Sequence[str]) -> str:
    def attribute(value: str) -> str:
        return xml_escape(value, {'"': "&quot;"})

    nodes = sorted(set(family_ids) | {node for edge in edges for node in (edge.get("shorter_family_id", ""), edge.get("longer_family_id", "")) if node})
    lines = ['<?xml version="1.0" encoding="UTF-8"?>', '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">', '<key id="edge_type" for="edge" attr.name="edge_type" attr.type="string"/>', '<key id="status" for="edge" attr.name="status" attr.type="string"/>', '<graph edgedefault="directed">']
    lines.extend(f'<node id="{attribute(node)}"/>' for node in nodes)
    for index, edge in enumerate(edges):
        source, target = edge.get("shorter_family_id", ""), edge.get("longer_family_id", "")
        if source and target:
            lines.append(f'<edge id="e{index}" source="{attribute(source)}" target="{attribute(target)}"><data key="edge_type">{xml_escape(edge.get("edge_type", ""))}</data><data key="status">{xml_escape(edge.get("status", ""))}</data></edge>')
    return "\n".join([*lines, "</graph>", "</graphml>", ""])


def write_report_data(outdir: Path, records: Sequence[object] = ()) -> dict[str, Any]:
    """Write family-level derivative data products and return their source data."""
    data = build_report_data(outdir, records)
    root = Path(outdir); family_dir = root / "families"
    _write_tsv(root / "summary.tsv", data["families"], SUMMARY_FIELDS)
    (root / "summary.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _write_tsv(family_dir / "families.tsv", data["families"], SUMMARY_FIELDS)
    fasta_lines = []
    for row in data["families"]:
        if row["monomer_sequence"]:
            fasta_lines.extend([
                f">family_id={row['family_id']};monomer_id={row['representative_monomer_id']};"
                f"length_bp={row['monomer_length_bp']};confidence={row['confidence']}",
                row["monomer_sequence"],
            ])
    (family_dir / "monomers.fa").write_text("\n".join(fasta_lines) + ("\n" if fasta_lines else ""), encoding="utf-8")
    _write_tsv(family_dir / "family_members.tsv", data["membership"], (
        "read_id", "candidate_id", "cluster_id", "family_id", "representative_sha256", "edit_distance_upper_bound", "similarity_lower_bound", "minimum_cluster_identity", "compatible_cluster_count", "alternative_cluster_ids", "status", "warning"))
    hierarchy_fields = tuple(data["architecture_edges"][0]) if data["architecture_edges"] else (
        "hierarchy_edge_id", "shorter_family_id", "longer_family_id", "shorter_length_bp", "longer_length_bp", "nearest_integer_multiple", "length_ratio", "multiple_error", "local_identity", "local_overlap_fraction_shorter", "shared_kmer_fraction", "orientation", "edge_type", "status", "warning")
    _write_tsv(family_dir / "family_hierarchy.tsv", data["architecture_edges"], hierarchy_fields)
    _write_tsv(family_dir / "repeat_architecture.tsv", data["architecture_edges"], hierarchy_fields)
    (family_dir / "family_network.graphml").write_text(_graphml(data["architecture_edges"], [row["family_id"] for row in data["families"]]), encoding="utf-8")
    return data
