"""Output schema validators for TandemX MVP files."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


class ValidationError(ValueError):
    """Raised when a TandemX output file does not match its schema."""


@dataclass(frozen=True)
class ValidationResult:
    path: Path
    record_count: int


TSV_SCHEMAS: dict[str, dict[str, set[str]]] = {
    "candidate_reads.tsv": {
        "required": {
            "read_id",
            "candidate_id",
            "read_start",
            "read_end",
            "strand",
            "period_bp",
            "repeat_span_bp",
            "unit_count",
            "score",
            "low_complexity_flag",
            "confidence",
            "warning",
        },
        "numeric": {"read_start", "read_end", "period_bp", "repeat_span_bp", "unit_count", "score"},
    },
    "families.tsv": {
        "required": {
            "family_id",
            "monomer_id",
            "monomer_length_bp",
            "consensus_md5",
            "gc_fraction",
            "support_read_count",
            "support_span_bp",
            "mean_identity",
            "low_complexity_flag",
            "confidence",
            "warning",
        },
        "numeric": {"monomer_length_bp", "gc_fraction", "support_read_count", "support_span_bp", "mean_identity"},
    },
    "copy_number.tsv": {
        "required": {
            "family_id",
            "monomer_length",
            "diagnostic_kmer_count",
            "median_kmer_depth",
            "haploid_depth",
            "estimated_copy_number",
            "estimated_bp",
            "confidence",
            "warning",
        },
        "numeric": {
            "monomer_length",
            "diagnostic_kmer_count",
            "median_kmer_depth",
            "haploid_depth",
            "estimated_copy_number",
            "estimated_bp",
        },
    },
    "assembly_vs_read_cn.tsv": {
        "required": {
            "family_id",
            "read_estimated_bp",
            "assembly_estimated_bp",
            "assembly_read_ratio",
            "status",
            "confidence",
            "warning",
        },
        "numeric": {"read_estimated_bp", "assembly_estimated_bp", "assembly_read_ratio"},
    },
    "probes.rank.tsv": {
        "required": {
            "probe_id",
            "family_id",
            "sequence_length",
            "gc_content",
            "tm",
            "estimated_copy_number",
            "arrayiness_score",
            "specificity_score",
            "off_target_hits",
            "predicted_regions",
            "probe_score",
            "confidence",
            "warning",
        },
        "numeric": {
            "sequence_length",
            "gc_content",
            "tm",
            "estimated_copy_number",
            "arrayiness_score",
            "specificity_score",
            "off_target_hits",
            "probe_score",
        },
    },
    "in_silico_fish.tsv": {
        "required": {"probe_id", "chrom", "start", "end", "predicted_signal", "confidence", "warning"},
        "numeric": {"start", "end", "predicted_signal"},
    },
}


FASTA_HEADER_PATTERNS = {
    "monomers.fa": re.compile(r"^family_id=[^;]+;monomer_id=[^;]+;length_bp=\d+;confidence=[^;]+$"),
    "probes.fa": re.compile(r"^probe_id=[^;]+;family_id=[^;]+;length_bp=\d+;probe_score=[0-9.]+;confidence=[^;]+$"),
}


def validate_project(project_dir: Path) -> list[ValidationResult]:
    results = []
    for path in sorted(project_dir.rglob("*")):
        if not path.is_file():
            continue
        name = path.name
        if name in TSV_SCHEMAS:
            results.append(validate_tsv(path, TSV_SCHEMAS[name]))
        elif name == "repeat_density.bedgraph":
            results.append(validate_bedgraph(path))
        elif name == "arrays.bed":
            results.append(validate_arrays_bed(path))
        elif name in FASTA_HEADER_PATTERNS:
            results.append(validate_tandemx_fasta(path, FASTA_HEADER_PATTERNS[name]))
    if not results:
        raise ValidationError(f"No recognized TandemX output files found under {project_dir}")
    return results


def validate_tsv(path: Path, schema: dict[str, set[str]]) -> ValidationResult:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValidationError(f"TSV file is empty: {path}")
    header = lines[0].split("\t")
    missing = sorted(schema["required"] - set(header))
    if missing:
        raise ValidationError(f"{path} is missing required field(s): {', '.join(missing)}")
    index = {name: header.index(name) for name in header}
    record_count = 0
    for line_number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValidationError(f"{path} line {line_number} has {len(fields)} fields, expected {len(header)}")
        for name in schema["numeric"]:
            parse_numeric(fields[index[name]], path, line_number, name)
        if "confidence" in index and not fields[index["confidence"]]:
            raise ValidationError(f"{path} line {line_number} has empty confidence")
        if "status" in index and not fields[index["status"]]:
            raise ValidationError(f"{path} line {line_number} has empty status")
        if "warning" in index and len(fields[index["warning"]]) == 0:
            # Empty warning is valid and means no warning for this record.
            pass
        record_count += 1
    if record_count == 0:
        raise ValidationError(f"TSV file has no records: {path}")
    return ValidationResult(path=path, record_count=record_count)


def parse_numeric(value: str, path: Path, line_number: int, field: str) -> None:
    try:
        float(value)
    except ValueError as exc:
        raise ValidationError(f"{path} line {line_number} field {field} is not numeric: {value}") from exc


def validate_bedgraph(path: Path) -> ValidationResult:
    record_count = 0
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 4:
            raise ValidationError(f"{path} line {line_number} must have 4 bedGraph fields")
        start = parse_int(fields[1], path, line_number, "start")
        end = parse_int(fields[2], path, line_number, "end")
        parse_numeric(fields[3], path, line_number, "score")
        validate_half_open(path, line_number, start, end)
        record_count += 1
    if record_count == 0:
        raise ValidationError(f"bedGraph file has no records: {path}")
    return ValidationResult(path=path, record_count=record_count)


def validate_arrays_bed(path: Path) -> ValidationResult:
    record_count = 0
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            raise ValidationError(f"{path} line {line_number} must have at least 8 fields")
        start = parse_int(fields[1], path, line_number, "start")
        end = parse_int(fields[2], path, line_number, "end")
        score = parse_int(fields[4], path, line_number, "score")
        validate_half_open(path, line_number, start, end)
        if not 0 <= score <= 1000:
            raise ValidationError(f"{path} line {line_number} BED score must be 0-1000")
        if fields[5] not in {"+", "-", "."}:
            raise ValidationError(f"{path} line {line_number} strand must be '+', '-', or '.'")
        if not fields[6]:
            raise ValidationError(f"{path} line {line_number} confidence is empty")
        record_count += 1
    if record_count == 0:
        raise ValidationError(f"BED file has no records: {path}")
    return ValidationResult(path=path, record_count=record_count)


def parse_int(value: str, path: Path, line_number: int, field: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValidationError(f"{path} line {line_number} field {field} is not an integer: {value}") from exc


def validate_half_open(path: Path, line_number: int, start: int, end: int) -> None:
    if start < 0:
        raise ValidationError(f"{path} line {line_number} start must be >= 0")
    if end <= start:
        raise ValidationError(f"{path} line {line_number} end must be greater than start")


def validate_tandemx_fasta(path: Path, pattern: re.Pattern[str]) -> ValidationResult:
    record_count = 0
    current_header: str | None = None
    current_has_sequence = False
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None and not current_has_sequence:
                raise ValidationError(f"{path} has a header without sequence: {current_header}")
            header = line[1:]
            if not pattern.match(header):
                raise ValidationError(f"{path} line {line_number} has invalid TandemX FASTA header: {header}")
            record_count += 1
            current_header = header
            current_has_sequence = False
        else:
            if record_count == 0:
                raise ValidationError(f"{path} line {line_number} has sequence before first FASTA header")
            current_has_sequence = True
    if record_count == 0:
        raise ValidationError(f"FASTA file is empty: {path}")
    if current_header is not None and not current_has_sequence:
        raise ValidationError(f"{path} has a header without sequence: {current_header}")
    return ValidationResult(path=path, record_count=record_count)
