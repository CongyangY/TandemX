"""Output manifests and human-readable reports for TandemX pipeline runs."""

from __future__ import annotations

import csv
import shlex
from pathlib import Path
from typing import Protocol, Sequence


class ReportConfig(Protocol):
    reads: Path
    assembly: Path | None
    outdir: Path
    steps: tuple[str, ...]


class ReportStep(Protocol):
    step: str
    runtime_seconds: float
    exit_status: int
    output_validated: bool
    notes: str


MANIFEST_FIELDS = (
    "step",
    "output_type",
    "file_path",
    "exists",
    "file_size_bytes",
    "description",
    "required_for_next_step",
    "notes",
)
OUTPUT_SPECS = {
    "discover": (
        ("candidate_reads", "discover/candidate_reads.tsv", "Read-local tandem repeat candidates.", "family review"),
        ("monomer_catalog", "discover/monomers.fa", "Discovered representative monomer sequences.", "quantify,locate,probe,visualize"),
        ("family_catalog", "discover/families.tsv", "Discovered repeat family summary.", "interpretation"),
    ),
    "quantify": (
        ("copy_number", "quantify/copy_number.tsv", "Read-based repeat copy-number estimates.", "locate,probe,visualize"),
    ),
    "locate": (
        ("repeat_density", "locate/repeat_density.bedgraph", "Assembly repeat-density track.", "assembly interpretation"),
        ("arrays", "locate/arrays.bed", "Candidate assembly repeat arrays.", "probe"),
        ("assembly_read_comparison", "locate/assembly_vs_read_cn.tsv", "Assembly-versus-read abundance comparison.", "visualize"),
    ),
    "probe": (
        ("probe_fasta", "probe/probes.fa", "Ranked probe sequences.", "probe synthesis review"),
        ("probe_ranking", "probe/probes.rank.tsv", "Probe scores and specificity fields.", "visualize"),
        ("in_silico_fish", "probe/in_silico_fish.tsv", "Predicted probe signal regions.", "visualize"),
    ),
    "visualize": (
        ("catalogue_svg", "visualize/catalogue_summary.svg", "Editable repeat catalogue summary.", "publication review"),
        ("catalogue_pdf", "visualize/catalogue_summary.pdf", "PDF repeat catalogue summary.", "publication review"),
        ("assembly_read_svg", "visualize/assembly_vs_read.svg", "Editable assembly-versus-read plot.", "publication review"),
        ("assembly_read_pdf", "visualize/assembly_vs_read.pdf", "PDF assembly-versus-read plot.", "publication review"),
        ("fish_svg", "visualize/in_silico_fish.svg", "Editable predicted FISH signal plot.", "publication review"),
        ("fish_pdf", "visualize/in_silico_fish.pdf", "PDF predicted FISH signal plot.", "publication review"),
    ),
}
PIPELINE_OUTPUT_SPECS = (
    ("pipeline_summary_tsv", "pipeline_summary.tsv", "Step-level pipeline summary for tabular analysis."),
    ("pipeline_summary_json", "pipeline_summary.json", "Step-level pipeline summary for programmatic use."),
    ("pipeline_log", "pipeline.log", "Pipeline step completion log."),
    ("run_report", "run_report.md", "Human-readable run report."),
)


def count_data_rows(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open("rt", encoding="utf-8") as handle:
        return max(0, sum(1 for line in handle if line.strip()) - 1)


def collect_output_warnings(config: ReportConfig) -> list[str]:
    warnings: set[str] = set()
    for relative_path in (
        "discover/candidate_reads.tsv",
        "discover/families.tsv",
        "quantify/copy_number.tsv",
        "locate/assembly_vs_read_cn.tsv",
        "probe/probes.rank.tsv",
        "probe/in_silico_fish.tsv",
    ):
        path = config.outdir / relative_path
        if not path.is_file():
            continue
        with path.open("rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                warning = row.get("warning", "").strip()
                if warning:
                    warnings.update(part for part in warning.split(";") if part)
    return sorted(warnings)


def write_run_report(config: ReportConfig, records: Sequence[ReportStep]) -> None:
    record_by_step = {record.step: record for record in records}
    completed = [
        record.step
        for record in records
        if record.exit_status == 0 and not record.notes.startswith("skipped_missing_")
    ]
    skipped = [
        f"{record.step}: {record.notes}"
        for record in records
        if record.notes.startswith("skipped_")
    ]
    failures = [
        f"{record.step}: {record.notes or f'exit_status={record.exit_status}'}"
        for record in records
        if record.exit_status != 0
    ]
    validation = record_by_step.get("validate")
    if validation is None:
        validation_status = "not requested"
    elif validation.exit_status == 0 and validation.output_validated:
        validation_status = "passed"
    else:
        validation_status = f"failed (exit status {validation.exit_status})"
    notes = skipped + failures + [
        f"output warning: {warning}" for warning in collect_output_warnings(config)
    ]
    total_runtime = sum(record.runtime_seconds for record in records)
    family_path = config.outdir / "discover" / "families.tsv"
    candidate_path = config.outdir / "discover" / "candidate_reads.tsv"
    copy_number_path = config.outdir / "quantify" / "copy_number.tsv"

    lines = [
        "# TandemX Run Report",
        "",
        "## Run overview",
        "",
        f"- Input reads: `{config.reads}`",
        f"- Input assembly: `{config.assembly}`" if config.assembly else "- Input assembly: not provided",
        f"- Steps requested: {', '.join(config.steps)}",
        f"- Steps completed: {', '.join(completed) if completed else 'none'}",
        f"- Total measured step runtime: {total_runtime:.3f} seconds",
        f"- Validation status: {validation_status}",
        "",
        "## Step runtimes",
        "",
        "| Step | Runtime (s) | Exit status | Output validated | Notes |",
        "|---|---:|---:|---|---|",
    ]
    lines.extend(
        f"| {record.step} | {record.runtime_seconds:.3f} | {record.exit_status} | "
        f"{str(record.output_validated).lower()} | {record.notes or ''} |"
        for record in records
    )
    lines.extend(
        [
            "",
            "## Output counts",
            "",
            f"- Discovered families: {count_data_rows(family_path)}",
            f"- Candidate reads: {count_data_rows(candidate_path)}",
            f"- Copy-number rows: {count_data_rows(copy_number_path)}",
            "",
            "## Main outputs",
            "",
            f"- Family catalogue: `{family_path}`",
            f"- Monomer FASTA: `{config.outdir / 'discover' / 'monomers.fa'}`",
            f"- Copy-number table: `{copy_number_path}`" if copy_number_path.is_file() else "- Copy-number table: not generated",
            f"- Output manifest: `{config.outdir / 'output_manifest.tsv'}`",
            f"- Pipeline summary: `{config.outdir / 'pipeline_summary.tsv'}`",
            "",
            "## Warnings and skipped steps",
            "",
        ]
    )
    lines.extend(f"- {note}" for note in notes)
    if not notes:
        lines.append("- None recorded.")
    lines.extend(
        [
            "",
            "## Suggested next commands",
            "",
            "```bash",
            f"tandemx validate --project {shlex.quote(str(config.outdir))}",
        ]
    )
    if (config.outdir / "discover" / "monomers.fa").is_file():
        lines.extend(
            [
                "python benchmarks/scripts/check_known_repeats_against_catalog.py \\",
                f"  --catalog {shlex.quote(str(config.outdir / 'discover' / 'monomers.fa'))} \\",
                "  --known known_repeats.fa \\",
                f"  --out {shlex.quote(str(config.outdir / 'known_repeat_matches.tsv'))}",
            ]
        )
    lines.extend(["```", ""])
    (config.outdir / "run_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_output_manifest(config: ReportConfig, records: Sequence[ReportStep]) -> None:
    record_by_step = {record.step: record for record in records}
    rows: list[dict[str, object]] = []
    for step in config.steps:
        if step not in OUTPUT_SPECS:
            continue
        record = record_by_step.get(step)
        notes = record.notes if record is not None else "step_not_recorded"
        for output_type, relative_path, description, required_for_next_step in OUTPUT_SPECS[step]:
            path = config.outdir / relative_path
            exists = path.is_file()
            row_notes = notes
            if not exists and not row_notes:
                row_notes = "output_missing"
            rows.append(
                {
                    "step": step,
                    "output_type": output_type,
                    "file_path": str(path),
                    "exists": str(exists).lower(),
                    "file_size_bytes": path.stat().st_size if exists else 0,
                    "description": description,
                    "required_for_next_step": required_for_next_step,
                    "notes": row_notes,
                }
            )
    for output_type, filename, description in PIPELINE_OUTPUT_SPECS:
        path = config.outdir / filename
        exists = path.is_file()
        rows.append(
            {
                "step": "pipeline",
                "output_type": output_type,
                "file_path": str(path),
                "exists": str(exists).lower(),
                "file_size_bytes": path.stat().st_size if exists else 0,
                "description": description,
                "required_for_next_step": "run review",
                "notes": "" if exists else "output_missing",
            }
        )
    with (config.outdir / "output_manifest.tsv").open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
