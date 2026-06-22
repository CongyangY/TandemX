"""Compare TandemX run directories for parameter and result consistency."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from tandemx.io.sequences import read_sequence_records


RESULT_AFFECTING_DISCOVER_PARAMETERS = (
    "reads",
    "max_reads",
    "max_read_bases",
    "sample_rate",
    "seed",
    "min_read_length",
    "min_monomer_len",
    "max_monomer_len",
    "min_support_reads",
    "min_repeat_span",
    "kmer_size",
    "top_periods",
    "min_seed_occurrences",
    "min_spacing_support",
    "max_pairs_per_kmer",
    "kmer_backend",
)


COMPARE_FIELDS = (
    "category",
    "item",
    "run_a_value",
    "run_b_value",
    "same",
    "direct_comparison_impact",
    "notes",
)


@dataclass(frozen=True)
class RunArtifacts:
    """Normalized paths and metrics for a TandemX discover or pipeline run."""

    label: str
    root_dir: Path
    discover_dir: Path
    run_config_path: Path
    pipeline_summary_path: Path | None
    families_path: Path
    candidate_reads_path: Path
    monomers_path: Path
    config: dict[str, Any]
    pipeline_rows: list[dict[str, str]]
    family_rows: list[dict[str, str]]
    candidate_count: int
    monomer_lengths: list[int]


def compare_run_directories(run_a: Path, run_b: Path, outdir: Path) -> list[dict[str, str]]:
    """Compare two TandemX run directories and write TSV/Markdown reports."""
    artifacts_a = load_run_artifacts(run_a, "run_a")
    artifacts_b = load_run_artifacts(run_b, "run_b")
    rows = build_comparison_rows(artifacts_a, artifacts_b)
    outdir.mkdir(parents=True, exist_ok=True)
    write_compare_tsv(outdir / "compare_runs.tsv", rows)
    write_compare_markdown(outdir / "compare_runs.md", artifacts_a, artifacts_b, rows)
    return rows


def load_run_artifacts(run_dir: Path, label: str) -> RunArtifacts:
    """Load discover outputs from either a pipeline root or a standalone discover output."""
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")
    if not run_dir.is_dir():
        raise NotADirectoryError(f"Run path is not a directory: {run_dir}")

    discover_dir = find_discover_dir(run_dir)
    run_config_path = discover_dir / "run_config.yaml"
    families_path = discover_dir / "families.tsv"
    candidate_reads_path = discover_dir / "candidate_reads.tsv"
    monomers_path = discover_dir / "monomers.fa"
    for path in (run_config_path, families_path, candidate_reads_path, monomers_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required TandemX discover output is missing: {path}")

    pipeline_summary_path = run_dir / "pipeline_summary.tsv"
    if not pipeline_summary_path.is_file():
        pipeline_summary_path = None

    family_rows = read_tsv(families_path)
    return RunArtifacts(
        label=label,
        root_dir=run_dir,
        discover_dir=discover_dir,
        run_config_path=run_config_path,
        pipeline_summary_path=pipeline_summary_path,
        families_path=families_path,
        candidate_reads_path=candidate_reads_path,
        monomers_path=monomers_path,
        config=parse_run_config(run_config_path),
        pipeline_rows=read_tsv(pipeline_summary_path) if pipeline_summary_path else [],
        family_rows=family_rows,
        candidate_count=count_tsv_records(candidate_reads_path),
        monomer_lengths=read_monomer_lengths(monomers_path, family_rows),
    )


def find_discover_dir(run_dir: Path) -> Path:
    """Return the directory containing discover run_config and catalog outputs."""
    candidates: list[Path] = []
    for candidate in (run_dir, run_dir / "discover"):
        if (candidate / "run_config.yaml").is_file() and (candidate / "families.tsv").is_file():
            candidates.append(candidate)
    candidates.extend(
        sorted(
            child
            for child in run_dir.iterdir()
            if child.is_dir()
            and child.name.startswith("discover")
            and (child / "run_config.yaml").is_file()
            and (child / "families.tsv").is_file()
        )
    )
    unique: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved not in seen:
            unique.append(candidate)
            seen.add(resolved)
    if len(unique) == 1:
        return unique[0]
    if not unique:
        raise FileNotFoundError(
            f"Could not find a discover output under {run_dir}. "
            "Expected run_config.yaml and families.tsv in the root, discover/, or discover* subdirectory."
        )
    joined = ", ".join(str(path) for path in unique)
    raise ValueError(f"Found multiple discover outputs under {run_dir}; choose one explicitly: {joined}")


def build_comparison_rows(a: RunArtifacts, b: RunArtifacts) -> list[dict[str, str]]:
    """Build row-wise parameter and result comparison records."""
    rows: list[dict[str, str]] = []
    params_a = a.config.get("parameters", {})
    params_b = b.config.get("parameters", {})

    def add(
        category: str,
        item: str,
        value_a: Any,
        value_b: Any,
        impact: str,
        notes: str = "",
    ) -> None:
        same = normalize_compare_value(value_a) == normalize_compare_value(value_b)
        rows.append(
            {
                "category": category,
                "item": item,
                "run_a_value": format_value(value_a),
                "run_b_value": format_value(value_b),
                "same": "true" if same else "false",
                "direct_comparison_impact": impact,
                "notes": notes,
            }
        )

    add("input", "reads", params_a.get("reads"), params_b.get("reads"), "blocking_if_different")
    add("input", "max_reads", params_a.get("max_reads"), params_b.get("max_reads"), "blocking_if_different")
    add(
        "input",
        "max_read_bases",
        params_a.get("max_read_bases"),
        params_b.get("max_read_bases"),
        "blocking_if_different",
    )
    add(
        "metadata",
        "pipeline_summary_present",
        bool(a.pipeline_summary_path),
        bool(b.pipeline_summary_path),
        "informational",
        "Standalone discover benchmarks may not have pipeline_summary.tsv.",
    )

    for key in RESULT_AFFECTING_DISCOVER_PARAMETERS:
        add(
            "discover_parameter",
            key,
            params_a.get(key),
            params_b.get(key),
            "blocking_if_different",
            discover_parameter_note(key),
        )

    add("result", "candidate_reads_count", a.candidate_count, b.candidate_count, "result_difference")
    add("result", "family_count", len(a.family_rows), len(b.family_rows), "result_difference")
    add(
        "result",
        "monomer_lengths_bp",
        sorted(a.monomer_lengths),
        sorted(b.monomer_lengths),
        "result_difference",
    )
    add(
        "result",
        "monomer_lengths_only_in_run_a",
        sorted(set(a.monomer_lengths) - set(b.monomer_lengths)),
        [],
        "result_difference",
    )
    add(
        "result",
        "monomer_lengths_only_in_run_b",
        [],
        sorted(set(b.monomer_lengths) - set(a.monomer_lengths)),
        "result_difference",
    )

    blocking = [
        row
        for row in rows
        if row["direct_comparison_impact"] == "blocking_if_different" and row["same"] == "false"
    ]
    reason = "; ".join(f"{row['item']}: {row['run_a_value']} vs {row['run_b_value']}" for row in blocking)
    rows.append(
        {
            "category": "interpretation",
            "item": "directly_comparable",
            "run_a_value": "true" if not blocking else "false",
            "run_b_value": "true" if not blocking else "false",
            "same": "true",
            "direct_comparison_impact": "summary",
            "notes": "Result catalogs are directly comparable." if not blocking else f"Not directly comparable: {reason}",
        }
    )
    if blocking:
        rows.append(
            {
                "category": "interpretation",
                "item": "reason_not_directly_comparable",
                "run_a_value": "",
                "run_b_value": "",
                "same": "false",
                "direct_comparison_impact": "summary",
                "notes": reason,
            }
        )
    return rows


def discover_parameter_note(key: str) -> str:
    if key == "min_support_reads":
        return "Different support thresholds can remove low-support families from one catalog."
    if key in {"min_monomer_len", "max_monomer_len", "top_periods"}:
        return "Different period-search settings can change which monomers are detected."
    if key in {"max_reads", "max_read_bases", "sample_rate", "seed"}:
        return "Different read subsets can change candidate and family counts."
    return ""


def parse_run_config(path: Path) -> dict[str, Any]:
    """Parse TandemX's simple generated run_config.yaml without a YAML dependency."""
    result: dict[str, Any] = {}
    section: str | None = None
    current_list_key: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip():
            continue
        if raw_line.startswith("  - ") and section and current_list_key:
            result.setdefault(current_list_key, []).append(parse_scalar(raw_line[4:]))
            continue
        if raw_line.startswith("  ") and section:
            key, value = split_key_value(raw_line.strip())
            result.setdefault(section, {})[key] = parse_scalar(value)
            continue
        key, value = split_key_value(raw_line)
        if value == "":
            section = key
            current_list_key = key
            if key == "parameters":
                result[key] = {}
            else:
                result[key] = []
            continue
        section = None
        current_list_key = None
        result[key] = parse_scalar(value)
    return result


def split_key_value(line: str) -> tuple[str, str]:
    if ":" not in line:
        raise ValueError(f"Invalid run_config.yaml line: {line}")
    key, value = line.split(":", 1)
    return key.strip(), value.strip()


def parse_scalar(value: str) -> Any:
    if value == "":
        return ""
    if value == "null":
        return None
    if value == "true":
        return True
    if value == "false":
        return False
    if value == "[]":
        return []
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1].replace('\\"', '"').replace("\\\\", "\\")
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_tsv_records(path: Path) -> int:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def read_monomer_lengths(path: Path, family_rows: Iterable[dict[str, str]]) -> list[int]:
    lengths: list[int] = []
    for row in family_rows:
        value = row.get("monomer_length_bp", "")
        if value:
            lengths.append(int(value))
    if lengths:
        return lengths
    return [len(record.sequence) for record in read_sequence_records(path)]


def normalize_compare_value(value: Any) -> Any:
    if isinstance(value, list):
        return [normalize_compare_value(item) for item in value]
    return value


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return ",".join(format_value(item) for item in value)
    return str(value)


def write_compare_tsv(path: Path, rows: Iterable[dict[str, str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=COMPARE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_compare_markdown(path: Path, a: RunArtifacts, b: RunArtifacts, rows: list[dict[str, str]]) -> None:
    row_by_item = {row["item"]: row for row in rows}
    direct = row_by_item["directly_comparable"]["run_a_value"]
    reason = row_by_item["directly_comparable"]["notes"]
    different_params = [
        row
        for row in rows
        if row["category"] == "discover_parameter" and row["same"] == "false"
    ]

    lines = [
        "# TandemX run comparison",
        "",
        "This report checks whether two TandemX outputs can be interpreted as directly comparable.",
        "",
        "## Inputs",
        "",
        f"- Run A: `{a.root_dir}`",
        f"- Run B: `{b.root_dir}`",
        f"- Run A discover directory: `{a.discover_dir}`",
        f"- Run B discover directory: `{b.discover_dir}`",
        f"- Same reads: `{row_by_item['reads']['same']}`",
        f"- Same max_reads: `{row_by_item['max_reads']['same']}`",
        f"- Same min-support-reads: `{row_by_item['min_support_reads']['same']}`",
        "",
        "## Results",
        "",
        f"- Candidate reads: {row_by_item['candidate_reads_count']['run_a_value']} vs {row_by_item['candidate_reads_count']['run_b_value']}",
        f"- Family count: {row_by_item['family_count']['run_a_value']} vs {row_by_item['family_count']['run_b_value']}",
        f"- Monomer lengths bp: `{row_by_item['monomer_lengths_bp']['run_a_value']}` vs `{row_by_item['monomer_lengths_bp']['run_b_value']}`",
        "",
        "## Discover parameter differences",
        "",
    ]
    if different_params:
        lines.extend(
            f"- `{row['item']}`: `{row['run_a_value']}` vs `{row['run_b_value']}`"
            + (f" — {row['notes']}" if row["notes"] else "")
            for row in different_params
        )
    else:
        lines.append("- No result-affecting discover parameter differences detected.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Directly comparable: `{direct}`",
            f"- Reason: {reason}",
            "",
            "If `directly_comparable` is false, treat result differences as parameter-driven until rerun with matched settings.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
