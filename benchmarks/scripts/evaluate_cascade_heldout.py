"""Apply the frozen cascade promotion gates to a completed held-out run."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import yaml

from benchmarks.challenge.schema import digest_file, write_table


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _true(value: object) -> bool:
    return str(value).lower() == "true"


def _geomean(values: list[float]) -> float | None:
    return math.exp(sum(math.log(value) for value in values) / len(values)) if values and all(value > 0 for value in values) else None


def _complete_min(rows: list[dict[str, str]], field: str) -> float | None:
    values = [_number(row.get(field)) for row in rows]
    return min(value for value in values if value is not None) if values and all(value is not None for value in values) else None


def _complete_max(rows: list[dict[str, str]], field: str) -> float | None:
    values = [_number(row.get(field)) for row in rows]
    return max(value for value in values if value is not None) if values and all(value is not None for value in values) else None


def _gate(name: str, observed: object, operator: str, threshold: float) -> dict[str, object]:
    value = _number(observed)
    passed = value is not None and (value <= threshold if operator == "<=" else value >= threshold)
    return {"name": name, "passed": passed, "observed": observed, "operator": operator, "threshold": threshold}


def validate_matrix(config: dict[str, Any], config_path: Path, run: Path, raw: list[dict[str, str]], summary: list[dict[str, str]]) -> dict[str, int]:
    required = ("environment.json", "run_config.yaml", "validation.json", "raw_runs.tsv", "summary.tsv")
    if any(not (run / name).is_file() for name in required):
        raise ValueError("Held-out run lacks required root files")
    environment = json.loads((run / "environment.json").read_text())
    run_config = yaml.safe_load((run / "run_config.yaml").read_text())
    validation = json.loads((run / "validation.json").read_text())
    if environment.get("config_sha256") != digest_file(config_path):
        raise ValueError("Run environment does not match the frozen config hash")
    if environment.get("split") != "heldout" or run_config.get("selected_split") != "heldout":
        raise ValueError("Run is not the held-out split")
    if run_config.get("selected_scenarios") is not None:
        raise ValueError("Held-out run omitted one or more predeclared scenarios")
    scenarios = {row["name"] for row in config["scenarios"]}
    seeds = {int(seed) for seed in config["seeds"]["heldout"]}
    tools = set(config["tools"])
    repetitions = int(config["repetitions"])
    expected_raw = len(scenarios) * len(seeds) * len(tools) * repetitions
    expected_summary = len(scenarios) * len(seeds) * len(tools)
    raw_keys = {(row["scenario"], int(row["seed"]), row["tool"], int(row["repetition"])) for row in raw}
    summary_keys = {(row["scenario"], int(row["seed"]), row["tool"]) for row in summary}
    expected_raw_keys = {(scenario, seed, tool, repetition) for scenario in scenarios for seed in seeds for tool in tools for repetition in range(1, repetitions + 1)}
    expected_summary_keys = {(scenario, seed, tool) for scenario in scenarios for seed in seeds for tool in tools}
    if raw_keys != expected_raw_keys or len(raw) != expected_raw:
        raise ValueError("Held-out raw run matrix is incomplete or duplicated")
    if summary_keys != expected_summary_keys or len(summary) != expected_summary:
        raise ValueError("Held-out summary matrix is incomplete or duplicated")
    if validation.get("complete") is not True or validation.get("total_runs") != expected_raw:
        raise ValueError("Held-out validation receipt disagrees with the matrix")
    return {"raw_rows": expected_raw, "summary_rows": expected_summary, "datasets": len(scenarios) * len(seeds)}


def evaluate(config_path: Path, run: Path, outdir: Path) -> dict[str, object]:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f"Choose a new empty output directory: {outdir}")
    config = yaml.safe_load(config_path.read_text())
    raw = _read_tsv(run / "raw_runs.tsv")
    summary = _read_tsv(run / "summary.tsv")
    matrix = validate_matrix(config, config_path, run, raw, summary)
    limits = config["acceptance_gates"]
    scenarios = {row["name"]: row for row in config["scenarios"]}
    positive = {name for name, row in scenarios.items() if float(row.get("positive_fraction", 1.0)) > 0}
    negative = set(scenarios) - positive
    tandemx_raw = [row for row in raw if row["tool"] == "tandemx"]
    comparator_raw = [row for row in raw if row["tool"] != "tandemx"]
    tandemx = [row for row in summary if row["tool"] == "tandemx"]
    tandemx_positive = [row for row in tandemx if row["scenario"] in positive]
    tandemx_negative = [row for row in tandemx if row["scenario"] in negative]
    related = [row for row in tandemx if row["scenario"] == "related_families"]

    gates = [
        _gate("tandemx_failed_runs", sum(row["status"] != "ok" for row in tandemx_raw), "<=", limits["tandemx_failed_runs_max"]),
        _gate("all_comparator_failed_runs", sum(row["status"] != "ok" for row in comparator_raw), "<=", limits["all_comparator_failed_runs_max"]),
        _gate("tandemx_deterministic_groups_fraction", sum(_true(row["deterministic"]) for row in tandemx) / len(tandemx), ">=", limits["tandemx_deterministic_groups_fraction_min"]),
        _gate("minimum_positive_array_recall", _complete_min(tandemx_positive, "array_recall"), ">=", limits["positive_array_recall_min"]),
        _gate("minimum_positive_array_precision", _complete_min(tandemx_positive, "array_precision"), ">=", limits["positive_array_precision_min"]),
        _gate("maximum_negative_read_call_rate", _complete_max(tandemx_negative, "negative_read_call_rate"), "<=", limits["negative_read_call_rate_max"]),
        _gate("minimum_related_family_cyclic_monomer_recall", _complete_min(related, "cyclic_monomer_recall"), ">=", limits["related_family_cyclic_monomer_recall_min"]),
    ]
    # A missing metric must fail rather than disappear from min/max calculations.
    for gate, rows, field in (
        (gates[3], tandemx_positive, "array_recall"),
        (gates[4], tandemx_positive, "array_precision"),
        (gates[5], tandemx_negative, "negative_read_call_rate"),
        (gates[6], related, "cyclic_monomer_recall"),
    ):
        if not rows or any(_number(row.get(field)) is None for row in rows):
            gate["passed"] = False
            gate["missing_group_count"] = sum(_number(row.get(field)) is None for row in rows)

    by_key = {(row["scenario"], row["seed"], row["tool"]): row for row in summary}
    paired_rows: list[dict[str, object]] = []
    for scenario in sorted(scenarios):
        for seed in map(str, config["seeds"]["heldout"]):
            tx = by_key[(scenario, seed, "tandemx")]
            th = by_key[(scenario, seed, "tidehunter")]
            valid = (
                _true(tx["deterministic"]) and _true(th["deterministic"])
                and int(tx["successful_runs"]) == int(tx["attempted_runs"])
                and int(th["successful_runs"]) == int(th["attempted_runs"])
            )
            row: dict[str, object] = {"scenario": scenario, "seed": seed, "positive": scenario in positive, "valid_pair": valid}
            for field in ("array_recall", "array_precision"):
                tx_value, th_value = _number(tx.get(field)), _number(th.get(field))
                row[f"tandemx_{field}"] = tx_value
                row[f"tidehunter_{field}"] = th_value
                row[f"{field}_difference"] = tx_value - th_value if tx_value is not None and th_value is not None else None
            for field in ("median_runtime_seconds", "median_peak_rss_mib"):
                tx_value, th_value = _number(tx.get(field)), _number(th.get(field))
                row[f"tandemx_{field}"] = tx_value
                row[f"tidehunter_{field}"] = th_value
                row[f"{field}_ratio"] = tx_value / th_value if tx_value is not None and th_value and th_value > 0 else None
            paired_rows.append(row)
    valid_pairs = [row for row in paired_rows if row["valid_pair"]]
    positive_pairs = [row for row in valid_pairs if row["positive"]]
    expected_pairs = matrix["datasets"]
    paired_fraction = len(valid_pairs) / expected_pairs
    recall_differences = [row["array_recall_difference"] for row in positive_pairs]
    precision_differences = [row["array_precision_difference"] for row in positive_pairs]
    runtime_ratios = [row["median_runtime_seconds_ratio"] for row in valid_pairs]
    rss_ratios = [row["median_peak_rss_mib_ratio"] for row in valid_pairs]
    gates.extend([
        _gate("paired_tidehunter_groups_fraction", paired_fraction, ">=", limits["paired_tidehunter_groups_fraction_min"]),
        _gate("worst_tidehunter_array_recall_difference", min(recall_differences) if recall_differences and all(value is not None for value in recall_differences) else None, ">=", -float(limits["tidehunter_array_recall_noninferiority_margin"])),
        _gate("worst_tidehunter_array_precision_difference", min(precision_differences) if precision_differences and all(value is not None for value in precision_differences) else None, ">=", -float(limits["tidehunter_array_precision_noninferiority_margin"])),
        _gate("tidehunter_runtime_geometric_mean_ratio", _geomean([float(value) for value in runtime_ratios if value is not None]) if len(runtime_ratios) == expected_pairs else None, "<=", limits["tidehunter_runtime_geometric_mean_ratio_max"]),
        _gate("tidehunter_peak_rss_geometric_mean_ratio", _geomean([float(value) for value in rss_ratios if value is not None]) if len(rss_ratios) == expected_pairs else None, "<=", limits["tidehunter_peak_rss_geometric_mean_ratio_max"]),
    ])
    result = {
        "benchmark_id": config["benchmark_id"],
        "status": "passed" if all(gate["passed"] for gate in gates) else "failed",
        "matrix": matrix,
        "config_sha256": digest_file(config_path),
        "run_source_digest": json.loads((run / "environment.json").read_text())["source_digest"],
        "gates": gates,
        "failed_gate_names": [gate["name"] for gate in gates if not gate["passed"]],
        "warning": "heldout_seeds_consumed_once;failed_gates_must_be_retained;timing_repetitions_are_not_biological_replicates",
    }
    outdir.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in paired_rows for key in row))
    write_table(outdir / "paired_tidehunter.tsv", ({key: row.get(key, "NA") for key in fields} for row in paired_rows), fields)
    (outdir / "gate_results.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = evaluate(args.config.resolve(), args.run.resolve(), args.outdir.resolve())
    except (OSError, ValueError, KeyError, ZeroDivisionError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
