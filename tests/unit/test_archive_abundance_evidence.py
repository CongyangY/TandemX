import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.archive_abundance_evidence import DETECTOR_FILES, ROOT_FILES, archive


def build_result(tmp_path: Path, scenario_count: int = 1) -> Path:
    source = tmp_path / "result"
    snapshot = source / "source_snapshot"
    snapshot.mkdir(parents=True)
    config = {
        "seeds": {"development": [1], "heldout": [2]},
        "periods": [11], "copies": [5], "coverages": [1],
        "substitution_rates": [0], "assembly_fractions": [1, 0],
    }
    if scenario_count == 2:
        config.update(unit_substitution_rates=[0, .1], array_fragment_counts=[1])
    (source / "run_config.json").write_text(json.dumps(config))
    benchmark_files = {"benchmarks/abundance/run.py"}
    hashes = {}
    benchmark_hashes = {}
    for name in sorted(set(DETECTOR_FILES) | benchmark_files):
        path = snapshot / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"source={name}\n")
        target = benchmark_hashes if name in benchmark_files else hashes
        target[name] = digest_file(path)
    environment = {
        "git_head": "a" * 40,
        "split": "heldout",
        "source_snapshot": str(snapshot),
        "file_hashes": hashes,
        "benchmark_source_sha256": benchmark_hashes,
        "source_digest": hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
        "config_sha256": digest_file(source / "run_config.json"),
    }
    (source / "environment.json").write_text(json.dumps(environment))
    validation = {"complete": True, "executions": 4 * scenario_count,
                  "successful": 4 * scenario_count,
                  "challenge_scenarios": scenario_count,
                  "copy_number_family_rows": scenario_count,
                  "localization_family_rows": 2 * scenario_count,
                  "comparison_family_rows": 2 * scenario_count}
    (source / "validation.json").write_text(json.dumps(validation))
    (source / "run.log").write_text("complete\n")
    write_table(source / "copy_number_metrics.tsv",
                [{"seed": 2, "family_id": "f1"} for _ in range(scenario_count)],
                ["seed", "family_id"])
    write_table(source / "localization_metrics.tsv",
                [{"seed": 2, "family_id": "f1"} for _ in range(2 * scenario_count)],
                ["seed", "family_id"])
    write_table(source / "comparison_metrics.tsv",
                [{"seed": 2, "family_id": "f1", "outcome": value}
                 for _ in range(scenario_count) for value in ("TN", "TP")],
                ["seed", "family_id", "outcome"])
    write_table(source / "copy_number_summary.tsv",
                [{"coverage": 1, "substitution_rate": 0} for _ in range(scenario_count)],
                ["coverage", "substitution_rate"])
    write_table(source / "comparison_summary.tsv",
                [{"coverage": 1, "substitution_rate": 0, "assembly_fraction": fraction,
                  "TP": int(fraction == 0), "FN": 0, "FP": 0, "TN": int(fraction == 1)}
                 for _ in range(scenario_count) for fraction in (1, 0)],
                ["coverage", "substitution_rate", "assembly_fraction", "TP", "FN", "FP", "TN"])
    labels = ("locate", "locate", "quantify", "compare", "compare") * scenario_count
    for index, label in enumerate(labels):
        # The configured matrix has five commands: two locate, one quantify and two compare.
        path = source / "runs" / str(index) / "receipt.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps({"label": label, "command": ["python", label], "exit_code": 0,
                                    "runtime_seconds": .1, "peak_rss_mib": 2,
                                    "cpu_user_seconds": .05, "cpu_system_seconds": .01,
                                    "timed_out": False}))
    validation["executions"] = validation["successful"] = 5 * scenario_count
    (source / "validation.json").write_text(json.dumps(validation))
    return source


def test_archive_validates_matrix_sources_and_receipts(tmp_path: Path) -> None:
    source = build_result(tmp_path)
    outdir = tmp_path / "archive"
    manifest = archive(source, outdir)
    assert len(manifest) == len(ROOT_FILES) + len(set(DETECTOR_FILES) | {"benchmarks/abundance/run.py"}) + 1
    resources = (outdir / "resource_metrics.tsv").read_text().splitlines()
    assert len(resources) == 6
    assert (outdir / "source_snapshot/tandemx/compare/mvp.py").is_file()
    for row in manifest:
        assert digest_file(outdir / row["file"]) == row["sha256"]


def test_archive_accepts_expanded_challenge_matrix(tmp_path: Path) -> None:
    source = build_result(tmp_path, scenario_count=2)
    outdir = tmp_path / "archive"
    archive(source, outdir)
    assert len((outdir / "resource_metrics.tsv").read_text().splitlines()) == 11


@pytest.mark.parametrize("fault", ["matrix", "source", "receipt", "seeds"])
def test_archive_rejects_incomplete_or_changed_evidence(tmp_path: Path, fault: str) -> None:
    source = build_result(tmp_path)
    if fault == "matrix":
        validation = json.loads((source / "validation.json").read_text())
        validation["comparison_family_rows"] = 1
        (source / "validation.json").write_text(json.dumps(validation))
    elif fault == "source":
        (source / "source_snapshot/tandemx/compare/mvp.py").write_text("changed\n")
    elif fault == "receipt":
        receipt = next((source / "runs").glob("**/receipt.json"))
        data = json.loads(receipt.read_text())
        data["exit_code"] = 1
        receipt.write_text(json.dumps(data))
    else:
        rows = (source / "copy_number_metrics.tsv").read_text().replace("\n2\t", "\n3\t")
        (source / "copy_number_metrics.tsv").write_text(rows)
    with pytest.raises(ValueError):
        archive(source, tmp_path / "archive")
