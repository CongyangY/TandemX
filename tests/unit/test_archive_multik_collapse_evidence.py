import csv
import hashlib
import json
from pathlib import Path

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_multik_collapse_evidence import SOURCE_FILES, archive


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fixture(tmp_path: Path) -> Path:
    source = tmp_path / "result"
    snapshot = source / "source_snapshot"
    snapshot.mkdir(parents=True)
    file_hashes = {}
    benchmark_hashes = {}
    for name in SOURCE_FILES:
        path = snapshot / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
        target = benchmark_hashes if name.startswith("benchmarks/") else file_hashes
        target[name] = digest_file(path)
    config = {
        "seeds": {"development": [1], "heldout": [2]},
        "collapse_model": {
            "method": "multik_loglinear_else_single_k21", "k_values": [15, 21, 27, 31],
            "low_depth_cutoff": 2.0, "low_depth_threshold": 0.5,
            "standard_depth_threshold": 0.6, "unavailable_rule": "fallback_single_k21",
            "calibration_seeds": [1],
        },
    }
    calibration = tmp_path / "calibration"
    calibration.mkdir()
    for name in ("validation.json", "calibration.tsv", "environment.json"):
        (calibration / name).write_text(name)
    config["collapse_model"]["calibration_result"] = str(calibration)
    for name, key in {
        "validation.json": "calibration_validation_sha256",
        "calibration.tsv": "calibration_table_sha256",
        "environment.json": "calibration_environment_sha256",
    }.items():
        config["collapse_model"][key] = digest_file(calibration / name)
    (source / "run_config.json").write_text(json.dumps(config))
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    for name in ("validation.json", "environment.json"):
        (baseline / name).write_text("baseline " + name)
    (baseline / "run_config.json").write_text((source / "run_config.json").read_text())
    environment = {
        "scope": "predeclared_heldout_known_catalogue_point_estimate_validation",
        "source_snapshot": str(snapshot), "file_hashes": file_hashes,
        "source_digest": hashlib.sha256(json.dumps(file_hashes, sort_keys=True).encode()).hexdigest(),
        "benchmark_source_sha256": benchmark_hashes, "previous_result": str(baseline),
        "previous_validation_sha256": digest_file(baseline / "validation.json"),
        "previous_environment_sha256": digest_file(baseline / "environment.json"),
        "previous_config_sha256": digest_file(source / "run_config.json"),
    }
    (source / "environment.json").write_text(json.dumps(environment))
    rows = []
    outcomes = {"single_k21": "FN", "multik_loglinear": "TP",
                "multik_fallback_depth_rule": "TP"}
    for method, outcome in outcomes.items():
        rows.append({"seed": 2, "coverage": 1, "substitution_rate": 0,
                     "assembly_fraction": 0.5, "family_id": "f1", "method": method,
                     "outcome": outcome, "decision_threshold": 0.5 if "fallback" in method else 0.6})
    write_tsv(source / "comparison_metrics.tsv", rows)
    write_tsv(source / "comparison_summary.tsv", [{"method": name} for name in outcomes])
    write_tsv(source / "calibration.tsv", [{
        "fold": "predeclared_heldout_validation", "training_seeds": "1",
        "evaluation_seeds": "2", "selected_low_depth_threshold": 0.5,
        "baseline_TP": 0, "baseline_FN": 1, "baseline_FP": 0, "baseline_TN": 0,
        "calibrated_TP": 1, "calibrated_FN": 0, "calibrated_FP": 0, "calibrated_TN": 0,
    }])
    confusion = {
        "single_k21": {"available": 1, "unavailable": 0, "TP": 0, "FN": 1, "FP": 0, "TN": 0},
        "multik_loglinear": {"available": 1, "unavailable": 0, "TP": 1, "FN": 0, "FP": 0, "TN": 0},
        "multik_fallback_depth_rule": {"available": 1, "unavailable": 0, "TP": 1, "FN": 0, "FP": 0, "TN": 0},
    }
    (source / "validation.json").write_text(json.dumps({
        "complete": True, "split": "heldout", "heldout_used": True,
        "leave_one_genome_out_folds": 0, "comparison_family_rows": 3,
        "expected_comparison_family_rows": 3, "method_confusion": confusion,
    }))
    return source


def test_archive_copies_and_hashes_compact_heldout_evidence(tmp_path: Path) -> None:
    source = fixture(tmp_path)
    outdir = tmp_path / "archive"
    manifest = archive(source, outdir)
    assert len(manifest) == 19
    assert all(digest_file(outdir / row["file"]) == row["sha256"] for row in manifest)
    assert (outdir / "calibration_source/calibration.tsv").is_file()
    assert (outdir / "baseline_receipt/validation.json").is_file()


def test_archive_rejects_changed_confusion_counts(tmp_path: Path) -> None:
    source = fixture(tmp_path)
    validation = json.loads((source / "validation.json").read_text())
    validation["method_confusion"]["single_k21"]["FN"] = 0
    (source / "validation.json").write_text(json.dumps(validation))
    try:
        archive(source, tmp_path / "archive")
    except ValueError as error:
        assert "confusion" in str(error)
    else:
        raise AssertionError("Changed held-out confusion counts were accepted")
