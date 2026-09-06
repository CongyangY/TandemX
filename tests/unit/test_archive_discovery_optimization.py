import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_discovery_optimization import validate_pair


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n")


def test_discovery_optimization_archive_requires_hash_checked_parity(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    replay = tmp_path / "replay"
    sequence_hash = "a" * 64
    baseline_environment = baseline / "environment.json"
    baseline_execution = baseline / "tandemx" / "execution.json"
    summary = {
        "processed_reads": 3, "processed_bases": 300,
        "candidate_count": 2, "family_count": 1,
    }
    _write_json(baseline_environment, {"input": {"fasta_sha256": sequence_hash}})
    _write_json(baseline_execution, {
        "exit_code": 0, "timed_out": False,
        "runtime_seconds": 4, "peak_rss_mib": 20,
    })
    _write_json(baseline / "tandemx" / "discover" / "discovery_summary.json", summary)
    _write_json(replay / "discover" / "discovery_summary.json", summary)
    products = {}
    for index in range(6):
        name = f"product_{index}.tsv"
        old = baseline / "tandemx" / "discover" / name
        new = replay / "discover" / name
        old.write_text(f"{index}\n")
        new.write_text(f"{index}\n")
        products[name] = {
            "expected_sha256": digest_file(old),
            "observed_sha256": digest_file(new),
            "byte_identical": True,
        }
    _write_json(replay / "environment.json", {
        "previous_run": str(baseline.resolve()),
        "baseline_environment_sha256": digest_file(baseline_environment),
        "baseline_execution_sha256": digest_file(baseline_execution),
        "input_sha256": sequence_hash,
    })
    validation_path = replay / "validation.json"
    _write_json(validation_path, {
        "complete": True,
        "execution": {
            "exit_code": 0, "timed_out": False,
            "runtime_seconds": 2, "peak_rss_mib": 10,
        },
        "products": products,
        "uncompared_new_products": [],
    })
    checked = validate_pair("toy", baseline, replay)
    assert checked["row"]["speedup"] == 2
    assert checked["row"]["peak_rss_change_percent"] == -50

    validation = json.loads(validation_path.read_text())
    validation["products"]["product_0.tsv"]["byte_identical"] = False
    _write_json(validation_path, validation)
    with pytest.raises(ValueError, match="Incomplete or mismatched"):
        validate_pair("toy", baseline, replay)
