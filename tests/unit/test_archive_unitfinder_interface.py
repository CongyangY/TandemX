import hashlib
import json

import pytest

from benchmarks.scripts.archive_unitfinder_interface import archive, validate


def _config(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(
        json.dumps(
            {
                "experiment_id": "unitfinder_smoke",
                "source_commit": "abc123",
                "boundary": "interface only",
            }
        )
    )
    return path


def test_archive_failed_build_retains_log_and_no_accuracy(tmp_path) -> None:
    source = tmp_path / "failed"
    source.mkdir()
    log = source / "docker_build.log"
    log.write_text("failed\n")
    (source / "run_receipt.json").write_text(
        json.dumps(
            {
                "complete": False,
                "fate": "container_build_failure",
                "accuracy_available": False,
                "smoke_execution_started": False,
                "docker_build_log": {
                    "file": log.name,
                    "bytes": log.stat().st_size,
                    "sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                },
            }
        )
    )
    manifest = archive(source, _config(tmp_path), tmp_path / "archive")
    assert manifest["experiment_complete"] is False
    headline = json.loads(
        (tmp_path / "archive/headline_summary.json").read_text()
    )
    assert headline["accuracy_available"] is False
    assert (tmp_path / "archive/results/docker_build.log").is_file()


def test_validate_failed_build_rejects_changed_log(tmp_path) -> None:
    source = tmp_path / "failed"
    source.mkdir()
    log = source / "docker_build.log"
    log.write_text("failed\n")
    (source / "run_receipt.json").write_text(
        json.dumps(
            {
                "complete": False,
                "fate": "container_build_failure",
                "smoke_execution_started": False,
                "docker_build_log": {
                    "file": log.name,
                    "bytes": log.stat().st_size,
                    "sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                },
            }
        )
    )
    log.write_text("changed\n")
    with pytest.raises(ValueError, match="log changed"):
        validate(source, _config(tmp_path))


def test_archive_successful_build_remains_build_only(tmp_path) -> None:
    source = tmp_path / "built"
    provenance = source / "image_provenance"
    provenance.mkdir(parents=True)
    log = source / "docker_build.log"
    log.write_text("built\n")
    tracked = provenance / "unitfinder_tracked_files.sha256"
    tracked.write_text("a" * 64 + "  README.md\n")
    (source / "run_receipt.json").write_text(
        json.dumps(
            {
                "complete": True,
                "fate": "container_build_passed",
                "accuracy_available": False,
                "smoke_execution_started": False,
                "docker_build_log": {
                    "file": log.name,
                    "bytes": log.stat().st_size,
                    "sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                },
                "image_provenance": {
                    tracked.name: {
                        "bytes": tracked.stat().st_size,
                        "sha256": hashlib.sha256(tracked.read_bytes()).hexdigest(),
                    }
                },
            }
        )
    )
    outdir = tmp_path / "archive"
    archive(source, _config(tmp_path), outdir)
    headline = json.loads((outdir / "headline_summary.json").read_text())
    assert headline["evidence_scope"] == "container_build_only"
    assert headline["accuracy_available"] is False
