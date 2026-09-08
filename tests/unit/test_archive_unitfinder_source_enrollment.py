import gzip
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.archive_unitfinder_source_enrollment import archive
from benchmarks.scripts.download_unitfinder_zh13_assembly import file_hashes


def test_archive_failed_source_keeps_large_partial_external(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    partial = source / "genome.fa.gz.partial"
    partial.write_bytes(b"partial")
    hashes = file_hashes(partial)
    receipt = source / "run_receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "complete": False,
                "fate": "source_enrollment_failure",
                "partial_file": {"file": partial.name, **hashes},
            }
        )
    )
    diagnosis = source / "transport_diagnosis.json"
    diagnosis.write_text(
        json.dumps(
            {
                "run_receipt_sha256": hashlib.sha256(receipt.read_bytes()).hexdigest()
            }
        )
    )
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "experiment_id": "test",
                "assembly": {"accession": "GWH-TEST"},
                "boundary": "source only",
            }
        )
    )
    outdir = tmp_path / "archive"
    manifest = archive(source, config, outdir)
    assert manifest["experiment_complete"] is False
    assert not (outdir / partial.name).exists()
    headline = json.loads((outdir / "headline_summary.json").read_text())
    assert headline["external_payload"]["sha256"] == hashes["sha256"]
    assert headline["external_payload"]["retained_outside_git"] is True
