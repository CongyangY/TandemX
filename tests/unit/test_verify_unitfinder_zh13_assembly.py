import gzip
import hashlib
import json
from pathlib import Path

from benchmarks.scripts.verify_unitfinder_zh13_assembly import inspect, verify


def test_inspect_and_verify_zh13_source_independently(tmp_path: Path) -> None:
    result = tmp_path / "result"
    result.mkdir()
    payload = result / "genome.fa.gz"
    with gzip.open(payload, "wt") as handle:
        handle.write(">chr01 note\nACGTNN\n>chr02\nTTAA\n")
    summary, lengths = inspect(payload)
    assert summary["record_count"] == 2
    assert summary["other_bases"] == 2
    assert lengths == [("chr01", 6), ("chr02", 4)]
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "assembly": {
                    "expected_bytes": summary["bytes"],
                    "expected_md5": summary["md5"],
                }
            }
        )
    )
    (result / "run_receipt.json").write_text(
        json.dumps(
            {
                "complete": True,
                "fate": "source_enrollment_passed",
                "compressed_file": {
                    "file": payload.name,
                    "bytes": summary["bytes"],
                    "md5": summary["md5"],
                    "sha256": summary["sha256"],
                },
                "fasta_qc": {
                    "record_count": summary["record_count"],
                    "base_count": summary["base_count"],
                    "acgt_bases": summary["acgt_bases"],
                    "other_bases": summary["other_bases"],
                },
            }
        )
    )
    checked = verify(config, result)
    assert checked["complete"] is True
    assert hashlib.sha256((result / "sequence_lengths.tsv").read_bytes()).hexdigest() == (
        checked["sequence_lengths_file"]["sha256"]
    )
