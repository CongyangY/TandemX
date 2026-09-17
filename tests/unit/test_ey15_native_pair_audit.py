"""Guard the small original-read Ey15-2 source archive against drift."""

import gzip
import hashlib
import json
from pathlib import Path


ARCHIVE = Path(__file__).parents[2] / "benchmarks/controlled_collapse/ey15_native_pair_audit_20260917"


def test_ey15_source_archive_is_structurally_and_hash_consistent():
    manifest = json.loads((ARCHIVE / "source_eligibility_manifest.json").read_text())
    extract = json.loads((ARCHIVE / "raw_read_extract_receipt.json").read_text())
    for name, expected in manifest["artifact_sha256"].items():
        assert hashlib.sha256((ARCHIVE / name).read_bytes()).hexdigest() == expected

    context = (ARCHIVE / "ey15_9994_context.fa").read_text().splitlines()[1]
    assert len(context) == 9244
    array = context[3000:6244]
    assert len(array) == 3244 and set(array) <= set("ACGT")
    identity = sum(a == b for a, b in zip(array[:-420], array[420:])) / (len(array) - 420)
    assert abs(identity - manifest["array"]["period_shift_identity"]) < 1e-12

    records = {}
    with gzip.open(ARCHIVE / "ey15_native_spanners.fastq.gz", "rb") as handle:
        while header := handle.readline():
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            assert header.startswith(b"@") and separator.startswith(b"+")
            assert len(sequence.strip()) == len(quality.strip())
            read_id = header[1:].split()[0].decode()
            assert read_id not in records
            records[read_id] = header + sequence + separator + quality
    assert len(records) == 7
    expected_records = {r["id"]: r["record_sha256"] for r in extract["reads"]}
    assert {rid: hashlib.sha256(record).hexdigest() for rid, record in records.items()} == expected_records
    zmws = {record.splitlines()[0].split()[1].split(b"/")[1] for record in records.values()}
    assert len(zmws) == 7

    full = [line.split("\t") for line in (ARCHIVE / "ey15_native_spanners_full_assembly.paf").read_text().splitlines()]
    assert len(full) == 7 and {row[0] for row in full} == set(records)
    for row in full:
        assert row[5] == "Chr1" and int(row[7]) <= 12828234 and int(row[8]) >= 12833478
        assert int(row[11]) >= 20 and int(row[9]) / int(row[10]) >= 0.99

    assert manifest["split"] == "development"
    assert manifest["final_heldout_enrolled"] is False
    assert manifest["benchmark_eligibility"]["physical_missing_copy_truth"] == "not_available"
