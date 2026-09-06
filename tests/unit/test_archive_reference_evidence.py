import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.scripts.archive_reference_evidence import archive


def make_reference(root: Path) -> Path:
    root.mkdir()
    metadata = root / "source_metadata.html"
    metadata.write_text("<p>published source</p>\n")
    reference = root / "reference.fa"
    reference.write_text(">chr1\nACGTNN\n>chr2\nACGT\n")
    reference_hash = hashlib.sha256(reference.read_bytes()).hexdigest()
    plan = {
        "filename": reference.name,
        "expected_sha256": reference_hash,
        "max_bytes": 1000,
        "metadata_sha256": digest_file(metadata),
    }
    plan_path = root / "reference_plan.json"
    plan_path.write_text(json.dumps(plan) + "\n")
    receipt = {
        "complete": True,
        "plan_sha256": digest_file(plan_path),
        "warning": "donor_identity_unresolved",
        "transfer": {"bytes": reference.stat().st_size, "sha256": reference_hash},
        "qc": {
            "complete": True,
            "input_sha256": reference_hash,
            "contig_count": 2,
            "total_bases": 10,
            "n_bases": 2,
            "other_ambiguous_bases": 0,
            "base_counts": {"A": 2, "C": 2, "G": 2, "T": 2, "N": 2},
            "contigs": [
                {"contig": "chr1", "length_bp": 6},
                {"contig": "chr2", "length_bp": 4},
            ],
        },
    }
    (root / "reference_receipt.json").write_text(json.dumps(receipt) + "\n")
    return reference


def test_archive_reference_copies_compact_evidence_and_rechecks_external_fasta(tmp_path: Path) -> None:
    source = tmp_path / "source"
    reference = make_reference(source)
    outdir = tmp_path / "archive"
    manifest = archive(source, outdir)
    assert manifest["complete"] is True
    assert len(manifest["archived_files"]) == 3
    assert manifest["external_reference"]["copied_to_repository"] is False
    assert manifest["external_reference"]["sha256"] == digest_file(reference)
    assert not (outdir / reference.name).exists()
    assert json.loads((outdir / "archive_manifest.json").read_text()) == manifest
    with pytest.raises(ValueError, match="already exists"):
        archive(source, outdir)


@pytest.mark.parametrize("fault", ["reference", "contigs", "metadata"])
def test_archive_reference_rejects_changed_or_inconsistent_evidence(tmp_path: Path, fault: str) -> None:
    source = tmp_path / "source"
    reference = make_reference(source)
    if fault == "reference":
        reference.write_text(">chr1\nAAAA\n")
    elif fault == "metadata":
        (source / "source_metadata.html").write_text("changed\n")
    else:
        path = source / "reference_receipt.json"
        receipt = json.loads(path.read_text())
        receipt["qc"]["contigs"][0]["length_bp"] = 5
        path.write_text(json.dumps(receipt) + "\n")
    with pytest.raises(ValueError):
        archive(source, tmp_path / "archive")
