"""Small source-independent tests for the preregistration input gate."""

from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "benchmarks/scripts/validate_input_manifest.py"
SPEC = importlib.util.spec_from_file_location("validate_input_manifest", SCRIPT)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def entry(path: Path, *, dataset: str = "d1", split: str = "development") -> dict:
    data = path.read_bytes()
    import hashlib

    return {
        "dataset": dataset, "sample": "sample-1", "species": "Species one",
        "benchmark_task": "raw-read-discovery",
        "source_provenance_status": "verified",
        "source_provenance_evidence": ["source receipt"],
        "source_pairing_status": "not-applicable",
        "source_pairing_evidence": [],
        "truth_type": "none", "family": [], "array_coordinates": [],
        "coverage": None, "read_source": "synthetic test", "assembly_source": "none",
        "split": split, "allowed_tuning": split == "development", "metrics": [],
        "competitors": [], "software_versions": {}, "donor_id": "donor-1",
        "accessions": ["RUN1"], "assembly_lineage": "assembly-1",
        "family_homology_group": "group-1", "status": "eligible",
        "files": [{"path": str(path), "kind": "fastq", "role": "reads", "bytes": len(data),
                   "sha256": hashlib.sha256(data).hexdigest()}],
    }


def test_valid_gzip_fastq_and_crc(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b"@read1\nACGT\n+\nIIII\n")
    result = module.validate({"schema_version": 2, "datasets": [entry(path)]})
    assert result["valid"]
    assert result["datasets_checked"][0]["files"][0]["records"] == 1
    path.write_bytes(path.read_bytes()[:-3])
    bad = entry(path)
    assert not module.validate({"schema_version": 2, "datasets": [bad]})["valid"]


def test_malformed_fastq_and_hash_mismatch(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@read1\nACGT\n+\nIII\n")
    assert not module.validate({"schema_version": 2, "datasets": [entry(path)]})["valid"]
    path.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    row = entry(path)
    row["files"][0]["sha256"] = "0" * 64
    assert "mismatch" in module.validate({"schema_version": 2, "datasets": [row]})["errors"][0]


@pytest.mark.parametrize("shared", ["donor_id", "sample", "accessions", "assembly_lineage", "family_homology_group", "species"])
def test_related_inputs_cannot_cross_splits(tmp_path: Path, shared: str) -> None:
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    first = entry(path)
    second = entry(path, dataset="d2", split="validation")
    for key, value in (("donor_id", "donor-2"), ("sample", "sample-2"),
                       ("accessions", ["RUN2"]), ("assembly_lineage", "assembly-2"),
                       ("family_homology_group", "group-2"), ("species", "Species two")):
        if key != shared:
            second[key] = value
    errors = module.validate({"schema_version": 2, "datasets": [first, second]})["errors"]
    assert any("split leakage" in error for error in errors)


def test_heldout_tuning_rejected(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    row = entry(path, split="final-heldout")
    row["allowed_tuning"] = True
    result = module.validate({"schema_version": 2, "datasets": [row]})
    assert not result["valid"]
    assert "tuning" in result["errors"][0]


@pytest.mark.parametrize("same_path", [True, False])
def test_renamed_metadata_same_file_cannot_cross_splits(tmp_path: Path, same_path: bool) -> None:
    first_path = tmp_path / "first.fastq"
    first_path.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    second_path = first_path if same_path else tmp_path / "copy.fastq"
    if second_path != first_path:
        second_path.write_bytes(first_path.read_bytes())
    first = entry(first_path)
    second = entry(second_path, dataset="renamed", split="final-heldout")
    second.update(sample="new-sample", species="Other species", donor_id="new-donor",
                  accessions=["NEW_RUN"], assembly_lineage="new-assembly",
                  family_homology_group="new-family-group")
    errors = module.validate({"schema_version": 2, "datasets": [first, second]})["errors"]
    assert any("file_sha256" in error for error in errors)
    assert any("file_path" in error for error in errors) is same_path


def test_eligible_requires_task_roles_and_pairing(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    row = entry(reads)
    row["benchmark_task"] = "read-assembly-deficit"
    result = module.validate({"schema_version": 2, "datasets": [row]})
    assert "missing file roles" in result["errors"][0]

    assembly = tmp_path / "assembly.fa"
    assembly.write_bytes(b">chr1\nACGT\n")
    import hashlib

    row["files"].append({"path": str(assembly), "kind": "fasta", "role": "assembly",
                         "bytes": assembly.stat().st_size,
                         "sha256": hashlib.sha256(assembly.read_bytes()).hexdigest()})
    row["source_pairing_status"] = "unresolved"
    result = module.validate({"schema_version": 2, "datasets": [row]})
    assert "verified pairing" in result["errors"][0]
    row["source_pairing_status"] = "verified"
    row["source_pairing_evidence"] = ["same donor and haplotype evidence receipt"]
    assert module.validate({"schema_version": 2, "datasets": [row]})["valid"]


def test_integrity_inventory_never_eligible(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    row = entry(path)
    row["benchmark_task"] = "input-integrity-only"
    assert "integrity-only" in module.validate({"schema_version": 2, "datasets": [row]})["errors"][0]
