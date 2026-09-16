import csv
import hashlib
import json
from pathlib import Path

import pytest

from benchmarks.scripts.controlled_bp_edit import generate, read_fragment


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "paper/evidence/legume_ysd56_candidate_v1"
SOURCE = ARCHIVE / "TXF000708_array.fasta"
EXTRACTION = ARCHIVE / "array_extraction_receipt.json"
MANIFEST = ARCHIVE / "MANIFEST.tsv"
CONFIRMATION = ARCHIVE / "array_confirmation_summary.tsv"


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_row(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


def test_real_archive_source_and_bp_only_truth():
    region_id, seq, region, period = read_fragment(SOURCE, EXTRACTION, MANIFEST, CONFIRMATION)
    assert (region_id, region["sequence_id"], region["start0"], region["end0"]) == (
        "TXF000708_array", "CP154579.1", 13966881, 14090910)
    assert len(seq) == 124029
    assert (period, len(seq) % period) == (785, 784)
    assert set(seq) == {"A", "C", "G", "T"}
    assert hash_file(SOURCE) == "0866dbf818cbcd81b7a768c974fb21f221b3de95cfd20e95b9314613b5b04df9"


def test_generated_gradient_independent_reconstruction(tmp_path):
    out = tmp_path / "generated"
    generate(SOURCE, EXTRACTION, MANIFEST, CONFIRMATION, out)
    source = "".join(SOURCE.read_text().splitlines()[1:])
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["source_fragment_length_bp"] == 124029
    assert receipt["complete_unit_count_truth"] == "unknown"
    assert receipt["read_pairing_status"] == "not_evaluated"
    assert len(receipt["cases"]) == 5
    for level in (0, 25, 50, 75, 100):
        case = f"bp_deletion_{level:03d}"
        row = one_row(out / f"{case}.ledger.tsv")
        retained = 124029 - (124029 * level // 100)
        assert int(row["injected_deleted_bp"]) == 124029 - retained
        assert (int(row["deleted_fragment_start0"]), int(row["deleted_fragment_end0"])) == (retained, 124029)
        assert (int(row["deleted_genomic_start0"]), int(row["deleted_genomic_end0"])) == (
            13966881 + retained, 14090910)
        assert row["complete_unit_count_truth"] == "unknown"
        assert row["truth_scope"] == "injected_bp_delta_only"
        lines = (out / f"{case}.fa").read_text().splitlines()
        assert lines[0].startswith(f">{case} ")
        assert "".join(lines[1:]) == source[:retained]
        with (out / f"{case}.chain.tsv").open(newline="") as handle:
            chain = list(csv.DictReader(handle, delimiter="\t"))
        assert int(chain[-1]["source_end0"]) == 124029
        assert int(chain[-1]["edited_end0"]) == retained
        if level == 100:
            assert len(lines) == 1
            assert receipt["cases"][level // 25]["empty_fasta_record"] is True
        case_receipt = next(item for item in receipt["cases"] if item["case_id"] == case)
        for suffix, key in (("fa", "fasta_sha256"), ("ledger.tsv", "ledger_sha256"),
                            ("chain.tsv", "chain_sha256")):
            assert hash_file(out / f"{case}.{suffix}") == case_receipt[key]


def test_source_hash_or_header_mismatch_blocks_edit(tmp_path):
    changed = tmp_path / "changed.fa"
    changed.write_text(SOURCE.read_text().replace("AATTTA", "CATTT A".replace(" ", ""), 1))
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        read_fragment(changed, EXTRACTION, MANIFEST, CONFIRMATION)
    bad_receipt = tmp_path / "bad_receipt.json"
    receipt = json.loads(EXTRACTION.read_text())
    receipt["extracted"][0]["start0"] += 1
    bad_receipt.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="disagree"):
        read_fragment(SOURCE, bad_receipt, MANIFEST, CONFIRMATION)
