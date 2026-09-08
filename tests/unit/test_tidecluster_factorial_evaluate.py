import json
from pathlib import Path

import pytest

from benchmarks.challenge.schema import digest_file
from benchmarks.tidecluster.factorial_evaluate import evaluate_factorial


def _genome(tmp_path: Path) -> Path:
    genome = tmp_path / "genome"
    genome.mkdir()
    (genome / "genome.fa").write_text(">chr1\n" + "A" * 10 + "ACGT" * 10 + "C" * 10 + "\n")
    (genome / "catalogue.fa").write_text(">f1\nACGT\n")
    (genome / "truth_copy_number.tsv").write_text(
        "chrom\tfamily_id\tstart\tend\tperiod\nchr1\tf1\t10\t50\t4\n"
    )
    files = {
        name: digest_file(genome / name)
        for name in ("genome.fa", "catalogue.fa", "truth_copy_number.tsv")
    }
    (genome / "manifest.json").write_text(
        json.dumps({"generator": "streamed_factorial_genome_v1", "seed": 1, "files": files})
    )
    return genome


def _native(tmp_path: Path, empty: bool = False) -> tuple[Path, Path, Path, Path]:
    tidehunter = tmp_path / "tidehunter.gff"
    intermediate = tmp_path / "intermediate.gff"
    final = tmp_path / "final.gff"
    consensus = tmp_path / "consensus.fa"
    header = "##gff-version 3\n"
    tidehunter.write_text(
        header if empty else header +
        "chr1\tTideHunter\ttandem_repeat\t11\t50\t1\t.\t.\t"
        "ID=rep1;consensus_sequence=ACGT;consensus_length=4;copy_number=10\n"
    )
    intermediate.write_text(
        header if empty else header +
        "chr1\tTideCluster\ttandem_repeat\t11\t50\t1\t.\t.\tName=rep1\n"
    )
    final.write_text(
        header if empty else header +
        "chr1\tTideCluster\ttandem_repeat\t11\t50\t1\t.\t.\tName=TRC_1\n"
    )
    consensus.write_text("" if empty else ">TRC_1_rep1\nACGT\n")
    return tidehunter, intermediate, final, consensus


def test_factorial_evaluator_scores_resolved_native_output(tmp_path: Path) -> None:
    metrics = evaluate_factorial(*_native(tmp_path), _genome(tmp_path), tmp_path / "out")
    assert metrics["array_recall"] == 1
    assert metrics["cyclic_monomer_recall"] == 1
    assert metrics["operational_family_count"] == 1


def test_factorial_evaluator_counts_valid_empty_output_as_zero_recall(tmp_path: Path) -> None:
    metrics = evaluate_factorial(*_native(tmp_path, empty=True), _genome(tmp_path), tmp_path / "out")
    assert metrics["predicted_array_count"] == 0
    assert metrics["array_recall"] == 0
    assert metrics["cyclic_monomer_recall"] == 0
    assert (tmp_path / "out/normalized_arrays.tsv").read_text().count("\n") == 1


def test_factorial_evaluator_rejects_changed_truth_and_overwrite(tmp_path: Path) -> None:
    native = _native(tmp_path)
    genome = _genome(tmp_path)
    (genome / "catalogue.fa").write_text(">f1\nAAAA\n")
    with pytest.raises(ValueError, match="frozen manifest"):
        evaluate_factorial(*native, genome, tmp_path / "out")
