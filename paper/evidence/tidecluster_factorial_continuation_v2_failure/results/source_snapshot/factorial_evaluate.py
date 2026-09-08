"""Score TideCluster against frozen planted factorial assembly truth."""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from benchmarks.challenge.adapters import read_fasta
from benchmarks.challenge.evaluate import score_arrays
from benchmarks.challenge.run import json_safe
from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.challenge.sequence_metrics import score_threshold_recovery
from benchmarks.tidecluster.normalize import (
    fasta_lengths,
    normalize_resolved_tidecluster,
    read_truth,
)


NORMALIZED_FIELDS = [
    "sequence_id", "start", "end", "family_id", "period", "consensus_sequence",
    "copy_number", "representative_tidehunter_id", "representative_selection_source",
    "supporting_tidehunter_count", "copy_number_source",
]
MATCH_FIELDS = [
    "prediction_index", "read_id", "start", "end", "period",
    "matched_truth_index", "truth_family_id", "iou", "status",
]
FAMILY_FIELDS = [
    "truth_id", "recovered", "assigned_sequence_index", "threshold", "criterion",
]


def evaluate_factorial(
    tidehunter_gff: Path,
    intermediate_clustering_gff: Path,
    clustering_gff: Path,
    family_consensus_fasta: Path,
    genome_dir: Path,
    outdir: Path,
) -> dict[str, object]:
    """Evaluate valid empty or non-empty native output without hiding failures."""
    if outdir.exists():
        raise FileExistsError(f"Choose a new output directory: {outdir}")
    manifest_path = genome_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("generator") != "streamed_factorial_genome_v1":
        raise ValueError("Require an independently generated factorial genome")
    assembly = genome_dir / "genome.fa"
    truth_path = genome_dir / "truth_copy_number.tsv"
    catalogue_path = genome_dir / "catalogue.fa"
    generated_inputs = {
        assembly: manifest["files"]["genome.fa"],
        truth_path: manifest["files"]["truth_copy_number.tsv"],
        catalogue_path: manifest["files"]["catalogue.fa"],
    }
    native_inputs = (
        tidehunter_gff,
        intermediate_clustering_gff,
        clustering_gff,
        family_consensus_fasta,
    )
    if any(not path.is_file() for path in (*generated_inputs, *native_inputs)):
        raise ValueError("A required generated or native input is missing")
    if any(digest_file(path) != expected for path, expected in generated_inputs.items()):
        raise ValueError("Factorial genome input differs from its frozen manifest")

    records = normalize_resolved_tidecluster(
        tidehunter_gff,
        intermediate_clustering_gff,
        clustering_gff,
        family_consensus_fasta,
        allow_empty=True,
    )
    founders = read_fasta(catalogue_path)
    truth = read_truth(truth_path, founders)
    metrics, matches = score_arrays(
        [record.array() for record in records], truth, fasta_lengths(assembly), 0.5
    )
    family_metrics, family_rows = score_threshold_recovery(
        [record.consensus_sequence for record in records], founders, 0.9
    )
    metrics.update(family_metrics)
    metrics.update(
        {
            "seed": int(manifest["seed"]),
            "operational_family_count": len({record.family_id for record in records}),
            "normalized_coordinate_system": "0-based_half-open",
            "source_coordinate_system": "GFF3_1-based_inclusive",
            "warning": (
                "known_planted_factorial_assembly_truth;same_process_validation_genome;"
                "not_an_independent_plant;TideCluster_settings_frozen_before_output_inspection"
            ),
        }
    )
    normalized_rows = [asdict(record) for record in records]
    outdir.mkdir(parents=True)
    write_table(outdir / "normalized_arrays.tsv", normalized_rows, NORMALIZED_FIELDS)
    write_table(outdir / "array_matches.tsv", matches, MATCH_FIELDS)
    write_table(outdir / "family_recovery.tsv", family_rows, FAMILY_FIELDS)
    (outdir / "metrics.json").write_text(
        json.dumps(json_safe(metrics), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    input_hashes = {
        str(path.resolve()): {"bytes": path.stat().st_size, "sha256": digest_file(path)}
        for path in (*generated_inputs, manifest_path, *native_inputs)
    }
    if any(digest_file(path) != expected for path, expected in generated_inputs.items()):
        raise ValueError("Factorial genome input changed during evaluation")
    outputs = {
        name: digest_file(outdir / name)
        for name in ("normalized_arrays.tsv", "array_matches.tsv", "family_recovery.tsv", "metrics.json")
    }
    receipt = {
        "schema_version": 1,
        "complete": True,
        "input_hashes": input_hashes,
        "outputs": outputs,
        "metrics": metrics,
    }
    (outdir / "evaluation_receipt.json").write_text(
        json.dumps(json_safe(receipt), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metrics
