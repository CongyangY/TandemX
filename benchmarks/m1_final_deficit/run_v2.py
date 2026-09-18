"""Correct the documented forward-only baseline defect without changing v1."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import tempfile

from benchmarks.abundance.simulate import GenomeSpec, build_genome, sample_reads
from benchmarks.m1_final_deficit.run import sha, summarize
from benchmarks.m1_shared_signature.occupancy_research import CompetitiveConfig, classify_window
from tandemx.io.sequences import read_sequence_records_many

ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "protocol_v2.json"
PROTOCOL_SHA256 = "c0e585559208e5b188fdc4f5e5e13d5927f75bd3b4d9c5dfc79dd9143fc1f5eb"


def run(v1_path: Path, outdir: Path) -> None:
    if outdir.exists():
        raise FileExistsError(outdir)
    raw_protocol = PROTOCOL.read_bytes()
    if sha(raw_protocol) != PROTOCOL_SHA256:
        raise ValueError("Frozen v2 baseline protocol changed")
    p = json.loads(raw_protocol)
    original = (ROOT / "protocol.json").read_bytes()
    if sha(original) != p["v1_protocol_sha256"]:
        raise ValueError("Original M1 panel protocol changed")
    v1_bytes = v1_path.read_bytes()
    if sha(v1_bytes) != p["v1_results_sha256"]:
        raise ValueError("v1 result archive differs")
    v1 = json.loads(v1_bytes)
    old = json.loads(original)
    outdir.mkdir(parents=True)
    rows, panels = [], []
    config = CompetitiveConfig(reads=Path("unused"), monomers=Path("unused"),
                               genome_size=1, outdir=Path("unused"),
                               window_bp=80, min_tail_bp=80,
                               max_edit_fraction=.225, min_margin_edits=1)
    with tempfile.TemporaryDirectory(prefix="m1_corrected_mapping_") as root:
        for panel in v1["panels"]:
            seed = panel["source_seed"]
            rate = panel["unit_substitution_rate"]
            coverage = panel["nominal_coverage"]
            spec = GenomeSpec(seed=seed, periods=tuple(old["periods_bp"]),
                              copies=tuple(old["source_copies"]), flank_bp=old["flank_bp"],
                              unit_substitution_rate=rate)
            genome, units, _ = build_genome(spec)
            sampling = sample_reads(genome, [], Path(root) / panel["panel_id"],
                                    seed=old["read_seed_offset"]+seed, coverage=coverage,
                                    read_length=old["read_length_bp"],
                                    substitution_rate=old["sequencing_substitution_rate"])
            if sampling["files"]["reads.fa"] != panel["reads_sha256"]:
                raise ValueError("Regenerated reads changed")
            assigned: Counter[str] = Counter()
            ambiguous = unknown = 0
            for record in read_sequence_records_many([Path(root) / panel["panel_id"] / "reads.fa"]):
                for offset in range(0, len(record.sequence), 80):
                    chunk = record.sequence[offset:offset+80]
                    state, family = classify_window(chunk, units, config)
                    if state == "assigned":
                        assigned[family] += len(chunk)
                    elif state == "ambiguous":
                        ambiguous += len(chunk)
                    else:
                        unknown += len(chunk)
            estimates = {family: assigned[family]/panel["actual_coverage"] for family in units}
            panels.append(dict(panel_id=panel["panel_id"], reads_sha256=panel["reads_sha256"],
                               estimates_bp=estimates, ambiguous_read_bp=ambiguous,
                               unknown_read_bp=unknown))
            for original_row in v1["rows"]:
                if original_row["panel_id"] != panel["panel_id"] or original_row["method"] != "production_quantify":
                    continue
                family = original_row["family"]
                estimate = estimates[family]
                predicted = estimate-original_row["edited_assembly_bp"]
                error = predicted-original_row["injected_missing_bp"]
                rows.append(dict(**{**original_row,
                                    "method": "strand_complete_competitive_mapping",
                                    "read_estimate_bp": estimate,
                                    "predicted_signed_deficit_bp": predicted,
                                    "signed_error_bp": error, "absolute_error_bp": abs(error),
                                    "relative_error_to_source": abs(error)/original_row["source_bp"]}))
    if len(rows) != 360:
        raise AssertionError("Corrected comparator denominator differs")
    result = dict(protocol_v2_sha256=sha(raw_protocol), v1_results_sha256=sha(v1_bytes),
                  v2_protocol_commit="7cfb915", rows=rows, panels=panels)
    (outdir / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    (outdir / "summary.json").write_text(json.dumps(summarize(v1["rows"] + rows, old), indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v1-results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.v1_results, arguments.outdir)
