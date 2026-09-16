"""Run a small exact-sequence M2 development challenge against a length baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from time import perf_counter

from benchmarks.architecture import prototype
from benchmarks.architecture.prototype import compare_anchored_reads


MONOMERS = {"A": "ACGTTGCA", "B": "GGAATCCA", "C": "CATGAGTC"}


def _array(labels: str) -> str:
    return "".join(MONOMERS[label] for label in labels)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _cases() -> list[tuple[str, str, str, bool]]:
    return [
        ("no_edit_ab", _array("ABAB"), _array("ABAB"), False),
        ("no_edit_abc", _array("ABCABC"), _array("ABCABC"), False),
        ("no_edit_variant", _array("ABAC"), _array("ABAC"), False),
        ("copy_deletion", _array("ABAB"), _array("ABABAB"), True),
        ("equal_length_variant", _array("ABAC"), _array("ABAB"), True),
        ("equal_length_label_swap", _array("ACAC"), _array("ABAB"), True),
        ("orientation_inversion", _array("ABAB"),
         MONOMERS["A"] + _reverse_complement(MONOMERS["B"]) + _array("AB"), True),
        ("hor_rearrangement", _array("ABAC"), _array("ACAB"), True),
        ("monomer_deletion", _array("ABA"), _array("ABAB"), True),
        ("monomer_duplication", _array("ABAAB"), _array("ABAB"), True),
        ("hor_compression", _array("ABAB"), _array("ABACABAC"), True),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = _cases()
    payload = json.dumps({"monomers": MONOMERS, "cases": cases},
                         sort_keys=True, separators=(",", ":")).encode()
    input_sha256 = hashlib.sha256(payload).hexdigest()
    rows = []
    for case_id, assembly, read, truth in cases:
        start = perf_counter()
        result = compare_anchored_reads(
            assembly, {"molecule-1": read, "molecule-2": read}, MONOMERS
        )
        elapsed = perf_counter() - start
        rows.append({
            "case": case_id,
            "truth_label_path_discordance": truth,
            "architecture_prediction": result.status == "candidate_discordance",
            "architecture_status": result.status,
            "length_only_prediction": len(read) != len(assembly),
            "assembly_bp": len(assembly),
            "read_bp": len(read),
            "elapsed_seconds": elapsed,
        })
    summary = {}
    for method in ("architecture", "length_only"):
        key = f"{method}_prediction"
        summary[method] = {
            "tp": sum(row[key] and row["truth_label_path_discordance"] for row in rows),
            "fn": sum(not row[key] and row["truth_label_path_discordance"] for row in rows),
            "fp": sum(row[key] and not row["truth_label_path_discordance"] for row in rows),
            "tn": sum(not row[key] and not row["truth_label_path_discordance"] for row in rows),
        }
    commit = subprocess.run(["git", "rev-parse", "HEAD"], check=True,
                            capture_output=True, text=True).stdout.strip()
    output = {
        "status": "development_exact_synthetic_only",
        "source_commit": commit,
        "source_sha256": {
            "prototype.py": hashlib.sha256(Path(prototype.__file__).read_bytes()).hexdigest(),
            "run_development.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "input_sha256": input_sha256,
        "baseline": "equal_input_anchored_read_vs_assembly_bp_length_only",
        "prototype": "same_inputs_monomer_decomposition_and_order",
        "warning": "identical_synthetic_read_sequences_are_not_biological_replication",
        "cases": rows,
        "summary": summary,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
