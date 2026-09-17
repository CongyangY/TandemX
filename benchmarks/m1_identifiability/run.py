"""Recalculate frozen M1 development-case and formal-catalogue diagnostics.

Example: conda run -n tandemx-dev python -m benchmarks.m1_identifiability.run
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from benchmarks.m1_shared_signature.run import PARAMETERS, SCENARIOS, make_case
from .diagnostics import catalogue_diagnostics, read_assignment_diagnostics


ROOT = Path(__file__).resolve().parents[2]
FORMAL = ROOT / "paper/evidence/srf_formal_unified_v1/datasets"
OUTPUT = ROOT / "benchmarks/m1_identifiability/evidence_20260917/results.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_small_fasta(path: Path) -> dict[str, str]:
    catalogue: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                catalogue[name] = "".join(chunks)
            name = line[1:].split()[0]
            if not name or name in catalogue:
                raise ValueError(f"Invalid or repeated FASTA identifier in {path}")
            chunks = []
        elif line.strip():
            if name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            chunks.append(line.strip().upper())
    if name is not None:
        catalogue[name] = "".join(chunks)
    if not catalogue:
        raise ValueError(f"Empty catalogue: {path}")
    return catalogue


def run() -> dict:
    cases = []
    for scenario in SCENARIOS:
        name, differences, error, family_error_shift = scenario
        for seed in PARAMETERS["seeds"]:
            catalogue, reads, labels = make_case(
                seed, differences, error, family_error_shift
            )
            cases.append({
                "case_id": f"{name}_s{seed}",
                "catalogue_sha256": hashlib.sha256(json.dumps(catalogue, sort_keys=True).encode()).hexdigest(),
                "reads_sha256": hashlib.sha256(json.dumps(reads, sort_keys=True).encode()).hexdigest(),
                "catalogue": catalogue_diagnostics(catalogue, PARAMETERS["k"]),
                "read_assignment": read_assignment_diagnostics(
                    catalogue, reads, labels,
                    round(PARAMETERS["unit_bp"] * PARAMETERS["gate_fraction"]),
                ),
            })
    formal = []
    for path in sorted(FORMAL.glob("*/truth_monomers.fa")):
        catalogue = parse_small_fasta(path)
        # These prior formal simulation inputs were already used elsewhere;
        # their sequence diagnostics are descriptive, not a new holdout.
        k = min(PARAMETERS["k"], min(map(len, catalogue.values())))
        formal.append({
            "dataset": path.parent.name,
            "fasta_sha256": sha256(path),
            "catalogue": catalogue_diagnostics(catalogue, k),
        })
    result = {
        "status": "development_diagnostic_not_independent_validation",
        "method": "noiseless_circular_kmer_linear_rank_plus_known_source_cyclic_hamming_assignment",
        "source_sha256": {
            path.name: sha256(path) for path in (
                Path(__file__), Path(__file__).with_name("diagnostics.py"),
                ROOT / "benchmarks/m1_shared_signature/run.py",
                ROOT / "benchmarks/m1_shared_signature/model.py",
            )
        },
        "k": PARAMETERS["k"],
        "gate_mismatches": round(PARAMETERS["unit_bp"] * PARAMETERS["gate_fraction"]),
        "m1_development_cases": cases,
        "prior_formal_simulation_catalogues": formal,
        "limits": [
            "Structural states condition on a supplied error-free family catalogue and circular k-mer features.",
            "Read attribution uses known synthetic source labels only for diagnosis, never to fit an estimator.",
            "Neither structural identifiability nor read assignment establishes physical copy-number truth.",
            "The 12 cases share only three simulation seeds and are not independent biological samples.",
            "Prior formal simulation catalogues are not independent held-out validation for a new gate.",
        ],
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    run()
