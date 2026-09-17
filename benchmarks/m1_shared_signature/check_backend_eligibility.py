"""Audit whether the frozen M1 mapping baseline can replace public quantify.

This checks evidence identity and interface semantics, without opening a new
estimator search or using any held-out M1/M2 material.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmarks.m1_shared_signature.model import align_and_gate

ROOT = Path(__file__).resolve().parents[2]
ACD = ROOT / "benchmarks/m1_shared_signature/evidence_acd_20260916/results.json"
OCCUPANCY = ROOT / "benchmarks/m1_shared_signature/evidence_occupancy_research_20260916/revised_12_case/results.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_identity(evidence: dict, prefix: str) -> bool:
    return all(
        sha256(ROOT / prefix / name) == digest
        for name, digest in evidence["source_sha256"].items()
    )


def rejects_full_read() -> bool:
    """The frozen unit-sized mapper cannot classify a two-unit read."""
    try:
        align_and_gate(["ACGT" * 40], {"f1": "ACGT" * 20}, 18)
    except ValueError as exc:
        return "equal-length" in str(exc)
    return False


def refuses_identical_families() -> bool:
    accepted, mapped, rejected, tied = align_and_gate(
        ["ACGT" * 20], {"f1": "ACGT" * 20, "f2": "ACGT" * 20}, 18
    )
    return len(accepted) == 1 and not rejected and tied == 1 and all(
        value == 0 for value in mapped.values()
    )


def audit() -> dict:
    acd = json.loads(ACD.read_text())
    occupancy = json.loads(OCCUPANCY.read_text())
    cases = acd["cases"]
    rows = occupancy["rows"]
    same_inputs = len(cases) == len(rows) == 12 and {
        (row["scenario"], row["seed"], row["reads_sha256"], row["catalogue_sha256"])
        for row in rows
    } == {
        (row["scenario"], row["seed"], row["reads_sha256"], row["catalogue_sha256"])
        for row in cases
    }
    scenarios = sorted({row["scenario"] for row in rows})
    mare = {
        scenario: {
            method: sum(row[field] for row in rows if row["scenario"] == scenario) / 3
            for method, field in (
                ("ordinary_mapping", "ordinary_positive_mare"),
                ("occupancy", "occupancy_positive_mare"),
            )
        }
        for scenario in scenarios
    }
    checks = {
        "frozen_acd_source_matches_receipt": source_identity(acd, "benchmarks/m1_shared_signature"),
        "frozen_occupancy_source_matches_receipt": source_identity(occupancy, ""),
        "same_12_development_inputs": same_inputs,
        "mapping_refuses_identical_family_attribution": refuses_identical_families(),
        "mapping_accepts_full_reads": not rejects_full_read(),
        "occupancy_dominates_mapping_all_identifiable_conditions": all(
            values["occupancy"] <= values["ordinary_mapping"]
            for scenario, values in mare.items() if not scenario.startswith("identical")
        ),
        "independent_full_read_copy_bp_validation_in_m1_receipts": any(
            "independent_full_read_copy_bp_validation" in record for record in (acd, occupancy)
        ),
        "paired_public_backend_runtime_rss_in_m1_receipts": any(
            "paired_public_backend_runtime_rss" in record for record in (acd, occupancy)
        ),
    }
    return {
        "status": "eligible_for_public_backend" if all(checks.values()) else "not_eligible",
        "evidence_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (ACD, OCCUPANCY)},
        "checks": checks,
        "development_positive_mare": mare,
        "scope": "Existing frozen M1 mapping and occupancy implementations only; no held-out data opened.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"])


if __name__ == "__main__":
    main()
