"""Run an explicitly engineered component challenge for the M2 graph route."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from time import perf_counter
import tracemalloc

from . import prototype
from .prototype import ReadPath, audit_transition_graph


MONOMERS = {"A": "ACGTTGCA", "B": "GGAATCCA", "C": "CATGAGTC",
            "A_snp": "ACGTTGCG"}


def _reads(labels: tuple[str, ...], number: int = 3, *,
           confidence: float = 1.0, profile: str | None = None,
           haplotype: str | None = None) -> tuple[ReadPath, ...]:
    return tuple(ReadPath(f"molecule-{i}", labels, (confidence,) * len(labels),
                          profile or f"profile-{i % 2}", haplotype=haplotype)
                 for i in range(number))


def _case(case_id: str, assembly: str, observed: str | tuple[str, ...],
          truth: str, *, number: int = 3, confidence: float = 1.0,
          profile: str | None = None) -> dict[str, object]:
    return {"case_id": case_id, "assembly": tuple(assembly),
            "reads": _reads(tuple(observed), number, confidence=confidence,
                            profile=profile), "truth_scope": truth}


def development_cases() -> list[dict[str, object]]:
    """All labels and error metadata are supplied; this is not raw-read validation."""
    cases = [
        _case("intact", "ABCABC", "ABCABC", "intact_transition_graph"),
        _case("copy_compression", "ABC", "ABCABC", "discordant_transition_graph"),
        _case("order_rewiring", "ABAC", "ACAB", "discordant_transition_graph"),
        _case("subclass_replacement", "ABC", "ABA", "discordant_transition_graph"),
        _case("equal_edges_different_order", "ABACA", "ACABA",
              "discordant_full_order_but_equal_transition_graph"),
        _case("close_monomer_variant_or_systematic_motif_error", "ABA",
              ("A", "A_snp", "A"), "identity_unresolved"),
        _case("one_molecule", "ABC", "ABCABC", "insufficient_molecules", number=1),
        _case("substitution_assignment_uncertain", "ABC", "ABCABC",
              "upstream_assignment_uncertain", confidence=0.75),
        _case("indel_assignment_uncertain", "ABC", "ABCABC",
              "upstream_assignment_uncertain", confidence=0.70),
        _case("homopolymer_assignment_uncertain", "ABC", "ABCABC",
              "upstream_assignment_uncertain", confidence=0.55),
        _case("correlated_error_profile", "ABC", "ABCABC",
              "profile_confounded", profile="one_basecaller"),
        _case("within_label_snp_unseen", "ABA", "ABA",
              "within_label_variant_unobservable"),
    ]
    mixed = _case("mixed_haplotypes", "ABC", "ABCABC",
                  "mixed_haplotypes")
    mixed["reads"] = _reads(tuple("ABCABC"), haplotype="h1") + tuple(
        ReadPath(f"h2-{i}", tuple("ABC"), (1.0,) * 3,
                 f"profile-{i % 2}", haplotype="h2") for i in range(3))
    cases.append(mixed)
    rare = _case("rare_variant_single_molecule", "ABC", "ABC",
                 "one_minor_molecule_unresolved_biology", number=4)
    rare["reads"] = rare["reads"] + (
        ReadPath("rare-molecule", tuple("ABA"), (1.0,) * 3, "profile-1"),)
    cases.append(rare)
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    cases = development_cases()
    case_payload = [{"case_id": case["case_id"], "assembly": case["assembly"],
                     "reads": [asdict(read) for read in case["reads"]],
                     "truth_scope": case["truth_scope"]} for case in cases]
    input_sha256 = hashlib.sha256(json.dumps(case_payload, sort_keys=True).encode()).hexdigest()
    rows = []
    for case in cases:
        tracemalloc.start()
        start = perf_counter()
        result = audit_transition_graph(case["assembly"], case["reads"],
                                        monomer_sequences=MONOMERS)
        elapsed = perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rows.append({"case_id": case["case_id"], "truth_scope": case["truth_scope"],
                     "result": asdict(result), "elapsed_seconds": elapsed,
                     "tracemalloc_peak_bytes": peak})
    output = {
        "status": "development_engineered_label_component_only",
        "input_sha256": input_sha256,
        "prototype_sha256": hashlib.sha256(Path(prototype.__file__).read_bytes()).hexdigest(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "case_count": len(rows),
        "counts_by_status": {status: sum(row["result"]["status"] == status for row in rows)
                             for status in ("SUPPORTED", "DISCORDANT", "AMBIGUOUS",
                                            "INSUFFICIENT_READ_SUPPORT")},
        "max_case_elapsed_seconds": max(row["elapsed_seconds"] for row in rows),
        "max_case_tracemalloc_peak_bytes": max(row["tracemalloc_peak_bytes"] for row in rows),
        "limitations": [
            "engineered labels and assignment confidence are supplied, not inferred from reads",
            "two error profiles are names in metadata, not independent sequencing validation",
            "source sequences are 8 bp motifs and cannot model realistic centromeric arrays",
            "tracemalloc measures Python allocations, not process RSS",
            "cases were developed alongside the algorithm and are not validation or held-out data",
        ],
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output["counts_by_status"], sort_keys=True))


if __name__ == "__main__":
    main()
