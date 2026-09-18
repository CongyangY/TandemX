"""Run the frozen M2 truth ladder without changing the alignment prototype."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

from benchmarks.m2_routes.alignment.prototype import decompose
from benchmarks.scripts.build_native_read_collapse import sha256_file


ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "benchmarks/m2_routes/truth_ladder_v1_20260918"
B1 = ROOT / "benchmarks/controlled_collapse/v3/development_bundle/inputs.jsonl"
L4 = ROOT / "benchmarks/m2_routes/native_equal_length_evidence_20260917/summary.json"
L5 = ROOT / "benchmarks/controlled_collapse/macadamia_bp_provisional_v1/m2_score/summary.json"
RC = str.maketrans("ACGT", "TGCA")
DNA = "ACGT"


def mutate_copy(sequence: str) -> str:
    """Apply the four position-fixed noise operations in protocol order."""
    bases = list(sequence)
    for position in (15, 70):
        bases[position] = DNA[(DNA.index(bases[position]) + 1) % 4]
    bases.insert(40, "A")
    bases.pop(60)
    return "".join(bases)


def assemble(pattern: str, monomers: dict[str, str], noisy: bool) -> tuple[str, list[str]]:
    sequence = []
    labels = []
    for name in pattern.split():
        orientation = "-" if name.endswith("-") else "+"
        label = name.rstrip("-")
        unit = monomers[label]
        if orientation == "-":
            unit = unit.translate(RC)[::-1]
        if noisy:
            unit = mutate_copy(unit)
        sequence.append(unit)
        labels.append(label + orientation)
    return "".join(sequence), labels


def evaluate(level: int, case_id: str, pattern: str, monomers: dict[str, str],
             noisy: bool) -> tuple[dict, dict]:
    sequence, truth_labels = assemble(pattern, monomers, noisy)
    result = decompose(sequence, monomers)
    recovered = [copy.label + copy.orientation for copy in result.copies]
    source = {"level": level, "case_id": case_id, "pattern": pattern,
              "sequence": sequence, "monomers": monomers,
              "truth_labels": truth_labels}
    score = {"level": level, "case_id": case_id,
             "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
             "length_bp": len(sequence), "state": result.state,
             "reason": result.reason, "score": result.score,
             "alternative_score": result.alternative_score,
             "truth_labels": truth_labels, "recovered_labels": recovered,
             "recovered_copy_count": len(result.copies),
             "pass": result.state == "RESOLVED" and recovered == truth_labels}
    return source, score


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError("Truth ladder output exists")
    protocol_path = FROZEN / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if protocol["status"] != "frozen_before_ladder_execution":
        raise ValueError("Truth ladder not frozen")
    for path, key in ((ROOT / "benchmarks/m2_routes/alignment/prototype.py", "prototype_sha256"),
                      (B1, "source_inputs_jsonl_sha256")):
        section = protocol if key == "prototype_sha256" else protocol["level_3"]
        if sha256_file(path) != section[key]:
            raise ValueError(f"Frozen source changed: {path}")
    for path, level in ((L4, "level_4"), (L5, "level_5")):
        if sha256_file(path) != protocol[level]["existing_summary_sha256"]:
            raise ValueError(f"Frozen upper-level evidence changed: {path}")
    rng = random.Random(protocol["level_1"]["random_seed"])
    length = protocol["level_1"]["template_length_bp"]
    synthetic = {label: "".join(rng.choice(DNA) for _ in range(length)) for label in "ABCD"}
    with B1.open() as stream:
        real = json.loads(stream.readline())["candidate_monomers"]
    if set(real) != {"C1", "C2", "C3"} or {len(v) for v in real.values()} != {178}:
        raise ValueError("Frozen real monomer catalogue changed")
    cases = [
        (1, "rule_ABC", "A B C A B C", synthetic, False),
        (1, "rule_ABCD", "A B C D A B C D", synthetic, False),
        (2, "noisy_ABC", "A B C A B C", synthetic, True),
        (2, "noisy_B_reverse", "A B C A B- C", synthetic, True),
        (3, "plant_consensus_HOR", "C1 C2 C3 C1 C2 C3", real, False),
    ]
    inputs, rows = [], []
    for level, case_id, pattern, monomers, noisy in cases:
        source, score = evaluate(level, case_id, pattern, monomers, noisy)
        inputs.append(source)
        rows.append(score)
    l4 = json.loads(L4.read_text())
    l5 = json.loads(L5.read_text())
    upper = {
        "level_4": {"source_summary_sha256": sha256_file(L4),
                    "positive_detected": l4["m2_technical"]["positive_detected"],
                    "positive_count": l4["positive_count"],
                    "negative_supported": l4["m2_technical"]["negative_supported"],
                    "negative_count": l4["negative_count"],
                    "pass": l4["m2_technical"]["positive_detected"] == 5 and
                            l4["m2_technical"]["negative_supported"] == 1},
        "level_5": {"source_summary_sha256": sha256_file(L5),
                    "technical_coverage": l5["technical_discordant_count"] +
                                          l5["technical_supported_count"],
                    "case_count": l5["case_count"],
                    "original_read_resolved_count": l5["original_read_resolved_count"],
                    "pass": l5["technical_abstain_count"] == 0 and
                            l5["technical_discordant_count"] == 8 and
                            l5["technical_supported_count"] == 1,
                    "biological_accuracy": "not_evaluated"},
    }
    output.mkdir(parents=True, exist_ok=False)
    with (output / "inputs.jsonl").open("w") as stream:
        for row in inputs:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    with (output / "per_case.jsonl").open("w") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    summary = {"protocol_sha256": sha256_file(protocol_path),
               "input_sha256": sha256_file(output / "inputs.jsonl"),
               "per_case_sha256": sha256_file(output / "per_case.jsonl"),
               "level_1_pass": all(row["pass"] for row in rows if row["level"] == 1),
               "level_2_pass": all(row["pass"] for row in rows if row["level"] == 2),
               "level_3_pass": all(row["pass"] for row in rows if row["level"] == 3),
               **upper, "status": "development_diagnostic_only"}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
