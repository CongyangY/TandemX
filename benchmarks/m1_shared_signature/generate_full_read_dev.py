"""Generate the prespecified small full-read M1 development challenge."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random

PROTOCOL = Path(__file__).with_name("full_read_protocol_v1.json")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(sequence: str, rng: random.Random, substitutions: float,
           indels: float = 0.0) -> str:
    output = []
    for base in sequence:
        if rng.random() < indels / 2:
            continue
        if rng.random() < substitutions:
            base = rng.choice([candidate for candidate in "ACGT" if candidate != base])
        output.append(base)
        if rng.random() < indels / 2:
            output.append(rng.choice("ACGT"))
    return "".join(output)


def generate(outdir: Path) -> dict:
    if outdir.exists():
        raise FileExistsError(outdir)
    protocol = json.loads(PROTOCOL.read_text())
    rng = random.Random(protocol["seed"])
    units = {name: random_dna(rng, protocol["unit_bp"])
             for name in ("f1", "f2", "twin_a", "decoy_zero")}
    units["twin_b"] = units["twin_a"]
    catalogue = outdir / "catalogue.fa"
    reads = outdir / "reads.fa"
    truth_path = outdir / "truth.json"
    outdir.mkdir(parents=True)
    catalogue.write_text("".join(f">{name}\n{units[name]}\n" for name in protocol["catalogue_families"]))
    truth = []
    fasta = []
    for scenario in protocol["scenarios"]:
        for replicate in range(protocol["read_repetitions_per_scenario"]):
            if scenario == "pure_f1":
                pieces = [("f1", units["f1"] * 5)]
            elif scenario == "mixed_f1_f2":
                pieces = [("background", random_dna(rng, 120)), ("f1", units["f1"] * 5),
                          ("f2", units["f2"] * 5), ("background", random_dna(rng, 120))]
            elif scenario == "reverse_f2":
                reverse = units["f2"].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                pieces = [("f2", reverse * 5)]
            elif scenario == "shared_twin":
                pieces = [("twin_a", units["twin_a"] * 5)]
            elif scenario == "background_only":
                pieces = [("background", random_dna(rng, 520))]
            else:
                pieces = [("background", random_dna(rng, 80)),
                          ("f1", units["f1"] * 5), ("f2", units["f2"] * 5),
                          ("background", random_dna(rng, 80))]
            sequence = ""
            spans = []
            for label, piece in pieces:
                altered = mutate(piece, rng, 0.03, 0.02 if scenario == "indel_mixed" else 0.0)
                start = len(sequence)
                sequence += altered
                spans.append(dict(start=start, end=len(sequence), label=label))
            read_id = f"{scenario}_{replicate}"
            truth.append(dict(read_id=read_id, length_bp=len(sequence), spans=spans))
            fasta.append(f">{read_id}\n{sequence}\n")
    reads.write_text("".join(fasta))
    truth_path.write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")
    manifest = dict(protocol_sha256=digest(PROTOCOL), files={
        name: digest(outdir / name) for name in ("catalogue.fa", "reads.fa", "truth.json")},
        read_count=len(truth), total_read_bp=sum(row["length_bp"] for row in truth))
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    generate(parser.parse_args().outdir)
