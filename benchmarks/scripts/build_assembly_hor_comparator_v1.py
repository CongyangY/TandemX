"""Build a fixed, sequence-level 171-bp HOR comparator with exact edit truth.

This is an engineered control inspired by the parameter range of the HiCAT
simulation, not a sampled biological centromere or held-out validation set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path


SEED = 17092026
MONOMER_BP = 171
FAMILY_COUNT = 5
HOR_COUNT = 40
FLANK_BP = 1000
BASES = "ACGT"
VARIANTS = {7: ("A", "B", "C", "E"), 18: ("A", "B", "C", "C", "D", "E"),
            29: ("A", "B", "D", "C", "E"), 35: ("A", "B", "C", "D", "E", "E")}


def mutate(sequence: str, count: int, rng: random.Random) -> str:
    value = list(sequence)
    for index in rng.sample(range(len(value)), count):
        value[index] = rng.choice([base for base in BASES if base != value[index]])
    return "".join(value)


def build(outdir: Path) -> dict[str, object]:
    rng = random.Random(SEED)
    ancestor = "".join(rng.choices(BASES, k=MONOMER_BP))
    monomers = {label: mutate(ancestor, 34, rng) for label in "ABCDE"}
    left = "".join(rng.choices(BASES, k=FLANK_BP))
    right = "".join(rng.choices(BASES, k=FLANK_BP))
    labels = []
    copies = []
    hors = []
    cursor = FLANK_BP
    for hor_index in range(HOR_COUNT):
        path = VARIANTS.get(hor_index, tuple("ABCDE"))
        start = cursor
        for within_index, label in enumerate(path):
            sequence = mutate(monomers[label], 1, rng)
            copies.append((hor_index, within_index, label, cursor, cursor + MONOMER_BP, sequence))
            labels.append(label)
            cursor += MONOMER_BP
        hors.append((hor_index, start, cursor, "".join(path), path == tuple("ABCDE")))
    assembly = left + "".join(row[5] for row in copies) + right
    outdir.mkdir(parents=True, exist_ok=False)
    fasta = outdir / "control.chr1.fasta"
    fasta.write_text(">control.chr1\n" + assembly + "\n")
    (outdir / "candidate_monomers.fa").write_text(
        "".join(f">{label}\n{sequence}\n" for label, sequence in monomers.items()))
    with (outdir / "truth_monomers.tsv").open("w") as handle:
        handle.write("sequence_id\tstart0\tend0\tlabel\thor_index\twithin_hor_index\torientation\n")
        for hor_index, within_index, label, start, end, _ in copies:
            handle.write(f"control.chr1\t{start}\t{end}\t{label}\t{hor_index}\t{within_index}\t+\n")
    with (outdir / "truth_hors.tsv").open("w") as handle:
        handle.write("sequence_id\tstart0\tend0\thor_index\tlabel_path\tcanonical\n")
        for index, start, end, path, canonical in hors:
            handle.write(f"control.chr1\t{start}\t{end}\t{index}\t{path}\t{int(canonical)}\n")
    receipt = {
        "schema_version": 1, "seed": SEED, "truth_type": "exact_engineered_only",
        "biological_source_status": "synthetic_centromere_like_not_biological_truth",
        "monomer_bp": MONOMER_BP, "family_count": FAMILY_COUNT,
        "hor_count": HOR_COUNT, "variant_hor_indices": sorted(VARIANTS),
        "array_interval_0based_halfopen": [FLANK_BP, cursor],
        "assembly_bp": len(assembly), "copy_count": len(copies),
        "input_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        "candidate_monomer_sha256": hashlib.sha256((outdir / "candidate_monomers.fa").read_bytes()).hexdigest(),
        "monomer_labels": labels,
    }
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.outdir), sort_keys=True))


if __name__ == "__main__":
    main()
