"""Generate independent, HiCAT-inspired synthetic HOR source lineages."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def substitute_fraction(sequence: str, fraction: float, rng: random.Random) -> str:
    if not 0 <= fraction < 1:
        raise ValueError("invalid substitution fraction")
    bases = list(sequence)
    for position in rng.sample(range(len(bases)), round(len(bases) * fraction)):
        bases[position] = rng.choice([base for base in "ACGT" if base != bases[position]])
    return "".join(bases)


def make_source(protocol_path: Path, out_manifest: Path) -> dict:
    protocol = json.loads(protocol_path.read_text())
    if (protocol.get("split") != "synthetic_held_out" or protocol.get("schema_version") != 2
            or protocol.get("read_pairing_status") != "synthetic_simulated"
            or protocol.get("monomer_lengths_bp") != [100, 200, 300, 400]):
        raise ValueError("B2 frozen source protocol gate failed")
    labels = protocol["monomer_labels"]
    path = list("".join(protocol["source_label_path_blocks"]))
    if len(path) != protocol["source_copy_count"] or any(label not in labels for label in path):
        raise ValueError("invalid synthetic nested HOR source path")
    lineages = []
    for length in protocol["monomer_lengths_bp"]:
        rng = random.Random(protocol["source_seed"] + length * 100003)
        founder = "".join(rng.choices("ACGT", k=length))
        templates = {label: substitute_fraction(founder, protocol["founder_divergence_fractions"][label], rng)
                     for label in labels}
        if len(set(templates.values())) != len(labels):
            raise ValueError("synthetic monomer templates collision")
        copies = []
        for index, label in enumerate(path):
            sequence = substitute_fraction(templates[label], protocol["copy_divergence_fraction"], rng)
            copies.append(dict(source_index=index, label=label, sequence=sequence,
                               sequence_sha256=sha256_text(sequence)))
        source_array = "".join(copy["sequence"] for copy in copies)
        lineages.append(dict(lineage_id=f"B2_L{length}", monomer_length_bp=length,
                             founder_sha256=sha256_text(founder),
                             candidate_monomers=templates, source_copies=copies,
                             source_label_path=path,
                             synthetic_canonical_HOR_period_monomers=protocol["canonical_HOR_period_monomers"],
                             nested_source_block=protocol["nested_source_block"],
                             source_array_sha256=sha256_text(source_array)))
    manifest = dict(schema_version=2, split="synthetic_held_out", source_type="independent_synthetic_HiCAT_inspired",
                    protocol_sha256=sha256_file(protocol_path), source_seed=protocol["source_seed"],
                    lineages=lineages, biological_validation_status="not_biological_validation")
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--out-manifest", type=Path, required=True)
    args = parser.parse_args()
    make_source(args.protocol, args.out_manifest)


if __name__ == "__main__":
    main()
