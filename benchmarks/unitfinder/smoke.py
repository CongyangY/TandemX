"""Generate a deterministic interface-only smoke input for unitFinder."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import random


SEED = 9101
DECOY_COPIES = (105, 110, 115, 120, 125, 130, 135, 140, 145, 150)
DECOY_PERIODS = (53, 59, 67, 73, 79, 83, 89, 97, 101, 109)
TARGET_PERIOD = 171
TARGET_COPIES = 600


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def wrapped(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate(outdir: Path) -> dict[str, object]:
    outdir.mkdir(parents=True, exist_ok=False)
    rng = random.Random(SEED)
    sequence_parts = [random_dna(rng, 25_000)]
    rows: list[dict[str, object]] = []
    position = len(sequence_parts[0])
    for index, (period, copies) in enumerate(zip(DECOY_PERIODS, DECOY_COPIES), start=1):
        monomer = random_dna(rng, period)
        array = monomer * copies
        start = position
        sequence_parts.append(array)
        position += len(array)
        rows.append({
            "array_id": f"decoy_{index:02d}",
            "start": start,
            "end": position,
            "period": period,
            "copies": copies,
            "role": "copy_number_background_for_unitfinder_outlier_rule",
        })
        spacer = random_dna(rng, 5_000)
        sequence_parts.append(spacer)
        position += len(spacer)
    target_monomer = random_dna(rng, TARGET_PERIOD)
    target_start = position
    target_array = target_monomer * TARGET_COPIES
    sequence_parts.append(target_array)
    position += len(target_array)
    target_end = position
    rows.append({
        "array_id": "target_high_copy",
        "start": target_start,
        "end": target_end,
        "period": TARGET_PERIOD,
        "copies": TARGET_COPIES,
        "role": "expected_unitfinder_copy_number_outlier",
    })
    sequence_parts.append(random_dna(rng, 25_000))
    sequence = "".join(sequence_parts)
    fasta = outdir / "unitfinder_smoke.fa"
    fasta.write_text(">unitfinder_smoke_chr\n" + wrapped(sequence) + "\n", encoding="ascii")
    truth = outdir / "truth_arrays.tsv"
    with truth.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("array_id", "start", "end", "period", "copies", "role"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": 1,
        "scope": "unitfinder_installation_and_interface_smoke_not_accuracy_evidence",
        "seed": SEED,
        "coordinate_system": "0_based_half_open",
        "sequence_length_bp": len(sequence),
        "decoy_copy_counts": list(DECOY_COPIES),
        "target": {
            "array_id": "target_high_copy",
            "start": target_start,
            "end": target_end,
            "period": TARGET_PERIOD,
            "copies": TARGET_COPIES,
        },
        "files": {
            fasta.name: {"bytes": fasta.stat().st_size, "sha256": digest(fasta)},
            truth.name: {"bytes": truth.stat().st_size, "sha256": digest(truth)},
        },
        "warning": "successful_recovery_does_not_establish_real_data_or_general_accuracy",
    }
    (outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest
