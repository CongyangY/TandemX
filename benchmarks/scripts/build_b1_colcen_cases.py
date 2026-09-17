"""Build frozen B1 v3 edited Col-CEN interval cases from source-qualified fragments."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path

from benchmarks.scripts.build_b1_structure_cases import (
    make_flanks, reverse_complement, sequence_hash, sha256_file, transform, write_jsonl,
)


def load_source(protocol_path: Path, source_fasta: Path, catalogue_path: Path) -> tuple[dict, dict[str, str], dict[str, str]]:
    protocol = json.loads(protocol_path.read_text())
    catalogue = json.loads(catalogue_path.read_text())
    if (protocol.get("schema_version") != 2 or protocol.get("split") != "development"
            or protocol.get("read_pairing_status") != "synthetic_simulated"
            or protocol.get("held_out_status") != "not_enrolled"):
        raise ValueError("v3 development/pairing gate failed")
    if (catalogue["source_assembly_sha256"] != protocol["source_assembly_sha256"]
            or catalogue["source_fasta_sha256"] != sha256_file(source_fasta)):
        raise ValueError("Col-CEN source catalogue/FASTA hash mismatch")
    arrays = {}
    header = None
    for line in source_fasta.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            header = line[1:].split("|")[0]
            if header in arrays:
                raise ValueError("duplicate source array ID")
            arrays[header] = ""
        elif header is None:
            raise ValueError("source FASTA sequence before header")
        else:
            arrays[header] += line.upper()
    expected = {row["array_id"] for row in protocol["source_intervals"]}
    if set(arrays) != expected:
        raise ValueError("source array IDs mismatch")
    period = protocol["operational_period_bp"]
    count = protocol["operational_tile_count"]
    if period != 178 or count != 20:
        raise ValueError("frozen 178bp/20-tile qualification changed")
    consensus = {}
    for row in protocol["source_intervals"]:
        array_id = row["array_id"]
        sequence = arrays[array_id]
        if (len(sequence) != period * count or set(sequence) - set("ACGT")
                or sequence_hash(sequence) != row["sequence_sha256"]):
            raise ValueError(f"source array sequence changed: {array_id}")
        tiles = [sequence[i:i + period] for i in range(0, len(sequence), period)]
        adjacent = [sum(a == b for a, b in zip(left, right)) / period
                    for left, right in zip(tiles, tiles[1:])]
        if min(adjacent) < protocol["minimum_adjacent_tile_identity"]:
            raise ValueError(f"native source periodicity gate failed: {array_id}")
        consensus[array_id] = "".join(Counter(tile[i] for tile in tiles).most_common(1)[0][0]
                                       for i in range(period))
    if len(set(consensus.values())) != len(consensus):
        raise ValueError("catalogue consensuses unexpectedly identical")
    return protocol, arrays, consensus


def make_case(protocol: dict, source_row: dict, source_array: str, candidate_monomers: dict[str, str],
              left: str, right: str, spec: dict) -> tuple[dict, dict]:
    array_id = source_row["array_id"]
    period = protocol["operational_period_bp"]
    tiles = [source_array[i:i + period] for i in range(0, len(source_array), period)]
    base = [dict(source_index=i, label=array_id, orientation="+", sequence=tile,
                 operation="retain", duplicate=False) for i, tile in enumerate(tiles)]
    edited = transform(base, spec)
    edited_array = "".join(copy["sequence"] for copy in edited)
    source_full = left + source_array + right
    edited_full = left + edited_array + right
    if any(sequence.count(left) != 1 or sequence.count(right) != 1 for sequence in (source_full, edited_full)):
        raise ValueError("engineered flank not unique")
    case_id = "B1R" + hashlib.sha256((protocol["case_id_salt"] + "/" + array_id + "/" +
                                      spec["name"]).encode()).hexdigest()[:12]
    flank = len(left)
    deleted_indices = sorted(set(range(len(base))) - {copy["source_index"] for copy in edited})
    deleted_spans = [[flank + i * period, flank + (i + 1) * period] for i in deleted_indices]
    if spec["operation"] == "trim_tail_each":
        trim = spec["trim_bp"]
        deleted_spans = [[flank + (i + 1) * period - trim, flank + (i + 1) * period]
                         for i in range(len(base))]
    chain = [dict(operation="retain_flank", source_start0=0, source_end0=flank,
                  edited_start0=0, edited_end0=flank, orientation="+")]
    insertions = []
    position = flank
    for edited_index, copy in enumerate(edited):
        start = position
        position += len(copy["sequence"])
        segment = dict(operation=copy["operation"], source_start0=flank + period * copy["source_index"],
                       source_end0=flank + period * (copy["source_index"] + 1),
                       edited_start0=start, edited_end0=position,
                       orientation=copy["orientation"], source_copy_index=copy["source_index"],
                       edited_copy_index=edited_index)
        chain.append(segment)
        if copy["duplicate"]:
            insertions.append([start, position])
    for start, end in deleted_spans:
        source_index = (end - 1 - flank) // period
        junction = flank + sum(len(copy["sequence"]) for copy in edited
                               if copy["source_index"] < source_index and not copy["duplicate"])
        if spec["operation"] == "trim_tail_each":
            junction += len(edited[source_index]["sequence"])
        chain.append(dict(operation="delete", source_start0=start, source_end0=end,
                          edited_start0=junction, edited_end0=junction, orientation="+"))
    chain.append(dict(operation="retain_flank", source_start0=flank + len(source_array),
                      source_end0=len(source_full), edited_start0=position,
                      edited_end0=len(edited_full), orientation="+"))
    removed_bp = sum(end - start for start, end in deleted_spans)
    inserted_bp = sum(end - start for start, end in insertions)
    if len(edited_full) - len(source_full) != inserted_bp - removed_bp:
        raise ValueError("edit bp conservation failed")
    reads = {f"simulated_read_{i:02d}": source_full for i in range(1, protocol["read_count"] + 1)}
    common = dict(schema_version=2, case_id=case_id, split="development",
                  material_id="Col_CEN_v1.2_assembly_lineage", array_id=array_id,
                  source_type="native_assembly_interval_with_synthetic_edits_and_reads",
                  donor_id=None, assembly_lineage_id="Col_CEN_v1.2",
                  haplotype_id=None, homology_group_id="Col_CEN_CEN180_related",
                  source_chromosome=source_row["chromosome"],
                  source_assembly_interval_bp=[source_row["start0"], source_row["end0"]])
    public = dict(common, assembly_sequence=edited_full, raw_read_sequences=reads,
                  candidate_monomers=candidate_monomers,
                  read_pairing_status="synthetic_simulated", error_profile=protocol["read_error_profile"],
                  haplotype_profile=protocol["haplotype_profile"], input_status="ok",
                  left_flank_sequence=left, right_flank_sequence=right,
                  array_window_policy="between_exact_unique_synthetic_flanks",
                  flank_source="engineered_not_native",
                  source_tile_status=protocol["tile_status"],
                  assembly_sequence_sha256=sequence_hash(edited_full),
                  raw_reads_sha256=hashlib.sha256(json.dumps(reads, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                  candidate_monomers_sha256=hashlib.sha256(json.dumps(candidate_monomers, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                  left_flank_sha256=sequence_hash(left), right_flank_sha256=sequence_hash(right))
    truth = dict(common, event_status="intact_negative" if spec["event_type"] == "intact" else "positive",
                 event_type=spec["event_type"], truth_scope="injected_edit_only",
                 source_label_path=[array_id] * len(base), edited_label_path=[array_id] * len(edited),
                 source_orientation_path=["+"] * len(base),
                 edited_orientation_path=[copy["orientation"] for copy in edited],
                 source_tile_id_path=[f"{array_id}:{i}" for i in range(len(base))],
                 edited_tile_id_path=[f"{array_id}:{copy['source_index']}" for copy in edited],
                 source_array_interval_bp=[flank, flank + len(source_array)],
                 edited_array_interval_bp=[flank, flank + len(edited_array)],
                 source_deleted_spans_bp=deleted_spans, edited_inserted_spans_bp=insertions,
                 removed_bp=removed_bp, inserted_bp=inserted_bp,
                 signed_bp_delta=inserted_bp - removed_bp,
                 source_copy_count=len(base), edited_copy_count=len(edited),
                 unit_truth_status="operational_178bp_tiles_not_native_monomer_truth",
                 hor_period_status="unverified", coordinate_chain=chain,
                 source_sequence_sha256=sequence_hash(source_full),
                 edited_sequence_sha256=sequence_hash(edited_full))
    return public, truth


def build(protocol_path: Path, source_fasta: Path, catalogue_path: Path, outdir: Path) -> dict:
    protocol, arrays, monomers = load_source(protocol_path, source_fasta, catalogue_path)
    pairs = []
    for source_row in protocol["source_intervals"]:
        array_id = source_row["array_id"]
        flank_protocol = dict(protocol, flank_seed=protocol["flank_seed"] + int(array_id[1:]))
        left, right = make_flanks(flank_protocol, arrays[array_id], monomers)
        for spec in protocol["case_specs"]:
            pairs.append(make_case(protocol, source_row, arrays[array_id], monomers,
                                   left, right, spec))
    if len(pairs) != 39 or len({public["case_id"] for public, _ in pairs}) != 39:
        raise ValueError("v3 case count/IDs failed")
    random.Random(protocol["case_order_seed"]).shuffle(pairs)
    outdir.mkdir(parents=True, exist_ok=False)
    inputs = outdir / "inputs.jsonl"
    truths = outdir / "truth.jsonl"
    write_jsonl(inputs, [public for public, _ in pairs])
    write_jsonl(truths, [truth for _, truth in pairs])
    assemblies = outdir / "assemblies"
    assemblies.mkdir()
    fasta_hashes = {}
    for public, _ in pairs:
        path = assemblies / f"{public['case_id']}.fa"
        path.write_text(f">{public['case_id']} source=Col_CEN_v1.2_development\n{public['assembly_sequence']}\n", encoding="ascii")
        fasta_hashes[path.name] = sha256_file(path)
    receipt = dict(schema_version=2, development_only=True,
                   source_assembly_sha256=protocol["source_assembly_sha256"],
                   source_fasta_sha256=sha256_file(source_fasta),
                   source_catalogue_sha256=sha256_file(catalogue_path),
                   protocol_sha256=sha256_file(protocol_path),
                   generator_sha256=sha256_file(Path(__file__)),
                   inputs_sha256=sha256_file(inputs), truth_sha256=sha256_file(truths),
                   assembly_fasta_sha256=fasta_hashes,
                   case_count=len(pairs), case_ids=[public["case_id"] for public, _ in pairs],
                   source_array_count=len(arrays), independent_donor_count=1,
                   split="development", held_out_status="not_enrolled",
                   read_pairing_status="synthetic_simulated",
                   original_raw_read_pairing_status="unverified",
                   biological_accuracy_status=protocol["biological_accuracy_status"],
                   masking_status="separate_truth_file_not_secure_blinding")
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--source-fasta", type=Path, required=True)
    parser.add_argument("--source-catalogue", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    build(args.protocol, args.source_fasta, args.source_catalogue, args.outdir)


if __name__ == "__main__":
    main()
