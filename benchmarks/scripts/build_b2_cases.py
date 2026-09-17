"""Build masked synthetic structural hold-out cases from independent B2 source lineages."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

from benchmarks.scripts.build_b1_structure_cases import make_flanks, sequence_hash, sha256_file, transform, write_jsonl


def mutate_read_array(source: str, rng: random.Random, substitution_rate: float,
                      indel_rate: float, source_offset0: int) -> tuple[str, list[dict]]:
    output = []
    events = []
    for position, base in enumerate(source):
        absolute = source_offset0 + position
        if rng.random() < indel_rate:
            if rng.random() < 0.5:
                events.append(dict(operation="delete", source_position0=absolute, source_base=base))
                continue
            inserted = rng.choice("ACGT")
            output.append(inserted)
            events.append(dict(operation="insert_before", source_position0=absolute, inserted_base=inserted))
        if rng.random() < substitution_rate:
            replacement = rng.choice([candidate for candidate in "ACGT" if candidate != base])
            output.append(replacement)
            events.append(dict(operation="substitute", source_position0=absolute,
                               source_base=base, replacement_base=replacement))
        else:
            output.append(base)
    return "".join(output), events


def load_source(protocol_path: Path, manifest_path: Path) -> tuple[dict, list[dict]]:
    protocol = json.loads(protocol_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if (protocol.get("schema_version") != 2 or protocol.get("split") != "synthetic_held_out"
            or manifest.get("protocol_sha256") != sha256_file(protocol_path)
            or manifest.get("split") != "synthetic_held_out"):
        raise ValueError("B2 frozen source/protocol gate failed")
    lineages = manifest["lineages"]
    if len(lineages) != 4 or {lineage["monomer_length_bp"] for lineage in lineages} != set(protocol["monomer_lengths_bp"]):
        raise ValueError("B2 source lineage set mismatch")
    for lineage in lineages:
        length = lineage["monomer_length_bp"]
        copies = lineage["source_copies"]
        if (len(copies) != protocol["source_copy_count"] or
                [copy["label"] for copy in copies] != lineage["source_label_path"] or
                any(len(copy["sequence"]) != length for copy in copies) or
                any(sequence_hash(copy["sequence"]) != copy["sequence_sha256"] for copy in copies) or
                sequence_hash("".join(copy["sequence"] for copy in copies)) != lineage["source_array_sha256"]):
            raise ValueError("B2 source copy integrity failed")
    return protocol, lineages


def make_case(protocol: dict, lineage: dict, spec: dict, left: str, right: str,
              case_index: int) -> tuple[dict, dict]:
    length = lineage["monomer_length_bp"]
    base = [dict(copy, orientation="+", operation="retain", duplicate=False)
            for copy in lineage["source_copies"]]
    source_array = "".join(copy["sequence"] for copy in base)
    edited = transform(base, spec)
    edited_array = "".join(copy["sequence"] for copy in edited)
    source_full = left + source_array + right
    edited_full = left + edited_array + right
    case_id = "B2" + hashlib.sha256((protocol["case_id_salt"] + "/" + lineage["lineage_id"] +
                                     "/" + spec["name"]).encode()).hexdigest()[:13]
    flank = len(left)
    deleted = sorted(set(range(len(base))) - {copy["source_index"] for copy in edited})
    deleted_spans = [[flank + i * length, flank + (i + 1) * length] for i in deleted]
    if spec["operation"] == "trim_tail_each":
        trim = spec["trim_bp"]
        deleted_spans = [[flank + (i + 1) * length - trim, flank + (i + 1) * length]
                         for i in range(len(base))]
    chain = [dict(operation="retain_flank", source_start0=0, source_end0=flank,
                  edited_start0=0, edited_end0=flank, orientation="+")]
    inserted_spans = []
    position = flank
    for edited_index, copy in enumerate(edited):
        start = position
        position += len(copy["sequence"])
        chain.append(dict(operation=copy["operation"], source_start0=flank + length * copy["source_index"],
                          source_end0=flank + length * (copy["source_index"] + 1),
                          edited_start0=start, edited_end0=position,
                          orientation=copy["orientation"], source_copy_index=copy["source_index"],
                          edited_copy_index=edited_index))
        if copy["duplicate"]:
            inserted_spans.append([start, position])
    for start, end in deleted_spans:
        source_index = (end - 1 - flank) // length
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
    inserted_bp = sum(end - start for start, end in inserted_spans)
    if len(edited_full) - len(source_full) != inserted_bp - removed_bp:
        raise ValueError("B2 edit length conservation failed")
    noisy = spec["error_mode"] == "noisy"
    sub_rate = protocol["noisy_substitution_rate"] if noisy else 0.0
    indel_rate = protocol["noisy_indel_rate"] if noisy else 0.0
    reads = {}
    error_ledger = {}
    for read_index in range(spec["read_support"]):
        rng = random.Random(protocol["read_error_seed"] + length * 100003 + case_index * 1009 + read_index)
        mutated, errors = mutate_read_array(source_array, rng, sub_rate, indel_rate, flank)
        read_id = f"simulated_read_{read_index + 1:02d}"
        reads[read_id] = left + mutated + right
        error_ledger[read_id] = errors
    if any(read.count(left) != 1 or read.count(right) != 1 for read in (edited_full, *reads.values())):
        raise ValueError("B2 synthetic flank not unique")
    common = dict(schema_version=2, case_id=case_id, split="synthetic_held_out",
                  material_id=lineage["lineage_id"], source_type="synthetic_HiCAT_inspired",
                  donor_id=None, assembly_lineage_id=lineage["lineage_id"],
                  haplotype_id=None, homology_group_id=lineage["lineage_id"],
                  monomer_length_bp=length)
    public = dict(common, assembly_sequence=edited_full, raw_read_sequences=reads,
                  candidate_monomers=lineage["candidate_monomers"],
                  read_pairing_status="synthetic_simulated", input_status="ok",
                  error_profile=dict(name="array_only_synthetic", substitution_rate=sub_rate,
                                     indel_rate=indel_rate, engineered_flanks_protected=True),
                  haplotype_profile=dict(name="one_synthetic_source_path", phased=True),
                  left_flank_sequence=left, right_flank_sequence=right,
                  array_window_policy="between_exact_unique_synthetic_flanks",
                  flank_source="engineered_not_native",
                  synthetic_read_support=len(reads),
                  assembly_sequence_sha256=sequence_hash(edited_full),
                  raw_reads_sha256=hashlib.sha256(json.dumps(reads, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                  candidate_monomers_sha256=hashlib.sha256(json.dumps(lineage["candidate_monomers"], sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                  left_flank_sha256=sequence_hash(left), right_flank_sha256=sequence_hash(right))
    truth = dict(common, event_status="intact_negative" if spec["event_type"] == "intact" else "positive",
                 event_type=spec["event_type"], truth_scope="exact_injected_synthetic_edits_only",
                 source_label_path=[copy["label"] for copy in base],
                 edited_label_path=[copy["label"] for copy in edited],
                 source_orientation_path=["+"] * len(base),
                 edited_orientation_path=[copy["orientation"] for copy in edited],
                 source_array_interval_bp=[flank, flank + len(source_array)],
                 edited_array_interval_bp=[flank, flank + len(edited_array)],
                 source_deleted_spans_bp=deleted_spans, edited_inserted_spans_bp=inserted_spans,
                 removed_bp=removed_bp, inserted_bp=inserted_bp,
                 signed_bp_delta=inserted_bp - removed_bp,
                 source_copy_count=len(base), edited_copy_count=len(edited),
                 source_HOR_period_monomers=protocol["canonical_HOR_period_monomers"],
                 source_nested_block=protocol["nested_source_block"],
                 unit_truth_status="exact_synthetic", coordinate_chain=chain,
                 read_error_events=error_ledger,
                 source_sequence_sha256=sequence_hash(source_full),
                 edited_sequence_sha256=sequence_hash(edited_full))
    return public, truth


def build(protocol_path: Path, manifest_path: Path, outdir: Path) -> dict:
    protocol, lineages = load_source(protocol_path, manifest_path)
    pairs = []
    for lineage in lineages:
        source_array = "".join(copy["sequence"] for copy in lineage["source_copies"])
        flank_protocol = dict(protocol, flank_seed=protocol["flank_seed"] + lineage["monomer_length_bp"])
        left, right = make_flanks(flank_protocol, source_array, lineage["candidate_monomers"])
        for index, spec in enumerate(protocol["case_specs"]):
            pairs.append(make_case(protocol, lineage, spec, left, right, index))
    if len(pairs) != 52 or len({public["case_id"] for public, _ in pairs}) != 52:
        raise ValueError("B2 expected 52 distinct case IDs")
    random.Random(protocol["case_seed"]).shuffle(pairs)
    outdir.mkdir(parents=True, exist_ok=False)
    inputs = outdir / "inputs.jsonl"
    truth = outdir / "truth.jsonl"
    write_jsonl(inputs, [public for public, _ in pairs])
    write_jsonl(truth, [target for _, target in pairs])
    assemblies = outdir / "assemblies"
    assemblies.mkdir()
    fasta_hashes = {}
    for public, _ in pairs:
        path = assemblies / f"{public['case_id']}.fa"
        path.write_text(f">{public['case_id']} source=synthetic_B2_held_out\n{public['assembly_sequence']}\n", encoding="ascii")
        fasta_hashes[path.name] = sha256_file(path)
    receipt = dict(schema_version=2, split="synthetic_held_out", development_only=False,
                   source_type="independent_synthetic_HiCAT_inspired",
                   protocol_sha256=sha256_file(protocol_path),
                   source_manifest_sha256=sha256_file(manifest_path),
                   source_generator_sha256=sha256_file(Path(__file__).with_name("generate_b2_source.py")),
                   case_generator_sha256=sha256_file(Path(__file__)),
                   inputs_sha256=sha256_file(inputs), truth_sha256=sha256_file(truth),
                   assembly_fasta_sha256=fasta_hashes, case_count=len(pairs),
                   case_ids=[public["case_id"] for public, _ in pairs],
                   synthetic_founder_count=len(lineages), biological_donor_count=0,
                   read_pairing_status="synthetic_simulated",
                   biological_validation_status="not_biological_validation",
                   masking_status="separate_truth_file_not_secure_blinding")
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    build(args.protocol, args.source_manifest, args.outdir)


if __name__ == "__main__":
    main()
