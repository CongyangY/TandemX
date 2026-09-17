"""Generate masked, real-sequence-derived B1 structure development cases.

The archived 785-bp YSD56 candidate contributes short sequence tiles only.
All array copies, HOR labels, flanks and reads in this bundle are engineered.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from pathlib import Path


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def sequence_hash(sequence: str) -> str:
    return sha256_bytes(sequence.encode("ascii"))


def object_hash(value: object) -> str:
    return sha256_bytes(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def load_protocol(path: Path, root: Path) -> tuple[dict, dict[str, str]]:
    protocol = json.loads(path.read_text())
    if (protocol.get("schema_version") != 2 or protocol.get("split") != "development"
            or protocol.get("read_pairing_status") != "synthetic_simulated"
            or protocol.get("held_out_status") != "not_enrolled"):
        raise ValueError("B1 v2 protocol split/pairing gate failed")
    source = root / protocol["source_template_path"]
    if sha256_file(source) != protocol["source_template_sha256"]:
        raise ValueError("archived source-template hash mismatch")
    lines = source.read_text(encoding="ascii").splitlines()
    if not lines or not lines[0].startswith(">") or any(line.startswith(">") for line in lines[1:]):
        raise ValueError("source template must be one FASTA record")
    template = "".join(lines[1:]).upper()
    if len(template) != 785 or set(template) - set("ACGT"):
        raise ValueError("archived 785-bp candidate template changed")
    size = protocol["monomer_tile_length_bp"]
    offsets = protocol["monomer_tile_offsets"]
    if not isinstance(size, int) or not 1 <= size <= 64 or set(offsets) != set(protocol["engineered_hor_unit_labels"]):
        raise ValueError("invalid engineered monomer tile plan")
    motifs = {}
    for label, offset in offsets.items():
        if not isinstance(offset, int) or offset < 0 or offset + size > len(template):
            raise ValueError("monomer tile outside archived template")
        motifs[label] = template[offset:offset + size]
    if len(set(motifs.values())) != len(motifs):
        raise ValueError("engineered monomer tiles are duplicate")
    return protocol, motifs


def make_flanks(protocol: dict, source_array: str, motifs: dict[str, str]) -> tuple[str, str]:
    rng = random.Random(protocol["flank_seed"])
    n = protocol["engineered_flank_length_bp"]
    if not isinstance(n, int) or not 24 <= n <= 256:
        raise ValueError("engineered flank length outside bounded range")
    flanks = []
    for _ in range(2):
        for _attempt in range(1000):
            candidate = "".join(rng.choices("ACGT", k=n))
            if candidate not in source_array and all(candidate != previous for previous in flanks) and all(
                candidate not in motif for motif in motifs.values()
            ):
                flanks.append(candidate)
                break
        else:
            raise ValueError("could not construct unique synthetic flanks")
    left, right = flanks
    full = left + source_array + right
    if full.count(left) != 1 or full.count(right) != 1:
        raise ValueError("engineered flank anchors are not unique")
    return left, right


def source_copies(protocol: dict, motifs: dict[str, str]) -> list[dict]:
    unit = protocol["engineered_hor_unit_labels"]
    count = protocol["engineered_hor_copies"]
    if not isinstance(count, int) or count < 2 or count > 12:
        raise ValueError("HOR copy count outside bounded range")
    labels = unit * count
    return [dict(source_index=i, label=label, orientation="+", sequence=motifs[label],
                 operation="retain", duplicate=False) for i, label in enumerate(labels)]


def transform(base: list[dict], spec: dict) -> list[dict]:
    operation = spec["operation"]
    indices = spec.get("indices", [])
    if any(not isinstance(i, int) or i < 0 or i >= len(base) for i in indices) or len(set(indices)) != len(indices):
        raise ValueError(f"invalid source copy indices: {spec['name']}")
    copies = [dict(copy) for copy in base]
    if operation == "identity":
        return copies
    if operation == "delete_indices":
        return [copy for copy in copies if copy["source_index"] not in indices]
    if operation == "duplicate_indices_after":
        if not indices or indices != list(range(indices[0], indices[-1] + 1)):
            raise ValueError("duplication requires one contiguous source block")
        duplicate = [dict(base[i], operation="duplicate", duplicate=True) for i in indices]
        return copies[:indices[-1] + 1] + duplicate + copies[indices[-1] + 1:]
    if operation == "swap_indices":
        if len(indices) != 2:
            raise ValueError("swap requires two indices")
        copies[indices[0]], copies[indices[1]] = copies[indices[1]], copies[indices[0]]
        for i in indices:
            copies[i]["operation"] = "reorder"
        return copies
    if operation == "invert_indices":
        for i in indices:
            copies[i]["sequence"] = reverse_complement(copies[i]["sequence"])
            copies[i]["orientation"] = "-"
            copies[i]["operation"] = "reverse_complement"
        return copies
    if operation == "trim_tail_each":
        trim = spec["trim_bp"]
        if not isinstance(trim, int) or trim < 1 or trim >= len(base[0]["sequence"]):
            raise ValueError("invalid per-copy compression")
        for copy in copies:
            copy["sequence"] = copy["sequence"][:-trim]
            copy["operation"] = "compress_tail"
        return copies
    raise ValueError(f"unsupported edit operation: {operation}")


def make_case(protocol: dict, spec: dict, base: list[dict], motifs: dict[str, str],
              left: str, right: str) -> tuple[dict, dict]:
    case_id = "B1" + sha256_bytes((protocol["case_id_salt"] + "/" + spec["name"]).encode())[:12]
    source_array = "".join(copy["sequence"] for copy in base)
    source_full = left + source_array + right
    edited_copies = transform(base, spec)
    edited_array = "".join(copy["sequence"] for copy in edited_copies)
    edited_full = left + edited_array + right
    left_bp = len(left)
    size = len(base[0]["sequence"])
    source_indices = {copy["source_index"] for copy in edited_copies}
    deleted_indices = sorted(set(range(len(base))) - source_indices)
    deleted_spans = [[left_bp + i * size, left_bp + (i + 1) * size] for i in deleted_indices]
    if spec["operation"] == "trim_tail_each":
        trim = spec["trim_bp"]
        deleted_spans = [[left_bp + (i + 1) * size - trim, left_bp + (i + 1) * size]
                         for i in range(len(base))]
    inserted_spans = []
    chain = [dict(operation="retain_flank", source_start0=0, source_end0=left_bp,
                  edited_start0=0, edited_end0=left_bp, orientation="+")]
    position = left_bp
    for edited_index, copy in enumerate(edited_copies):
        start = position
        position += len(copy["sequence"])
        source_start = left_bp + copy["source_index"] * size
        chain.append(dict(operation=copy["operation"], source_start0=source_start,
                          source_end0=source_start + size, edited_start0=start,
                          edited_end0=position, orientation=copy["orientation"],
                          source_copy_index=copy["source_index"], edited_copy_index=edited_index))
        if copy["duplicate"]:
            inserted_spans.append([start, position])
    for start, end in deleted_spans:
        chain.append(dict(operation="delete", source_start0=start, source_end0=end,
                          edited_start0=position, edited_end0=position, orientation="+"))
    chain.append(dict(operation="retain_flank", source_start0=left_bp + len(source_array),
                      source_end0=len(source_full), edited_start0=position,
                      edited_end0=len(edited_full), orientation="+"))
    removed_bp = sum(end - start for start, end in deleted_spans)
    inserted_bp = sum(end - start for start, end in inserted_spans)
    if len(edited_full) - len(source_full) != inserted_bp - removed_bp:
        raise ValueError(f"edit length conservation failed: {spec['name']}")
    reads = {f"simulated_read_{i:02d}": source_full for i in range(1, protocol["read_count"] + 1)}
    if any(sequence.count(left) != 1 or sequence.count(right) != 1 for sequence in (edited_full, *reads.values())):
        raise ValueError("synthetic flank anchor not unique in case")
    common = dict(schema_version=2, case_id=case_id, material_id=protocol["material_id"],
                  split="development", source_type="semi_synthetic_real_sequence_derived",
                  donor_id=None, assembly_lineage_id=protocol["assembly_lineage_id"],
                  haplotype_id=None, homology_group_id=protocol["homology_group_id"])
    public = dict(common, assembly_sequence=edited_full, raw_read_sequences=reads,
                  candidate_monomers=motifs, read_pairing_status="synthetic_simulated",
                  error_profile=protocol["read_error_profile"],
                  haplotype_profile=protocol["haplotype_profile"], input_status="ok",
                  left_flank_sequence=left, right_flank_sequence=right,
                  array_window_policy="between_exact_unique_synthetic_flanks",
                  flank_source="engineered_not_native",
                  assembly_sequence_sha256=sequence_hash(edited_full),
                  raw_reads_sha256=object_hash(reads), candidate_monomers_sha256=object_hash(motifs),
                  left_flank_sha256=sequence_hash(left), right_flank_sha256=sequence_hash(right))
    edited_labels = [copy["label"] for copy in edited_copies]
    edited_orientations = [copy["orientation"] for copy in edited_copies]
    hor_unit = protocol["engineered_hor_unit_labels"]
    edited_hor_count = (len(edited_labels) // len(hor_unit) if edited_labels == hor_unit *
                        (len(edited_labels) // len(hor_unit)) and all(x == "+" for x in edited_orientations)
                        and spec["operation"] != "trim_tail_each" else None)
    truth = dict(common, event_status="intact_negative" if spec["event_type"] == "intact" else "positive",
                 event_type=spec["event_type"], source_label_path=[copy["label"] for copy in base],
                 edited_label_path=edited_labels, source_orientation_path=["+"] * len(base),
                 edited_orientation_path=edited_orientations,
                 source_array_interval_bp=[left_bp, left_bp + len(source_array)],
                 edited_array_interval_bp=[left_bp, left_bp + len(edited_array)],
                 source_deleted_spans_bp=deleted_spans, edited_inserted_spans_bp=inserted_spans,
                 removed_bp=removed_bp, inserted_bp=inserted_bp,
                 signed_bp_delta=inserted_bp - removed_bp,
                 source_hor_unit_labels=hor_unit, source_hor_count=protocol["engineered_hor_copies"],
                 edited_hor_count=edited_hor_count,
                 source_copy_count=len(base), edited_copy_count=len(edited_copies),
                 unit_truth_status="exact_engineered", truth_scope="injected_edit_only",
                 coordinate_chain=chain, source_sequence_sha256=sequence_hash(source_full),
                 edited_sequence_sha256=sequence_hash(edited_full))
    return public, truth


def write_jsonl(path: Path, rows: list[dict]) -> None:
    path.write_text("".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n"
                            for row in rows), encoding="utf-8")


def build(protocol_path: Path, root: Path, outdir: Path) -> None:
    protocol, motifs = load_protocol(protocol_path, root)
    base = source_copies(protocol, motifs)
    source_array = "".join(copy["sequence"] for copy in base)
    left, right = make_flanks(protocol, source_array, motifs)
    pairs = [make_case(protocol, spec, base, motifs, left, right) for spec in protocol["case_specs"]]
    if len({public["case_id"] for public, _ in pairs}) != len(pairs):
        raise ValueError("case ID collision")
    random.Random(protocol["case_order_seed"]).shuffle(pairs)
    outdir.mkdir(parents=True, exist_ok=False)
    public_rows = [public for public, _ in pairs]
    truth_rows = [truth for _, truth in pairs]
    inputs = outdir / "inputs.jsonl"
    truths = outdir / "truth.jsonl"
    write_jsonl(inputs, public_rows)
    write_jsonl(truths, truth_rows)
    fasta_dir = outdir / "assemblies"
    fasta_dir.mkdir()
    assembly_hashes = {}
    for public in public_rows:
        path = fasta_dir / f"{public['case_id']}.fa"
        path.write_text(f">{public['case_id']} source=engineered_development\n"
                        f"{public['assembly_sequence']}\n", encoding="ascii")
        assembly_hashes[path.name] = sha256_file(path)
    receipt = dict(schema_version=2, development_only=True,
                   source_template_sha256=protocol["source_template_sha256"],
                   protocol_sha256=sha256_file(protocol_path), generator_sha256=sha256_file(Path(__file__)),
                   inputs_sha256=sha256_file(inputs), truth_sha256=sha256_file(truths),
                   assembly_fasta_sha256=assembly_hashes, case_count=len(pairs),
                   case_ids=[row["case_id"] for row in public_rows],
                   material_count=1, split="development", held_out_status="not_enrolled",
                   read_pairing_status="synthetic_simulated",
                   biological_accuracy_status="blocked_no_original_reads_or_native_copy_truth",
                   masking_status="separate_truth_file_not_secure_blinding")
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    build(args.protocol, args.repo_root, args.outdir)


if __name__ == "__main__":
    main()
