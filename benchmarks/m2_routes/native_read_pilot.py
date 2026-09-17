"""Research-only C3 native-read feasibility pilot for frozen M2 alignment.

The assembly edits are controlled development inputs. Original HiFi molecules
are never altered. Donor/haplotype pairing and genome-wide uniqueness remain
unverified, so results below are technical paths, not biological calls.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

import edlib

from benchmarks.m2_routes.alignment.prototype import decompose


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/controlled_collapse/native_pair_audit_20260917"
EDIT = ROOT / "benchmarks/controlled_collapse/native_edit_development_v1"
B1 = ROOT / "benchmarks/controlled_collapse/v3/development_bundle/inputs.jsonl"
FLANK_BP = 1024  # Development pilot choice after observing shorter repeated flanks.
MAX_FLANK_ERROR_FRACTION = 0.05


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError("duplicate FASTA ID")
            records[name] = ""
        elif name is None:
            raise ValueError("FASTA sequence before header")
        else:
            records[name] += line.strip().upper()
    if not records:
        raise ValueError("empty FASTA")
    return records


def fastq_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    with gzip.open(path, "rt") as handle:
        while header := handle.readline():
            sequence, plus, quality = (handle.readline() for _ in range(3))
            if not header.startswith("@") or not plus.startswith("+") or not quality:
                raise ValueError("invalid FASTQ record")
            sequence, quality = sequence.strip().upper(), quality.strip()
            if len(sequence) != len(quality):
                raise ValueError("FASTQ sequence and quality lengths differ")
            read_id = header[1:].split()[0]
            if read_id in records:
                raise ValueError("duplicate FASTQ ID")
            records[read_id] = sequence
    if not records:
        raise ValueError("empty FASTQ")
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def unique_flank_hit(flank: str, read: str) -> tuple[int, int, int] | None:
    if not flank or not read:
        return None
    limit = max(2, math.ceil(len(flank) * MAX_FLANK_ERROR_FRACTION))
    result = edlib.align(flank, read, mode="HW", task="locations", k=limit)
    if result["editDistance"] < 0 or len(result["locations"]) != 1:
        return None
    start, end_inclusive = result["locations"][0]
    return start, end_inclusive + 1, int(result["editDistance"])


def trim_between_flanks(read: str, left: str, right: str) -> dict[str, object]:
    a, b = unique_flank_hit(left, read), unique_flank_hit(right, read)
    if a is None or b is None:
        return {"status": "ABSTAIN", "reason": "missing_or_nonunique_natural_flank"}
    if a[1] > b[0]:
        return {"status": "ABSTAIN", "reason": "flank_order_or_overlap_invalid"}
    start, end = a[1], b[0]
    return {"status": "ELIGIBLE", "reason": "two_unique_natural_flanks_in_read",
            "trim_interval_0based": [start, end], "array_sequence": read[start:end],
            "left_hit_0based": [a[0], a[1]], "right_hit_0based": [b[0], b[1]],
            "left_edit_cost": a[2], "right_edit_cost": b[2]}


def primary_paf_array_interval(paf_path: Path, read_id: str, context_id: str,
                               array_start: int, array_end: int) -> list[int]:
    """Project reference boundaries onto an oriented query using the primary CIGAR."""
    matches: list[list[int]] = []
    for line in paf_path.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] != read_id or not fields[5].startswith(context_id + "|") \
                or "tp:A:P" not in fields[12:]:
            continue
        cigars = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
        if len(cigars) != 1:
            raise ValueError("primary PAF CIGAR absent or repeated")
        query_length, query_start, query_end = map(int, fields[1:4])
        strand = fields[4]
        if strand not in {"+", "-"}:
            raise ValueError("invalid PAF strand")
        query = query_start if strand == "+" else query_length - query_end
        reference = int(fields[7])
        boundary_map: dict[int, int] = {}
        for count_text, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigars[0]):
            count = int(count_text)
            if operation in "M=X":
                for boundary in (array_start, array_end):
                    if reference <= boundary <= reference + count:
                        candidate = query + boundary - reference
                        if boundary in boundary_map and boundary_map[boundary] != candidate:
                            raise ValueError("ambiguous reference boundary in CIGAR")
                        boundary_map[boundary] = candidate
                query += count
                reference += count
            elif operation in "DN":
                reference += count
            elif operation == "I":
                query += count
            else:
                raise ValueError("unsupported clipping in PAF CIGAR")
        if len(boundary_map) != 2:
            raise ValueError("array boundary not in aligned match block")
        matches.append([boundary_map[array_start], boundary_map[array_end]])
    if len(matches) != 1:
        raise ValueError("expected exactly one C3 primary PAF alignment")
    return matches[0]


def path_summary(sequence: str, monomers: dict[str, str]) -> dict[str, object]:
    result = decompose(sequence, monomers)
    return {"state": result.state, "reason": result.reason,
            "score": result.score, "alternative_score": result.alternative_score,
            "labels": [copy.label + copy.orientation for copy in result.copies],
            "copy_count": len(result.copies)}


def run(output: Path) -> dict[str, object]:
    if output.exists():
        raise FileExistsError("pilot output exists; choose a new path to preserve prior evidence")
    source_receipt = json.loads((SOURCE / "source_and_alignment_receipt.json").read_text())
    source_config = json.loads((EDIT / "source_config.json").read_text())
    edit_receipt = json.loads((EDIT / "generated/receipt.json").read_text())
    context_path = SOURCE / "colcen_native_contexts.fa"
    read_path = SOURCE / "colcen_native_spanners.fastq.gz"
    if sha256(context_path) != source_receipt["context_fasta_sha256"] or \
            sha256(read_path) != source_receipt["native_spanners_fastq_sha256"]:
        raise ValueError("frozen native source hash mismatch")
    paf_path = SOURCE / "colcen_native_spanners.paf"
    if sha256(paf_path) != source_receipt["native_spanners_paf_sha256"]:
        raise ValueError("frozen native PAF hash mismatch")
    if source_config["read_pairing_status"] != "lineage_context_supported_donor_unverified" \
            or edit_receipt["read_pairing_status"] != source_config["read_pairing_status"]:
        raise ValueError("unsupported donor pairing status")
    contexts = fasta_records(context_path)
    context = next(c for c in source_config["contexts"] if c["array_id"] == "C3")
    original = contexts[context["reference_record_id"]]
    start, end = context["array_start0"], context["array_end0"]
    if hashlib.sha256(original[start:end].encode()).hexdigest() != context["array_sequence_sha256"]:
        raise ValueError("C3 source array hash mismatch")
    left, right = original[start - FLANK_BP:start], original[end:end + FLANK_BP]
    reads = fastq_records(read_path)
    alignments = [x for x in source_receipt["native_molecules"] if x["context"] == "C3"]
    expected_ids = set(context["supporting_read_ids"])
    if len(alignments) != 3 or {x["read_id"] for x in alignments} != expected_ids:
        raise ValueError("C3 read IDs or alignment count changed")
    trims: list[dict[str, object]] = []
    for alignment in sorted(alignments, key=lambda x: x["read_id"]):
        read_id = alignment["read_id"]
        read = reads[read_id]
        if len(read) != alignment["read_length_bp"] or alignment["identity"] < 0.99 or \
                alignment["local_context_mapq"] < 20 or \
                min(alignment["natural_left_flank_bp"], alignment["natural_right_flank_bp"]) < 1000:
            raise ValueError("frozen alignment eligibility changed")
        oriented = read if alignment["strand"] == "+" else reverse_complement(read)
        record = trim_between_flanks(oriented, left, right)
        projected = primary_paf_array_interval(paf_path, read_id, "C3", start, end)
        if record["status"] == "ELIGIBLE":
            observed = record["trim_interval_0based"]
            if any(abs(a-b) > 10 for a, b in zip(observed, projected)):
                record = {"status": "ABSTAIN", "reason": "flank_PAF_boundary_disagreement",
                          "observed_trim_interval_0based": observed}
        record["primary_PAF_projected_interval_0based"] = projected
        record.update(read_id=read_id, original_read_sha256=hashlib.sha256(read.encode()).hexdigest(),
                      strand=alignment["strand"], oriented_read_length_bp=len(oriented),
                      original_read_length_bp=len(read), local_context_mapq=alignment["local_context_mapq"],
                      local_context_identity=alignment["identity"],
                      other_context_max_identity=alignment["other_context_max_identity"])
        trims.append(record)
    eligible = {r["read_id"]: r for r in trims if r["status"] == "ELIGIBLE"}
    with B1.open() as handle:
        first = json.loads(handle.readline())
    monomers = first["candidate_monomers"]
    if set(monomers) != {"C1", "C2", "C3"}:
        raise ValueError("frozen B1 monomer catalogue changed")
    cases = [("C3_source_original", original, start, end, 0, "source_original")]
    for item in edit_receipt["cases"]:
        if item["array_id"] != "C3":
            continue
        edited_path = EDIT / "generated" / (item["case_id"] + ".fa")
        if sha256(edited_path) != item["edited_fasta_sha256"]:
            raise ValueError("edited C3 FASTA hash mismatch")
        edited = next(iter(fasta_records(edited_path).values()))
        s, e = item["edited_array_interval"]
        if edited[:s] != original[:start] or edited[e:] != original[end:]:
            raise ValueError("edited C3 natural flanks changed")
        cases.append((item["case_id"], edited, s, e, item["injected_deleted_bp"], item["edit_type"]))
    if len(cases) != 10:
        raise ValueError("expected original plus nine C3 edits")
    read_paths = {rid: path_summary(str(r["array_sequence"]), monomers) for rid, r in eligible.items()}
    rows = []
    for case_id, sequence, s, e, deleted, edit_type in cases:
        assembly_trim = trim_between_flanks(sequence, left, right)
        if assembly_trim["status"] != "ELIGIBLE" or assembly_trim["trim_interval_0based"] != [s, e]:
            raise ValueError("edited assembly flank coordinates ambiguous")
        assembly_path = path_summary(str(assembly_trim["array_sequence"]), monomers)
        read_lengths = [len(str(r["array_sequence"])) for r in eligible.values()]
        baseline = {"state": "INSUFFICIENT_READ_SUPPORT", "reason": "fewer_than_three_eligible_reads",
                    "median_read_span_bp": None, "signed_bp_delta_assembly_minus_read": None}
        if len(read_lengths) >= 3:
            median_length = statistics.median(read_lengths)
            delta = (e - s) - median_length
            threshold = max(2, math.ceil(max(e - s, median_length) * 0.05))
            baseline = {"state": "DISCORDANT" if abs(delta) > threshold else "SUPPORTED",
                        "reason": "fixed_five_percent_span_threshold", "median_read_span_bp": median_length,
                        "signed_bp_delta_assembly_minus_read": delta, "threshold_bp": threshold}
        technical = {"state": "ABSTAIN", "reason": "source_pairing_unverified_for_biological_audit"}
        if assembly_path["state"] != "RESOLVED":
            technical["reason"] = "assembly_" + str(assembly_path["reason"])
        elif len(eligible) < 3:
            technical["reason"] = "fewer_than_three_eligible_reads"
        elif any(p["state"] != "RESOLVED" for p in read_paths.values()):
            technical["reason"] = "one_or_more_native_read_paths_unresolved"
        else:
            paths = [tuple(p["labels"]) for p in read_paths.values()]
            if len(set(paths)) != 1:
                technical["reason"] = "mixed_native_read_paths"
            else:
                technical = {"state": "DISCORDANT" if paths[0] != tuple(assembly_path["labels"])
                             else "SUPPORTED", "reason": "technical_label_path_comparison_only"}
        rows.append({"case_id": case_id, "edit_type": edit_type, "injected_deleted_bp": deleted,
                     "assembly_array_interval_0based": [s, e], "assembly_array_bp": e-s,
                     "assembly_path": assembly_path, "native_read_paths": read_paths,
                     "m2_technical": technical, "length_baseline": baseline,
                     "biological_audit_state": "NOT_EVALUATED_PAIRING_UNVERIFIED"})
    output.mkdir(parents=True, exist_ok=False)
    for record in trims:
        record.pop("array_sequence", None)
    (output / "read_trims.json").write_text(json.dumps(trims, indent=2, sort_keys=True) + "\n")
    with (output / "per_case.jsonl").open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    receipt = {"case_count": len(rows), "eligible_read_count": len(eligible),
               "read_ids": sorted(expected_ids), "flank_bp": FLANK_BP,
               "flank_max_edit_fraction": MAX_FLANK_ERROR_FRACTION,
               "source_context_sha256": sha256(context_path), "native_fastq_sha256": sha256(read_path),
               "native_PAF_sha256": sha256(paf_path),
               "source_receipt_sha256": sha256(SOURCE / "source_and_alignment_receipt.json"),
               "source_config_sha256": sha256(EDIT / "source_config.json"),
               "edit_receipt_sha256": sha256(EDIT / "generated/receipt.json"),
               "frozen_m2_prototype_sha256": sha256(ROOT / "benchmarks/m2_routes/alignment/prototype.py"),
               "frozen_b1_inputs_sha256": sha256(B1),
               "read_trims_sha256": sha256(output / "read_trims.json"),
               "per_case_sha256": sha256(output / "per_case.jsonl"),
               "pairing_status": source_config["read_pairing_status"],
               "biological_accuracy_status": "not_evaluated",
               "truth_scope": "injected_edit_delta_only",
               "pilot_scope": "single_C3_development_lineage"}
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
