"""Small, independent controlled-array-edit benchmark (not a TandemX command).

The editor streams FASTA records and buffers only a nominated array or control
interval. Coordinates are zero-based, half-open. No reads or inference algorithm
are modified here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterator


LEVELS = (0, 25, 50, 75, 100)
PREDICTION_STATES = {"ok", "abstain", "unmapped", "failed", "not_run"}
LEDGER_FIELDS = (
    "case_id", "family_id", "locus_id", "contig", "edit_type", "level_percent",
    "source_start", "source_end", "edited_start", "edited_end", "unit_bp",
    "source_full_units", "retained_full_units", "removed_full_units",
    "source_residual_bp", "retained_residual_bp", "removed_bp",
    "inserted_bp", "source_array_sha256", "edited_array_sha256",
    "left_breakpoint_source", "right_breakpoint_source", "breakpoint_edited",
    "source_unit_labels", "edited_unit_labels", "truth_scope", "pairing_status",
)
CHAIN_FIELDS = ("case_id", "contig", "source_start", "source_end", "edited_start", "edited_end", "operation")
EVENT_FIELDS = (
    "case_id", "event_id", "event_class", "contig", "source_start", "source_end",
    "edited_start", "edited_end", "removed_bp", "inserted_bp",
    "left_breakpoint_source", "right_breakpoint_source", "breakpoint_edited",
    "source_sequence_sha256", "edited_sequence_sha256", "source_unit_labels",
    "edited_unit_labels", "truth_scope", "pairing_status",
)
SCORE_FIELDS = (
    "case_id", "family_id", "locus_id", "truth_positive", "truth_missing_bp",
    "prediction_status", "prediction_score", "predicted_missing_bp", "absolute_error_bp",
    "pairing_status", "truth_scope",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def fasta_records(path: Path) -> Iterator[tuple[str, Iterator[str]]]:
    """Yield one record at a time; sequence iterator is consumed before next record."""
    with path.open(encoding="ascii") as handle:
        pending = None
        while True:
            line = pending or handle.readline()
            pending = None
            if not line:
                return
            if not line.startswith(">"):
                raise ValueError("FASTA must start with a header")
            name = line[1:].strip().split()[0]
            if not name:
                raise ValueError("empty FASTA record name")

            def sequence_lines() -> Iterator[str]:
                nonlocal pending
                for raw in handle:
                    if raw.startswith(">"):
                        pending = raw
                        return
                    seq = raw.strip().upper()
                    if not seq or any(base not in "ACGTN" for base in seq):
                        raise ValueError(f"invalid sequence line in {name}")
                    yield seq

            yield name, sequence_lines()


def load_plan(path: Path) -> dict:
    plan = json.loads(path.read_text())
    required = {"assembly_donor_id", "reads_donor_id", "assembly_haplotype_id",
                "reads_haplotype_id", "baseline_read_consistency", "arrays", "nonrepeat_interval"}
    if not required <= plan.keys():
        raise ValueError(f"plan missing fields: {sorted(required - plan.keys())}")
    arrays = plan["arrays"]
    if not arrays:
        raise ValueError("plan needs at least one array")
    seen = set()
    intervals: dict[str, list[tuple[int, int]]] = {}
    for item in arrays:
        for key in ("family_id", "locus_id", "contig", "start", "end", "unit_bp"):
            if key not in item:
                raise ValueError(f"array missing {key}")
        key = (item["family_id"], item["locus_id"])
        if key in seen:
            raise ValueError(f"duplicate family/locus: {key}")
        seen.add(key)
        start, end, unit = item["start"], item["end"], item["unit_bp"]
        if not all(isinstance(x, int) for x in (start, end, unit)) or start < 0 or end <= start or unit <= 0:
            raise ValueError(f"invalid array coordinates: {key}")
        if (end - start) % unit or ((end - start) // unit) % 4:
            raise ValueError(f"array must contain a multiple of four complete units: {key}")
        if "unit_labels" in item:
            labels = item["unit_labels"]
            monomers = item.get("unit_monomers", {})
            if (not isinstance(labels, list) or len(labels) != (end - start) // unit
                    or not all(isinstance(x, str) and x in monomers and len(monomers[x]) == unit for x in labels)):
                raise ValueError(f"invalid unit labels/monomers: {key}")
        intervals.setdefault(item["contig"], []).append((start, end))
    control = plan["nonrepeat_interval"]
    if not {"contig", "start", "end"} <= control.keys():
        raise ValueError("nonrepeat_interval requires contig/start/end")
    if not isinstance(control["start"], int) or not isinstance(control["end"], int) or control["start"] < 0 or control["end"] <= control["start"]:
        raise ValueError("invalid nonrepeat interval")
    intervals.setdefault(control["contig"], []).append((control["start"], control["end"]))
    for contig, spans in intervals.items():
        ordered = sorted(spans)
        if any(a[1] > b[0] for a, b in zip(ordered, ordered[1:])):
            raise ValueError(f"overlapping edited intervals on {contig}")
    return plan


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def case_specs(plan: dict) -> list[tuple[str, int, str]]:
    return [(f"contraction_{level:03d}", level, "contraction") for level in LEVELS] + [
        ("internal_deletion_025", 25, "internal_deletion"),
        ("boundary_left_025", 25, "boundary_left"),
        ("nonrepeat_deletion", 0, "nonrepeat_deletion"),
        ("array_expansion", 25, "expansion"),
        ("donor_swap", 0, "donor_swap"),
    ]


def render_case(source: Path, output: Path, plan: dict, case_id: str, level: int, kind: str) -> tuple[list[dict], list[dict], list[dict]]:
    """Render one case with constant memory outside nominated edit intervals."""
    edits = {item["contig"]: [] for item in plan["arrays"]}
    for item in plan["arrays"]:
        edits[item["contig"]].append((item["start"], item["end"], item))
    if kind == "nonrepeat_deletion":
        control = plan["nonrepeat_interval"]
        edits.setdefault(control["contig"], []).append((control["start"], control["end"], control))
    for spans in edits.values():
        spans.sort(key=lambda x: x[0])
    ledgers, chain, events = [], [], []
    seen_contigs = set()
    seen_edits = set()
    with output.open("w", encoding="ascii") as dest:
        for contig, lines in fasta_records(source):
            if contig in seen_contigs:
                raise ValueError(f"duplicate FASTA record: {contig}")
            seen_contigs.add(contig)
            dest.write(f">{contig}\n")
            spans = edits.get(contig, [])
            source_pos = edited_pos = 0
            chunk_buffer = ""
            span_index = 0

            def emit_segment(start: int, end: int, value: str, operation: str) -> None:
                nonlocal edited_pos
                if value:
                    dest.write(value + "\n")
                if (operation == "match" and chain and chain[-1]["case_id"] == case_id
                        and chain[-1]["contig"] == contig and chain[-1]["operation"] == "match"
                        and chain[-1]["source_end"] == start and chain[-1]["edited_end"] == edited_pos):
                    chain[-1]["source_end"] = end
                    chain[-1]["edited_end"] += len(value)
                else:
                    chain.append(dict(case_id=case_id, contig=contig, source_start=start,
                                      source_end=end, edited_start=edited_pos,
                                      edited_end=edited_pos + len(value), operation=operation))
                edited_pos += len(value)

            def consume_until(target: int) -> str:
                nonlocal chunk_buffer, source_pos
                chunks = []
                while source_pos < target:
                    if not chunk_buffer:
                        try:
                            chunk_buffer = next(line_iter)
                        except StopIteration as exc:
                            raise ValueError(f"interval exceeds contig {contig}") from exc
                    take = min(target - source_pos, len(chunk_buffer))
                    chunks.append(chunk_buffer[:take])
                    chunk_buffer = chunk_buffer[take:]
                    source_pos += take
                return "".join(chunks)

            line_iter = iter(lines)
            while span_index < len(spans):
                start, end, item = spans[span_index]
                if start > source_pos:
                    # The inter-interval sequence is emitted in bounded chunks.
                    while source_pos < start:
                        upto = min(start, source_pos + 65536)
                        old = source_pos
                        emit_segment(old, upto, consume_until(upto), "match")
                original = consume_until(end)
                source_labels = item.get("unit_labels", [])
                if source_labels:
                    unit = item["unit_bp"]
                    expected_sequence = "".join(item["unit_monomers"][label] for label in source_labels)
                    if original != expected_sequence:
                        raise ValueError(f"unit-label sequence mismatch: {contig}:{start}-{end}")
                if item is plan["nonrepeat_interval"]:
                    edited = ""
                    removed_units = 0
                    removed_bp = len(original)
                    edit_type = "nonrepeat_deletion"
                    deleted_start = start
                    retained_prefix = ""
                    retained_suffix = ""
                    edited_labels = []
                else:
                    unit = item["unit_bp"]
                    units = len(original) // unit
                    if kind in ("contraction", "internal_deletion", "boundary_left"):
                        removed_units = units * level // 100
                        first_removed = (units - removed_units if kind == "contraction" else
                                         (units - removed_units) // 2 if kind == "internal_deletion" else 0)
                        deleted_start = start + first_removed * unit
                        retained_prefix = original[:first_removed * unit]
                        retained_suffix = original[(first_removed + removed_units) * unit:]
                        edited = retained_prefix + retained_suffix
                        edited_labels = source_labels[:first_removed] + source_labels[first_removed + removed_units:] if source_labels else []
                        edit_type = kind if level else "unedited"
                    elif kind == "expansion":
                        removed_units = 0
                        edited = original + original[:units * level // 100 * unit]
                        deleted_start = end
                        retained_prefix = original
                        retained_suffix = ""
                        edited_labels = source_labels + source_labels[:units * level // 100] if source_labels else []
                        edit_type = "expansion"
                    else:
                        removed_units = 0
                        edited = original
                        deleted_start = end
                        retained_prefix = original
                        retained_suffix = ""
                        edited_labels = source_labels.copy()
                        edit_type = "donor_swap" if kind == "donor_swap" else "unedited"
                    removed_bp = removed_units * unit
                new_start = edited_pos
                breakpoint = new_start + len(retained_prefix)
                if edit_type == "expansion":
                    emit_segment(start, end, original, "match")
                    emit_segment(end, end, edited[len(original):], "insert")
                else:
                    if retained_prefix:
                        emit_segment(start, deleted_start, retained_prefix, "match")
                    if removed_bp:
                        emit_segment(deleted_start, deleted_start + removed_bp, "", "delete")
                    if retained_suffix:
                        emit_segment(deleted_start + removed_bp, end, retained_suffix, "match")
                    if not retained_prefix and not retained_suffix and not removed_bp:
                        emit_segment(start, end, edited, "match")
                pairing = "invalid_donor_pair" if kind == "donor_swap" else (
                    "declared_same_donor_baseline_consistent" if plan["assembly_donor_id"] == plan["reads_donor_id"]
                    and plan["assembly_haplotype_id"] == plan["reads_haplotype_id"]
                    and plan["baseline_read_consistency"] is True else "unverified_pair")
                events.append(dict(case_id=case_id,
                                   event_id=("nonrepeat_control" if item is plan["nonrepeat_interval"]
                                             else f"{item['family_id']}/{item['locus_id']}"),
                                   event_class=edit_type, contig=contig, source_start=start, source_end=end,
                                   edited_start=new_start, edited_end=edited_pos,
                                   removed_bp=removed_bp, inserted_bp=max(0, len(edited) - len(original)),
                                   left_breakpoint_source=deleted_start,
                                   right_breakpoint_source=deleted_start + removed_bp,
                                   breakpoint_edited=breakpoint,
                                   source_sequence_sha256=sha256_text(original),
                                   edited_sequence_sha256=sha256_text(edited),
                                   source_unit_labels=json.dumps(source_labels, separators=(",", ":")),
                                   edited_unit_labels=json.dumps(edited_labels, separators=(",", ":")),
                                   truth_scope="injected_edit_delta_only", pairing_status=pairing))
                if item is not plan["nonrepeat_interval"]:
                    ledgers.append(dict(case_id=case_id, family_id=item["family_id"], locus_id=item["locus_id"],
                                        contig=contig, edit_type=edit_type, level_percent=level,
                                        source_start=start, source_end=end, edited_start=new_start,
                                        edited_end=edited_pos, unit_bp=item["unit_bp"],
                                        source_full_units=len(original) // item["unit_bp"],
                                        retained_full_units=len(edited) // item["unit_bp"],
                                        removed_full_units=removed_units, source_residual_bp=0,
                                        retained_residual_bp=0, removed_bp=removed_bp,
                                        inserted_bp=max(0, len(edited) - len(original)),
                                        source_array_sha256=sha256_text(original), edited_array_sha256=sha256_text(edited),
                                        left_breakpoint_source=deleted_start,
                                        right_breakpoint_source=deleted_start + removed_bp,
                                        breakpoint_edited=breakpoint,
                                        source_unit_labels=json.dumps(source_labels, separators=(",", ":")),
                                        edited_unit_labels=json.dumps(edited_labels, separators=(",", ":")),
                                        truth_scope="injected_edit_delta_only", pairing_status=pairing))
                seen_edits.add((contig, start, end))
                span_index += 1
            # Copy tail without holding entire contig.
            if chunk_buffer:
                old = source_pos
                emit_segment(old, old + len(chunk_buffer), chunk_buffer, "match")
                source_pos += len(chunk_buffer)
            for chunk in line_iter:
                old = source_pos
                emit_segment(old, old + len(chunk), chunk, "match")
                source_pos += len(chunk)
    expected = {(contig, a, b) for contig, spans in edits.items() for a, b, _ in spans}
    if seen_edits != expected:
        raise ValueError(f"missing edit target(s): {sorted(expected - seen_edits)}")
    return ledgers, chain, events


def generate(source: Path, plan_path: Path, outdir: Path) -> None:
    if not source.is_file() or not plan_path.is_file():
        raise ValueError("source FASTA and plan JSON must exist")
    plan = load_plan(plan_path)
    outdir.mkdir(parents=True, exist_ok=False)
    source_hash = sha256_file(source)
    plan_hash = sha256_file(plan_path)
    manifests = []
    for case_id, level, kind in case_specs(plan):
        fasta = outdir / f"{case_id}.fa"
        ledger, chain, events = render_case(source, fasta, plan, case_id, level, kind)
        if len(ledger) != len(plan["arrays"]):
            raise ValueError(f"incomplete array ledger: {case_id}")
        write_tsv(outdir / f"{case_id}.ledger.tsv", LEDGER_FIELDS, ledger)
        write_tsv(outdir / f"{case_id}.chain.tsv", CHAIN_FIELDS, chain)
        write_tsv(outdir / f"{case_id}.events.tsv", EVENT_FIELDS, events)
        manifests.append(dict(case_id=case_id, edit_type=kind, level_percent=level,
                              fasta_sha256=sha256_file(fasta),
                              ledger_sha256=sha256_file(outdir / f"{case_id}.ledger.tsv"),
                              chain_sha256=sha256_file(outdir / f"{case_id}.chain.tsv"),
                              events_sha256=sha256_file(outdir / f"{case_id}.events.tsv")))
    receipt = dict(source_fasta=str(source.resolve()), source_sha256=source_hash,
                   plan_json=str(plan_path.resolve()), plan_sha256=plan_hash,
                   generator_sha256=sha256_file(Path(__file__)), cases=manifests,
                   truth_scope="injected_edit_delta_only",
                   note="Input assembly is a reference proxy; no physical copy-number truth is asserted.")
    (outdir / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def score(outdir: Path, predictions: Path, scoredir: Path, decision_threshold: float) -> None:
    """Score a fixed denominator. Missing calls remain missing, never zero."""
    if not math.isfinite(decision_threshold) or not 0 <= decision_threshold <= 1:
        raise ValueError("decision threshold must be finite and in [0,1]")
    receipt = json.loads((outdir / "receipt.json").read_text())
    with predictions.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        need = {"case_id", "family_id", "locus_id", "status", "score", "predicted_missing_bp"}
        if reader.fieldnames is None or not need <= set(reader.fieldnames):
            raise ValueError("prediction TSV missing required columns")
        predictions_by_key = {}
        for row in reader:
            key = (row["case_id"], row["family_id"], row["locus_id"])
            if key in predictions_by_key or row["status"] not in PREDICTION_STATES:
                raise ValueError(f"duplicate or invalid prediction: {key}")
            if row["status"] == "ok":
                value = float(row["score"])
                bp = int(row["predicted_missing_bp"])
                if not math.isfinite(value) or not 0 <= value <= 1 or bp < 0:
                    raise ValueError(f"invalid score or bp: {key}")
            elif row["score"] or row["predicted_missing_bp"]:
                raise ValueError(f"non-ok prediction must leave score and bp blank: {key}")
            predictions_by_key[key] = row
    rows = []
    for case in receipt["cases"]:
        for suffix, expected in (("fa", "fasta_sha256"), ("ledger.tsv", "ledger_sha256"),
                                 ("chain.tsv", "chain_sha256"), ("events.tsv", "events_sha256")):
            artifact = outdir / f"{case['case_id']}.{suffix}"
            if sha256_file(artifact) != case[expected]:
                raise ValueError(f"generated artifact hash mismatch: {artifact}")
        with (outdir / f"{case['case_id']}.ledger.tsv").open(newline="") as handle:
            for truth in csv.DictReader(handle, delimiter="\t"):
                key = (truth["case_id"], truth["family_id"], truth["locus_id"])
                pred = predictions_by_key.pop(key, None)
                state = pred["status"] if pred else "not_reported"
                bp = int(truth["removed_bp"])
                ok = state == "ok"
                rows.append(dict(case_id=key[0], family_id=key[1], locus_id=key[2],
                                 truth_positive=int(bp > 0), truth_missing_bp=bp,
                                 prediction_status=state, prediction_score=pred["score"] if ok else "",
                                 predicted_missing_bp=pred["predicted_missing_bp"] if ok else "",
                                 absolute_error_bp=(abs(int(pred["predicted_missing_bp"]) - bp)
                                                    if ok and truth["pairing_status"] != "invalid_donor_pair" else ""),
                                 pairing_status=truth["pairing_status"], truth_scope=truth["truth_scope"]))
    if predictions_by_key:
        raise ValueError(f"unknown prediction keys: {sorted(predictions_by_key)}")
    scoredir.mkdir(parents=True, exist_ok=False)
    write_tsv(scoredir / "scored.tsv", SCORE_FIELDS, rows)
    technical = [row for row in rows if row["pairing_status"] != "invalid_donor_pair"]
    complete = bool(technical) and all(row["prediction_status"] == "ok" for row in technical)
    positives = sum(row["truth_positive"] for row in technical)
    negatives = len(technical) - positives
    metric = dict(status="ok" if complete and positives and negatives else "blocked",
                  reason="" if complete and positives and negatives else
                  "incomplete_predictions_or_single_truth_class",
                  threshold=decision_threshold, technical_injected_edit_denominator=len(technical),
                  scope="injected_edit_technical_only",
                  tp=None, fn=None, fp=None, tn=None, sensitivity=None, fpr=None,
                  missing_bp_mae=None)
    if metric["status"] == "ok":
        counts = dict(tp=0, fn=0, fp=0, tn=0)
        for row in technical:
            true_positive = bool(row["truth_positive"])
            predicted_positive = float(row["prediction_score"]) >= decision_threshold
            label = ("tp" if true_positive else "fp") if predicted_positive else (
                "fn" if true_positive else "tn")
            counts[label] += 1
        metric.update(counts)
        metric["sensitivity"] = counts["tp"] / positives
        metric["fpr"] = counts["fp"] / negatives
        metric["missing_bp_mae"] = sum(row["absolute_error_bp"] for row in technical) / len(technical)
    summary = dict(denominator=len(rows), status_counts={state: sum(x["prediction_status"] == state for x in rows)
                                                         for state in sorted({x["prediction_status"] for x in rows})},
                   declared_consistent_pair_rows=sum(x["pairing_status"] == "declared_same_donor_baseline_consistent" for x in rows),
                   positive_rows=sum(x["truth_positive"] for x in rows),
                   zero_representation_rows=sum(x["case_id"] == "contraction_100" for x in rows),
                   injected_edit_metrics=metric,
                   verified_same_donor_read_baseline_rows=0,
                   read_assembly_accuracy_status="blocked_no_independent_pairing_and_baseline_evidence",
                   primary_auprc=None,
                   primary_auprc_reason="Requires preregistered eligible subset, baseline pairing, and full scored predictions; not computed by generator.",
                   truth_scope="injected_edit_delta_only", prediction_sha256=sha256_file(predictions),
                   receipt_sha256=sha256_file(outdir / "receipt.json"))
    (scoredir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    gen = sub.add_parser("generate", help="Generate controlled edits and exact edit ledgers")
    gen.add_argument("--source", type=Path, required=True)
    gen.add_argument("--plan", type=Path, required=True)
    gen.add_argument("--outdir", type=Path, required=True)
    scoring = sub.add_parser("score", help="Join predictions to the fixed ledger denominator")
    scoring.add_argument("--generated", type=Path, required=True)
    scoring.add_argument("--predictions", type=Path, required=True)
    scoring.add_argument("--outdir", type=Path, required=True)
    scoring.add_argument("--decision-threshold", type=float, required=True,
                         help="Frozen inclusive score threshold for binary injected-edit detection")
    args = parser.parse_args()
    if args.command == "generate":
        generate(args.source, args.plan, args.outdir)
    else:
        score(args.generated, args.predictions, args.outdir, args.decision_threshold)


if __name__ == "__main__":
    main()
