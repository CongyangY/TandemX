import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.scripts.controlled_collapse import generate, load_plan, score


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/controlled_collapse/toy_assembly.fa"
PLAN = ROOT / "benchmarks/controlled_collapse/toy_plan.json"


def records(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sequence(path):
    return "".join(line for line in path.read_text().splitlines() if not line.startswith(">"))


def fasta_dict(path):
    result = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            result[name] = ""
        else:
            result[name] += line
    return result


def test_generate_exact_levels_controls_and_chains(tmp_path):
    out = tmp_path / "generated"
    generate(SOURCE, PLAN, out)
    source = sequence(SOURCE)
    assert len(source) == 196
    assert sequence(out / "contraction_000.fa") == source
    assert sequence(out / "donor_swap.fa") == source
    assert len(sequence(out / "nonrepeat_deletion.fa")) == len(source) - 24
    assert len(sequence(out / "array_expansion.fa")) == len(source) + 24
    for level in (0, 25, 50, 75, 100):
        case = f"contraction_{level:03d}"
        ledger = records(out / f"{case}.ledger.tsv")
        assert len(ledger) == 2
        assert [int(x["removed_bp"]) for x in ledger] == [48 * level // 100, 48 * level // 100]
        assert len(sequence(out / f"{case}.fa")) == len(source) - sum(int(x["removed_bp"]) for x in ledger)
    gone = records(out / "contraction_100.ledger.tsv")
    assert all(x["edited_start"] == x["edited_end"] for x in gone)
    assert all(x["retained_full_units"] == "0" for x in gone)
    assert {x["operation"] for x in records(out / "contraction_100.chain.tsv")} == {"match", "delete"}
    source_records = fasta_dict(SOURCE)
    for case in ("contraction_000", "contraction_025", "contraction_100",
                 "internal_deletion_025", "boundary_left_025",
                 "nonrepeat_deletion", "array_expansion", "donor_swap"):
        edited_records = fasta_dict(out / f"{case}.fa")
        for contig in source_records:
            chain = [row for row in records(out / f"{case}.chain.tsv") if row["contig"] == contig]
            assert int(chain[0]["source_start"]) == int(chain[0]["edited_start"]) == 0
            assert int(chain[-1]["source_end"]) == len(source_records[contig])
            assert int(chain[-1]["edited_end"]) == len(edited_records[contig])
            for row in chain:
                old = source_records[contig][int(row["source_start"]):int(row["source_end"])]
                new = edited_records[contig][int(row["edited_start"]):int(row["edited_end"])]
                if row["operation"] == "match":
                    assert old == new
                elif row["operation"] == "delete":
                    assert old and not new
                else:
                    assert row["operation"] == "insert" and new and not old
    assert all(x["pairing_status"] == "invalid_donor_pair" for x in records(out / "donor_swap.ledger.tsv"))
    assert all(x["removed_bp"] == "0" for x in records(out / "nonrepeat_deletion.ledger.tsv"))
    nonrepeat_events = records(out / "nonrepeat_deletion.events.tsv")
    assert len(nonrepeat_events) == 3
    assert [int(x["removed_bp"]) for x in nonrepeat_events if x["event_id"] == "nonrepeat_control"] == [24]
    assert all(x["event_class"] == "expansion" and x["inserted_bp"] == "12"
               for x in records(out / "array_expansion.events.tsv"))
    assert all(x["event_class"] == "donor_swap" and x["pairing_status"] == "invalid_donor_pair"
               for x in records(out / "donor_swap.events.tsv"))
    internal = records(out / "internal_deletion_025.ledger.tsv")
    assert [(int(x["left_breakpoint_source"]), int(x["right_breakpoint_source"])) for x in internal] == [(38, 50), (106, 118)]
    assert all(x["edited_unit_labels"] for x in internal)
    boundary = records(out / "boundary_left_025.ledger.tsv")
    assert [(int(x["left_breakpoint_source"]), int(x["right_breakpoint_source"])) for x in boundary] == [(20, 32), (88, 100)]
    assert json.loads(records(out / "contraction_000.ledger.tsv")[1]["source_unit_labels"]) == ["A", "A", "B", "B", "A", "A", "B", "B"]
    receipt = json.loads((out / "receipt.json").read_text())
    assert len(receipt["cases"]) == 10
    assert len(receipt["source_sha256"]) == 64


def test_score_keeps_missing_and_zero_representation(tmp_path):
    out = tmp_path / "generated"
    generate(SOURCE, PLAN, out)
    predictions = tmp_path / "predictions.tsv"
    predictions.write_text("case_id\tfamily_id\tlocus_id\tstatus\tscore\tpredicted_missing_bp\n"
                           "contraction_100\ttoy_family_1\tlocus_1\tok\t0.9\t48\n"
                           "contraction_000\ttoy_family_1\tlocus_1\tabstain\t\t\n")
    scored = tmp_path / "scored"
    score(out, predictions, scored, 0.5)
    rows = records(scored / "scored.tsv")
    assert len(rows) == 20
    assert rows[0]["prediction_status"] == "abstain"
    assert rows[0]["absolute_error_bp"] == ""
    full = [row for row in rows if row["case_id"] == "contraction_100"]
    assert len(full) == 2
    assert full[0]["absolute_error_bp"] == "0"
    assert full[1]["prediction_status"] == "not_reported"
    assert json.loads((scored / "summary.json").read_text())["primary_auprc"] is None
    assert json.loads((scored / "summary.json").read_text())["injected_edit_metrics"]["tp"] is None


def test_complete_scoring_has_binary_and_continuous_metrics(tmp_path):
    out = tmp_path / "generated"
    generate(SOURCE, PLAN, out)
    predictions = tmp_path / "predictions.tsv"
    with predictions.open("w") as handle:
        handle.write("case_id\tfamily_id\tlocus_id\tstatus\tscore\tpredicted_missing_bp\n")
        for case in json.loads((out / "receipt.json").read_text())["cases"]:
            if case["case_id"] == "donor_swap":
                continue
            for row in records(out / f"{case['case_id']}.ledger.tsv"):
                missing = int(row["removed_bp"])
                handle.write(f"{row['case_id']}\t{row['family_id']}\t{row['locus_id']}\tok\t"
                             f"{0.9 if missing else 0.1}\t{missing}\n")
    scored = tmp_path / "scored"
    score(out, predictions, scored, 0.5)
    summary = json.loads((scored / "summary.json").read_text())
    metric = summary["injected_edit_metrics"]
    assert metric["status"] == "ok"
    assert (metric["denominator"], metric["tp"], metric["fn"], metric["fp"], metric["tn"]) == (18, 12, 0, 0, 6)
    assert (metric["sensitivity"], metric["fpr"], metric["missing_bp_mae"]) == (1.0, 0.0, 0.0)
    assert summary["denominator"] == 20
    assert summary["status_counts"] == {"not_reported": 2, "ok": 18}
    assert summary["primary_auprc"] is None


def test_invalid_or_overlapping_plan(tmp_path):
    plan = json.loads(PLAN.read_text())
    plan["arrays"][1]["start"] = 52
    plan["arrays"][1]["end"] = 100
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="overlapping"):
        load_plan(bad)
    plan["arrays"][1]["start"] = 88
    plan["arrays"][1]["end"] = 136
    plan["arrays"][0]["end"] = 67
    bad.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="multiple of four"):
        load_plan(bad)


def test_cli_and_reproducibility(tmp_path):
    generated = []
    for index in (1, 2):
        out = tmp_path / f"run{index}"
        completed = subprocess.run([sys.executable, "-m", "benchmarks.scripts.controlled_collapse",
                                    "generate", "--source", str(SOURCE), "--plan", str(PLAN),
                                    "--outdir", str(out)], cwd=ROOT, capture_output=True, text=True)
        assert completed.returncode == 0, completed.stderr
        generated.append(json.loads((out / "receipt.json").read_text()))
    assert generated[0] == generated[1]


def test_empty_fasta_and_reverse_complement_interval(tmp_path):
    plan = json.loads(PLAN.read_text())
    source = tmp_path / "reverse.fa"
    original = sequence(SOURCE)
    first = original[:156]
    control = original[156:]
    complement = str.maketrans("ACGT", "TGCA")
    reversed_first = first.translate(complement)[::-1]
    # Coordinates move under reverse complement; the two loci swap order.
    plan["arrays"] = [
        {"family_id": "toy_mixed_A_B", "locus_id": "rc_2", "contig": "chrToy", "start": 20, "end": 68, "unit_bp": 6},
        {"family_id": "toy_family_1", "locus_id": "rc_1", "contig": "chrToy", "start": 88, "end": 136, "unit_bp": 6},
    ]
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    source.write_text(f">chrToy\n{reversed_first}\n>chrControl\n{control}\n")
    out = tmp_path / "rc"
    generate(source, path, out)
    assert len(sequence(out / "contraction_100.fa")) == len(original) - 96
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    with pytest.raises(ValueError, match="missing edit target"):
        generate(empty, path, tmp_path / "empty_out")


def test_independent_eight_bp_unit_plan(tmp_path):
    """Exercise a second unit length without changing the shared-monomer fixture."""
    source = tmp_path / "multi_period.fa"
    source.write_text(SOURCE.read_text() + ">chrEight\n" + "TCGAGTAA" * 8 + "\n")
    plan = json.loads(PLAN.read_text())
    plan["arrays"].append({"family_id": "toy_eight_bp", "locus_id": "locus_3",
                           "contig": "chrEight", "start": 0, "end": 64, "unit_bp": 8,
                           "unit_labels": ["C"] * 8, "unit_monomers": {"C": "TCGAGTAA"}})
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    out = tmp_path / "generated"
    generate(source, path, out)
    for level, expected in ((0, 0), (25, 16), (50, 32), (75, 48), (100, 64)):
        ledger = records(out / f"contraction_{level:03d}.ledger.tsv")
        row = next(x for x in ledger if x["family_id"] == "toy_eight_bp")
        assert int(row["removed_bp"]) == expected
        assert int(row["retained_full_units"]) == (64 - expected) // 8
        assert row["unit_bp"] == "8"
    internal = next(x for x in records(out / "internal_deletion_025.ledger.tsv")
                    if x["family_id"] == "toy_eight_bp")
    assert (internal["removed_bp"], internal["left_breakpoint_source"],
            internal["right_breakpoint_source"]) == ("16", "24", "40")


def test_scoring_detects_generated_artifact_tamper(tmp_path):
    out = tmp_path / "generated"
    generate(SOURCE, PLAN, out)
    ledger = out / "contraction_100.ledger.tsv"
    ledger.write_text(ledger.read_text().replace("toy_family_1", "altered_family"))
    preds = tmp_path / "predictions.tsv"
    preds.write_text("case_id\tfamily_id\tlocus_id\tstatus\tscore\tpredicted_missing_bp\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        score(out, preds, tmp_path / "scored", 0.5)


@pytest.mark.parametrize("status,score_value,bp", [
    ("ok", "", "12"), ("ok", "NaN", "12"), ("ok", "0.5", ""),
    ("not_run", "0", ""), ("failed", "", "0"),
])
def test_prediction_state_validation(tmp_path, status, score_value, bp):
    out = tmp_path / "generated"
    generate(SOURCE, PLAN, out)
    preds = tmp_path / "predictions.tsv"
    preds.write_text("case_id\tfamily_id\tlocus_id\tstatus\tscore\tpredicted_missing_bp\n"
                     f"contraction_025\ttoy_family_1\tlocus_1\t{status}\t{score_value}\t{bp}\n")
    with pytest.raises(ValueError):
        score(out, preds, tmp_path / "scored", 0.5)
