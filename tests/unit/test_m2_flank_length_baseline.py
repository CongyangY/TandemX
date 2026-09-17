from benchmarks.m2_routes.baselines.flank_length import predict_case


def _row(assembly_array: str, reads: list[str]) -> dict:
    left = "ACGTTGCATGCATGCACTGATCGTACGTGCTA"
    right = "TTGCAACGATCGACCTAGGCTAACGTTAGGCA"
    return {
        "case_id": "case",
        "left_flank_sequence": left,
        "right_flank_sequence": right,
        "assembly_sequence": left + assembly_array + right,
        "raw_read_sequences": {
            str(index): left + read + right for index, read in enumerate(reads)
        },
    }


def test_length_baseline_detects_large_span_loss() -> None:
    row = _row("ACGT" * 15, ["ACGT" * 20] * 3)
    result = predict_case(row)
    assert result["audit_state"] == "DISCORDANT"
    assert result["predicted_signed_bp_delta"] == 20


def test_length_baseline_does_not_detect_order_only_edit() -> None:
    row = _row("ACGT" * 15, ["TGCA" * 15] * 3)
    assert predict_case(row)["audit_state"] == "SUPPORTED"


def test_length_baseline_rejects_partial_reads() -> None:
    row = _row("ACGT" * 15, ["ACGT" * 15] * 3)
    row["raw_read_sequences"] = {str(i): row["left_flank_sequence"] + "ACGT" for i in range(3)}
    assert predict_case(row)["audit_state"] == "INSUFFICIENT_READ_SUPPORT"
