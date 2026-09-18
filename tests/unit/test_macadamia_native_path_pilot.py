from benchmarks.m2_routes.macadamia_native_path_pilot import technical_decision


def test_technical_decision_abstains_on_unresolved_or_mixed_reads():
    assembly = {"state": "RESOLVED", "labels": ["A+", "B+"]}
    agreed = {str(i): {"state": "RESOLVED", "labels": ["A+", "B+"]}
              for i in range(3)}
    assert technical_decision(assembly, agreed)["state"] == "SUPPORTED"
    agreed["2"] = {"state": "AMBIGUOUS", "labels": []}
    assert technical_decision(assembly, agreed)["reason"] == "one_or_more_original_read_paths_unresolved"
    agreed["2"] = {"state": "RESOLVED", "labels": ["B+", "A+"]}
    assert technical_decision(assembly, agreed)["reason"] == "mixed_original_read_paths"
    assert technical_decision(assembly, {"0": agreed["0"]})["state"] == "ABSTAIN"
    assert technical_decision({"state": "AMBIGUOUS", "reason": "no_full_monomer_tiling"},
                              agreed)["reason"] == "assembly_no_full_monomer_tiling"
