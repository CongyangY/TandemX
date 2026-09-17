"""Score only native assembly-only endpoints supported by the frozen control."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from benchmarks.tidecluster.normalize import read_gff


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmarks/assembly_hor_comparator_v1"
TIDE = ROOT / "benchmarks/competitor_envs/tidecluster/work/assembly_hor_v1"
CEN = ROOT / "benchmarks/competitor_envs/cendetecthor/work/assembly_hor_v1"
CEN_5000 = ROOT / "benchmarks/competitor_envs/cendetecthor/work/assembly_hor_v1_window5000"
CEN_5000_KNOWN = ROOT / "benchmarks/competitor_envs/cendetecthor/work/assembly_hor_v1_window5000_known171"
EXPECTED = "f2a968383dd9bef0af38ee5891ce97b7228f67cbe71495e36340adc20d956e8e"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage(profile: str) -> dict[str, object]:
    path = BASE / "profiles" / profile / "stages.tsv"
    lines = path.read_text().splitlines()
    if len(lines) != 2:
        raise ValueError(f"Expected exactly one stage in {path}")
    keys = lines[0].split("\t")
    values = lines[1].split("\t")
    return dict(zip(keys, values, strict=True))


def union(intervals: list[tuple[int, int]]) -> int:
    total = 0
    end = 0
    for start, stop in sorted(intervals):
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def main() -> None:
    receipt = json.loads((BASE / "control/receipt.json").read_text())
    if digest(BASE / "control/control.chr1.fasta") != EXPECTED or receipt["input_sha256"] != EXPECTED:
        raise ValueError("Frozen control changed")
    truth_start, truth_end = receipt["array_interval_0based_halfopen"]
    source_files = {
        "tidehunter.gff3": TIDE / "control_tidehunter.gff3",
        "tidecluster_clustering.gff3": TIDE / "control_clustering.gff3",
        "tidecluster_TRC_1.fasta": TIDE / "control_consensus/TRC_1.fasta",
        "cendetecthor_default_window_summary.tsv": CEN / "results/wind_summary/control.chr1.txt",
        "cendetecthor_default_selected_windows.bed": CEN / "results/wind2analize/control.chr1.windows.bed",
        "cendetecthor_window5000_periodic.bed": CEN_5000 / "results/wind2analize/control.chr1.windows.bed",
        "cendetecthor_window5000_filtered.bed": CEN_5000 / "results/wind2analize/control.FULLchr.windows.filtered.bed",
        "cendetecthor_window5000_consensus.fa": CEN_5000 / "monomers/control_chr0cons.fasta",
        "cendetecthor_window5000_known171_periodic.bed": CEN_5000_KNOWN / "results/wind2analize/control.chr1.windows.bed",
        "cendetecthor_window5000_known171_candidate.bed": CEN_5000_KNOWN / "results/wind2analize/filtering/control/mon_1710bps.bed",
    }
    native = BASE / "native"
    native.mkdir(exist_ok=True)
    hashes = {}
    for name, source in source_files.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        target = native / name
        if target.exists() and digest(target) != digest(source):
            raise ValueError(f"Native archive differs: {target}")
        shutil.copyfile(source, target)
        hashes[name] = digest(target)
    rows = read_gff(native / "tidecluster_clustering.gff3")
    if any(row["sequence_id"] != f"control:0-{receipt['assembly_bp']}" for row in rows):
        raise ValueError("Unexpected TideCluster sequence ID")
    intervals = [(int(row["start"]), int(row["end"])) for row in rows]
    predicted = union(intervals)
    overlapping = union([(max(start, truth_start), min(end, truth_end))
                         for start, end in intervals if min(end, truth_end) > max(start, truth_start)])
    hunter = read_gff(native / "tidehunter.gff3")
    periods = [int(row["attributes"]["consensus_length"]) for row in hunter]
    if any(period <= 0 for period in periods):
        raise ValueError("Invalid native TideHunter period")
    cen_profile = stage("cendetecthor_corrected")
    cen_5000_profile = stage("cendetecthor_window5000_exploratory")
    cen_5000_known_profile = stage("cendetecthor_window5000_known171_exploratory")
    tidehunter_profile = stage("tidehunter")
    tidecluster_profile = stage("tidecluster")
    if any(profile["exit_code"] == "0" for profile in
           (cen_profile, cen_5000_profile, cen_5000_known_profile)) or any(
        profile["exit_code"] != "0" for profile in (tidehunter_profile, tidecluster_profile)
    ):
        raise ValueError("Run states differ from frozen scorer expectations")
    output = {
        "schema_version": 1,
        "truth_scope": receipt["truth_type"],
        "shared_input_sha256": EXPECTED,
        "TideCluster": {
            "status": "completed",
            "array_count": len(rows),
            "array_recall_bp": overlapping / (truth_end - truth_start),
            "array_precision_bp": overlapping / predicted if predicted else None,
            "predicted_union_bp": predicted,
            "truth_array_bp": truth_end - truth_start,
            "tidehunter_native_periods_bp": periods,
            "canonical_hor_bp": 855,
            "native_period_exact_canonical_hor_count": sum(x == 855 for x in periods),
            "monomer_decomposition": None,
            "hor_label_order": None,
            "wall_seconds_stage_sum": sum(float(profile["wall_seconds"]) for profile in
                                          (tidehunter_profile, tidecluster_profile)),
            "peak_process_tree_rss_mib": max(float(profile["peak_process_tree_rss_mib"])
                                              for profile in (tidehunter_profile, tidecluster_profile)),
            "rss_interpretation": "sampled process-tree peak; stage-specific max, not additive",
        },
        "CENdetectHOR": {
            "status": "technical_failure_no_selected_periodic_windows",
            "accuracy_endpoints": None,
            "wall_seconds": float(cen_profile["wall_seconds"]),
            "peak_process_tree_rss_mib": float(cen_profile["peak_process_tree_rss_mib"]),
        },
        "CENdetectHOR_exploratory": [
            {"setting": "windowSize_5000_no_monomer_prior",
             "status": "technical_failure_empty_extracted_consensus_after_1710bp_period_candidate",
             "accuracy_endpoints": None,
             "wall_seconds": float(cen_5000_profile["wall_seconds"]),
             "peak_process_tree_rss_mib": float(cen_5000_profile["peak_process_tree_rss_mib"])},
            {"setting": "windowSize_5000_expMonSize_171_supplied",
             "status": "technical_failure_no_filtered_main_bed_for_1710bp_period_candidate",
             "accuracy_endpoints": None,
             "wall_seconds": float(cen_5000_known_profile["wall_seconds"]),
             "peak_process_tree_rss_mib": float(cen_5000_known_profile["peak_process_tree_rss_mib"])},
        ],
        "TandemX": {"assembly_only_hor_endpoint": None, "status": "N/A_no_same_prior_assembly_only_run"},
        "native_sha256": hashes,
    }
    (BASE / "score_primary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, sort_keys=True))


if __name__ == "__main__":
    main()
