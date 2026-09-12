#!/usr/bin/env python3
"""Score frozen historical-assembly family prioritization without changing callers."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


THRESHOLD = 0.6


def spearman(xs: list[float], ys: list[float]) -> float | str:
    if len(xs) < 2: return "N/A"
    def ranks(values: list[float]) -> list[float]:
        order=sorted(range(len(values)), key=values.__getitem__); out=[0.0]*len(values); i=0
        while i < len(order):
            j=i
            while j+1 < len(order) and values[order[j+1]] == values[order[i]]: j += 1
            rank=(i+j+2)/2
            for k in range(i,j+1): out[order[k]]=rank
            i=j+1
        return out
    a,b=ranks(xs),ranks(ys); ma=sum(a)/len(a); mb=sum(b)/len(b)
    den=math.sqrt(sum((v-ma)**2 for v in a)*sum((v-mb)**2 for v in b))
    return sum((u-ma)*(v-mb) for u,v in zip(a,b))/den if den else "N/A"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evaluate(reference: list[dict[str, str]], estimates: dict[str, dict[str, str]] | None, method: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    scored = []
    for r in reference:
        if r["eligibility"] != "eligible":
            continue
        if r.get("reference_state") not in {"reference_collapse", "reference_retained"}:
            scored.append({"family_id": r["family_id"], "reference_proxy_positive": "N/A", "method": method,
                           "state": "unresolved_reference_proxy_state", "predicted_positive": "N/A", "outcome": "N/A"})
            continue
        truth = r["reference_state"] == "reference_collapse"
        e = (estimates or {}).get(r["family_id"])
        if e is None or e.get("state") != "ok":
            state = "N/A" if e is None else str(e.get("state"))
            scored.append({"family_id": r["family_id"], "reference_proxy_positive": truth, "method": method, "state": state, "predicted_positive": "N/A", "outcome": "N/A"})
            continue
        bp = float(e["read_estimated_bp"])
        if bp <= 0:
            scored.append({"family_id": r["family_id"], "reference_proxy_positive": truth, "method": method,
                           "state": "nonpositive_read_estimate", "predicted_positive": "N/A", "outcome": "N/A"})
            continue
        old = float(r["old_assembly_bp"])
        call = bp > 0 and old / bp < THRESHOLD
        outcome = "TP" if truth and call else "FN" if truth else "FP" if call else "TN"
        scored.append({"family_id": r["family_id"], "reference_proxy_positive": truth, "method": method, "state": "ok", "read_estimated_bp": bp, "old_assembly_bp": old, "old_read_ratio": old / bp, "predicted_positive": call, "outcome": outcome})
    c = Counter(x["outcome"] for x in scored)
    available = sum(x["outcome"] != "N/A" for x in scored)
    counts: object = {"TP": c["TP"], "FN": c["FN"], "FP": c["FP"], "TN": c["TN"]} if available else "N/A"
    paired=[(max(float(x["read_estimated_bp"])-float(x["old_assembly_bp"]),0), float(next(r["observed_gain_bp"] for r in reference if r["family_id"] == x["family_id"]))) for x in scored if x["state"] == "ok" and x["outcome"] != "N/A"]
    rho=spearman([x[0] for x in paired],[x[1] for x in paired]) if paired else "N/A"
    return scored, {"method": method, "eligible_families": len(scored), "available_families": available, "unavailable_families": c["N/A"], "reference_proxy_confusion_matrix": counts, "reference_proxy_recall": c["TP"]/(c["TP"]+c["FN"]) if available and c["TP"]+c["FN"] else "N/A", "reference_proxy_precision": c["TP"]/(c["TP"]+c["FP"]) if available and c["TP"]+c["FP"] else "N/A", "ranking_concordance_spearman": rho, "ranking_concordance_state": "assessed_read_deficit_vs_new_old_gain" if rho != "N/A" else "N/A_no_available_estimates", "interpretation_boundary": "retrospective agreement with a newer-assembly reference proxy; not independent physical truth or a population accuracy estimate"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--ey15", type=Path, required=True); p.add_argument("--macadamia", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True); args = p.parse_args()
    if args.outdir.exists(): raise FileExistsError(args.outdir)
    args.outdir.mkdir(parents=True)
    refs = {"ey15": rows(args.ey15), "macadamia": rows(args.macadamia)}
    all_rows=[]; summaries=[]
    for species, reference in refs.items():
        tx = {r["family_id"]: {"state":"ok", "read_estimated_bp":r["read_estimated_bp"]} for r in reference}
        for method, estimate in (("tandemx_frozen", tx), ("srf_k151", None), ("srf_k101", None), ("competitive_mapping", None)):
            out, summary=evaluate(reference, estimate, method); all_rows += [{"species":species, **x} for x in out]; summaries.append({"species":species, **summary})
    fields=sorted({k for x in all_rows for k in x})
    with (args.outdir/"family_prioritization.tsv").open("w", newline="") as h: w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(all_rows)
    with (args.outdir/"summary.json").open("w") as h: json.dump({"endpoint":"historical_assembly_under_representation_prioritization","threshold":THRESHOLD,"reference_proxy_label":"frozen old/new ratio below 0.6 among predeclared >=15kb-new-assembly eligible families","competitive_mapping_discovery":"N/A_shared_TandemX_catalogue","formal_native_comparator_status":"blocked_t7_unverified","interpretation_boundary":"TandemX values are frozen bookkeeping against the same retrospective newer-assembly proxy used to define the label; SRF k151, SRF k101 and competitive mapping remain separate N/A native-comparator rows until their inputs run.","summaries":summaries},h,indent=2); h.write("\n")


if __name__ == "__main__": main()
