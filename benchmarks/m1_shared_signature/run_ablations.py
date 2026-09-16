"""Bounded A/C/D development tournament on the 12 frozen M1 toy cases.

Example: conda run -n tandemx-dev python -m \
  benchmarks.m1_shared_signature.run_ablations --outdir /tmp/tandemx-m1-acd
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import time

from . import ablations, model, probabilistic, run as original, run_probabilistic


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def source_label_confusion(labels: list[str], assignments: list[str]) -> dict[str, int]:
    return dict(sorted(Counter(f"{truth}->{call}" for truth, call in zip(labels, assignments)).items()))


def evaluate_case(seed: int, scenario: tuple, b_prior: dict, e_prior: dict) -> dict:
    catalogue, reads, labels = original.make_case(seed, *scenario[1:])
    for prior in (b_prior, e_prior):
        if prior["catalogue_sha256"] != digest(catalogue) or prior["reads_sha256"] != digest(reads):
            raise ValueError("Frozen M1 input hash changed")
    gate = round(original.PARAMETERS["unit_bp"] * original.PARAMETERS["gate_fraction"])
    start = time.perf_counter()
    accepted, mapping, rejected, ties = model.align_and_gate(reads, catalogue, gate)
    mapping_seconds = time.perf_counter() - start
    if mapping != b_prior["ordinary_mapping"] or mapping != e_prior["ordinary_mapping"]:
        raise ValueError("Ordinary mapping changed since frozen evidence")
    accepted_counts = Counter(accepted)
    accepted_labels = []
    for read, label in zip(reads, labels):
        if accepted_counts[read]:
            accepted_labels.append(label)
            accepted_counts[read] -= 1
    if len(accepted_labels) != len(accepted):
        raise ValueError("Lost accepted read provenance")

    start = time.perf_counter()
    b_family, b_group, _ = model.shared_signature_fit(accepted, catalogue, ablations.K)
    b_seconds = time.perf_counter() - start
    if b_family != b_prior["shared_signature"] or b_family != e_prior["shared_signature"]:
        raise ValueError("B shared-signature estimate changed")
    start = time.perf_counter()
    e = probabilistic.infer(reads, catalogue)
    e_seconds = time.perf_counter() - start
    e.pop("read_assignment")
    if e["family"] != e_prior["probabilistic"]["family"]:
        raise ValueError("E per-read estimate changed")

    outputs: dict[str, dict] = {}
    seconds: dict[str, float] = {"ordinary_mapping": mapping_seconds,
                                 "B_shared_signature": b_seconds,
                                 "E_read_probability": e_seconds}
    for name, implementation in (("A_discriminative", ablations.discriminative),
                                 ("C_poisson", ablations.poisson_counts),
                                 ("D_kmer_em", ablations.em_mixture)):
        start = time.perf_counter()
        result = implementation(accepted, catalogue)
        seconds[name] = time.perf_counter() - start
        assignments = result.pop("read_assignment", None)
        if assignments is not None:
            result["source_to_call"] = source_label_confusion(accepted_labels, assignments)
            result["wrong_family_attribution"] = sum(
                call not in (truth, "ambiguous", "f1+f2")
                for truth, call in zip(accepted_labels, assignments) if truth != "unknown")
            result["negative_false_attribution"] = sum(
                call != "ambiguous" for truth, call in zip(accepted_labels, assignments)
                if truth == "unknown")
        result["score"] = run_probabilistic.score(result["family"])
        result["total_rejected_reads"] = rejected + result.get("rejected_ambiguous", 0)
        outputs[name] = result
    return dict(seed=seed, scenario=scenario[0], reads_sha256=digest(reads),
                catalogue_sha256=digest(catalogue), truth_reads=run_probabilistic.TRUTH,
                n_reads=len(reads), gate_rejected=rejected, mapping_tied=ties,
                ordinary_mapping=dict(family=mapping, score=run_probabilistic.score(mapping),
                                      total_rejected_reads=rejected + ties),
                B_shared_signature=dict(family=b_family, groups=b_group,
                                        score=run_probabilistic.score(b_family),
                                        total_rejected_reads=rejected),
                E_read_probability=dict(family=e["family"], groups=e["groups"],
                                        score=run_probabilistic.score(e["family"]),
                                        total_rejected_reads=e["rejected_gate"]
                                        + e["rejected_unknown"] + e["rejected_ambiguous"],
                                        source_to_call=e_prior["probabilistic_confusion"]),
                F_prior_read_bootstrap=dict(B_coverage=b_prior["scores"]["shared_signature"]["positive_interval_covered"],
                                            E_coverage=e_prior["bootstrap_positive_coverage"],
                                            note="40 read resamples; no new intervals for A/C/D"),
                methods=outputs, seconds=seconds)


def run(outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    directory = Path(__file__).parent
    b_path = directory / "evidence_20260916" / "results.json"
    e_path = directory / "evidence_probabilistic_20260916" / "results.json"
    b_prior, e_prior = (json.loads(path.read_text()) for path in (b_path, e_path))
    source_map = {"run.py": original, "model.py": model,
                  "probabilistic.py": probabilistic,
                  "run_probabilistic.py": run_probabilistic}
    for filename, module in source_map.items():
        expected = (b_prior if filename in ("run.py", "model.py") else e_prior)["source_sha256"][filename]
        if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Frozen source changed: {filename}")
    b_cases = {(case["scenario"], case["seed"]): case for case in b_prior["cases"]}
    e_cases = {(case["scenario"], case["seed"]): case for case in e_prior["cases"]}
    cases = [evaluate_case(seed, scenario, b_cases[(scenario[0], seed)],
                           e_cases[(scenario[0], seed)])
             for scenario in original.SCENARIOS for seed in original.PARAMETERS["seeds"]]
    source_sha256 = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                     for name in ("run_ablations.py", "ablations.py", *source_map)}
    summary = dict(status="development_only_no_heldout", cases=cases,
                   scenarios=original.SCENARIOS, parameters=original.PARAMETERS,
                   A_C_D_settings=dict(k=ablations.K, assumed_substitution=ablations.ASSUMED_SUBSTITUTION,
                                       max_iterations=ablations.MAX_ITERATIONS,
                                       posterior_min=ablations.POSTERIOR_MIN),
                   prior_evidence_sha256={path.name + ":" + path.parent.name:
                                          hashlib.sha256(path.read_bytes()).hexdigest()
                                          for path in (b_path, e_path)},
                   source_sha256=source_sha256,
                   process_max_rss_raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   resource_note="macOS ru_maxrss bytes; one process; fit times include design building; no T7 I/O")
    (outdir / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    run(parser.parse_args().outdir)
