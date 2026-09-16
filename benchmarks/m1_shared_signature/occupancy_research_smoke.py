"""Same-input development smoke for the experimental streaming occupancy backend.

One predeclared seed from each frozen M1 scenario. This is a feasibility and
negative-result screen, not a held-out accuracy benchmark.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import time

from tandemx.io.sequences import SequenceRecord
from .occupancy_research import CompetitiveConfig, classify_window, summarize_records

from . import run as original


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def run(outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=False)
    prior_path = Path(__file__).with_name("evidence_20260916") / "results.json"
    prior = json.loads(prior_path.read_text())
    prior_cases = {(case["scenario"], case["seed"]): case for case in prior["cases"]}
    rows = []
    for scenario in original.SCENARIOS:
        for seed in original.PARAMETERS["seeds"]:
            catalogue, reads, labels = original.make_case(seed, *scenario[1:])
            previous = prior_cases[(scenario[0], seed)]
            if digest(reads) != previous["reads_sha256"] or digest(catalogue) != previous["catalogue_sha256"]:
                raise ValueError("Frozen M1 development input changed")
            config = CompetitiveConfig(reads=Path("not_read"), monomers=Path("not_read"),
                                       genome_size=48_000, outdir=outdir, haploid_depth=1.0)
            records = (SequenceRecord(id=str(index), sequence=read)
                       for index, read in enumerate(reads))
            start = time.perf_counter()
            summary = summarize_records(records, catalogue, config)
            elapsed = time.perf_counter() - start
            scoring_start = time.perf_counter()
            calls = [classify_window(read, catalogue, config) for read in reads]
            confusion = dict(sorted(Counter(
                f"{source}->{identity if outcome == 'assigned' else outcome}"
                for source, (outcome, identity) in zip(labels, calls)
            ).items()))
            wrong_family_calls = sum(
                outcome == "assigned" and identity != source
                for source, (outcome, identity) in zip(labels, calls) if source != "unknown"
            )
            negative_false_calls = sum(
                outcome == "assigned" for source, (outcome, _identity) in zip(labels, calls)
                if source == "unknown"
            )
            scoring_seconds = time.perf_counter() - scoring_start
            candidate = {name: summary["assigned_by_family"][name] / 80
                         for name in catalogue}
            truth = {"f1": 320, "f2": 180, "decoy_zero": 0}
            mare = sum(abs(candidate[name] / truth[name] - 1) for name in ("f1", "f2")) / 2
            rows.append(dict(scenario=scenario[0], seed=seed,
                             reads_sha256=digest(reads), catalogue_sha256=digest(catalogue),
                             truth_unit_sized_reads=truth,
                             ordinary_mapping=previous["ordinary_mapping"],
                             ordinary_positive_mare=previous["scores"]["ordinary_mapping"]["all_positive_mare_refusal_penalized"],
                             occupancy_read_equivalents=candidate,
                             occupancy_positive_mare=mare,
                             assigned_bp=summary["assigned_read_bp"],
                             ambiguous_bp=summary["ambiguous_read_bp"],
                             unknown_bp=summary["unknown_read_bp"],
                             total_bp=summary["total_read_bp"],
                             source_to_call=confusion,
                             wrong_family_calls=wrong_family_calls,
                             negative_false_calls=negative_false_calls,
                             wall_seconds=elapsed,
                             truth_scoring_seconds=scoring_seconds))
    source_files = [Path(__file__), Path(original.__file__),
                    Path(__file__).with_name("occupancy_research.py")]
    result = dict(status="development_smoke_not_validation", seeds=original.PARAMETERS["seeds"], rows=rows,
                  prior_results_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest(),
                  source_sha256={str(path.relative_to(Path(__file__).parents[2])):
                                 hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in source_files},
                  process_max_rss_raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  resource_note="macOS ru_maxrss bytes; one Python process; includes simulation and benchmarking; no T7 I/O")
    (outdir / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    run(parser.parse_args().outdir)
