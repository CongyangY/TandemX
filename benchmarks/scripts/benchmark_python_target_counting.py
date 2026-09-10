"""Reproducibly compare legacy and rolling-code Python quantify target counts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
import time


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def worker(mode: str, reads: Path, catalogue: Path, genome_size: int, k: int) -> dict:
    from tandemx.quantify import mvp
    if mode == "legacy":
        def legacy(counts, sequence, word_k, targets, _code_targets):
            counts.update(word for word in mvp.iter_kmers(sequence, word_k) if word in targets)
        mvp.update_selected_python_kmer_counts = legacy
    root = Path(tempfile.mkdtemp(prefix="tandemx-python-target-count-"))
    try:
        start = time.perf_counter()
        mvp.quantify_toy_copy_number(
            mvp.QuantifyConfig(reads=reads, monomers=catalogue, genome_size=genome_size,
                                outdir=root / "out", k=k, haploid_depth=1.0,
                                kmer_backend="python")
        )
        elapsed = time.perf_counter() - start
        output = root / "out" / "copy_number.tsv"
        return {"mode": mode, "quantify_call_seconds": elapsed,
                "ru_maxrss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                "copy_number_sha256": digest(output)}
    finally:
        shutil.rmtree(root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("legacy", "rolling"))
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--catalogue", type=Path, required=True)
    parser.add_argument("--genome-size", type=int, required=True)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--result-json", type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker, args.reads, args.catalogue, args.genome_size, args.k)))
        return
    if args.repetitions < 1 or args.result_json is None:
        parser.error("parent mode requires positive --repetitions and --result-json")
    command_base = [sys.executable, str(Path(__file__).resolve()), "--reads", str(args.reads.resolve()),
                    "--catalogue", str(args.catalogue.resolve()), "--genome-size", str(args.genome_size),
                    "--k", str(args.k)]
    rows = []
    for mode in ("legacy", "rolling"):
        for _ in range(args.repetitions):
            run = subprocess.run(command_base + ["--worker", mode], text=True, capture_output=True, check=True)
            rows.append(json.loads(run.stdout))
    if len({row["copy_number_sha256"] for row in rows}) != 1:
        raise RuntimeError("Legacy and rolling outputs differ")
    result = {
        "schema_version": 1,
        "complete": True,
        "timer_scope": "time.perf_counter around quantify_toy_copy_number only; interpreter startup and imports excluded",
        "reads": {"path": str(args.reads.resolve()), "sha256": digest(args.reads)},
        "catalogue": {"path": str(args.catalogue.resolve()), "sha256": digest(args.catalogue)},
        "genome_size": args.genome_size, "k": args.k, "repetitions": args.repetitions,
        "python": sys.version, "platform": platform.platform(),
        "benchmark_script_sha256": digest(Path(__file__)), "measurements": rows,
        "boundary": "small consumed development input; not a release performance claim",
    }
    args.result_json.parent.mkdir(parents=True, exist_ok=True)
    if args.result_json.exists():
        raise FileExistsError(args.result_json)
    args.result_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
