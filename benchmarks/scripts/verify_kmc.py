"""Check a built KMC against independent exact canonical counts on bounded data."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import random

from benchmarks.challenge.run import run_process
from benchmarks.challenge.schema import digest_file


def verify(kmc: Path, dump: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=False)
    rng = random.Random(981782)
    a = "".join(rng.choices("ACGT", k=1200))
    reverse = a.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    reads = [a, a, reverse, a[:800] + "N" * 10 + a[810:], "AT" * 400]
    fasta = outdir / "reads.fa"
    fasta.write_text("".join(f">r{i}\n{read[:600]}\n{read[600:]}\n" for i, read in enumerate(reads)))
    results = []
    for k in (17, 151):
        expected = Counter()
        for read in reads:
            for i in range(len(read) - k + 1):
                word = read[i:i + k]
                if "N" not in word:
                    expected[min(word, word.translate(str.maketrans("ACGT", "TGCA"))[::-1])] += 1
        folder = outdir / str(k)
        folder.mkdir()
        temporary = folder / "tmp"
        temporary.mkdir()
        commands = [[str(kmc), "-fm", f"-k{k}", "-t1", "-m2", "-sm", "-ci1", "-cs100000",
                     str(fasta), str(folder / "counts"), str(temporary)],
                    [str(dump), str(folder / "counts"), str(folder / "counts.tsv")]]
        for i, command in enumerate(commands):
            measured = run_process(command, folder / f"stage{i}.stdout", folder / f"stage{i}.stderr", 120)
            results.append({"k": k, "command": command, **measured})
            if measured["exit_code"]:
                (outdir / "failure.json").write_text(json.dumps(results, indent=2) + "\n")
                raise ValueError(f"KMC validation subprocess failed; see {folder}")
        observed = {}
        for row in (folder / "counts.tsv").read_text().splitlines():
            word, count = row.split()
            if word in observed:
                raise ValueError("Duplicate dump row")
            observed[word] = int(count)
        if observed != expected:
            raise ValueError(f"KMC differs from independent Counter at k={k}")
        results[-1].update(exact_count_agreement=True, distinct_kmers=len(expected), total_kmers=sum(expected.values()))
    receipt = {"results": results, "input_sha256": digest_file(fasta),
               "executable_sha256": {str(p): digest_file(p) for p in (kmc, dump)},
               "script_sha256": digest_file(Path(__file__)), "complete": True}
    (outdir / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("kmc", "dump", "outdir"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    verify(args.kmc, args.dump, args.outdir)
