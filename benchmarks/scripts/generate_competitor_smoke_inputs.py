"""Generate tiny deterministic comparator availability inputs.

These are synthetic interface fixtures, never biological truth or a formal
cross-tool benchmark. Usage: python generate_competitor_smoke_inputs.py OUTDIR
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path


MONOMER_A = "ACGTGCAATGTCAGTACCGTACGATCGTTA"
MONOMER_B = "TTGACCGATGCTAGACCTGAGTCATCGTAC"


def generate(outdir: Path) -> None:
    assert len(MONOMER_A) == len(MONOMER_B) == 30
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "trf_tidehunter.fa").write_text(
        ">synthetic_repeat\n" + MONOMER_A * 12 + "\n"
    )
    (outdir / "cendetecthor_ab.fa").write_text(
        ">toy.chr1\n" + (MONOMER_A + MONOMER_B) * 12 + "\n"
    )
    # The pipeline's periodicity script parses chromosome/start/end from
    # a FASTA header of the form >chrom:start-end. A periodic A insertion
    # makes the toy's k-mer spacing nonconstant for its window filter.
    pipeline_sequence = ((MONOMER_A + MONOMER_B) * 9 + MONOMER_A) * 10
    (outdir / "cendetecthor_pipeline.fa").write_text(
        f">toy:0-{len(pipeline_sequence)}\n" + pipeline_sequence + "\n"
    )
    # SRF consumes KMC-style canonical k-mer/count text. This deterministic
    # core smoke bypasses KMC; it is not a whole read-first SRF workflow.
    complement = str.maketrans("ACGT", "TGCA")
    counts = Counter()
    for start in range(len(pipeline_sequence) - 151 + 1):
        kmer = pipeline_sequence[start : start + 151]
        counts[min(kmer, kmer.translate(complement)[::-1])] += 1
    (outdir / "srf_k151_counts.txt").write_text(
        "".join(f"{kmer}\t{count}\n" for kmer, count in sorted(counts.items()))
    )
    (outdir / "cendetecthor_ab.bed").write_text(
        "".join(
            f"toy.chr1\t{i * 30}\t{(i + 1) * 30}\tmonomer_{i + 1:02d}\t0\t+\n"
            for i in range(24)
        )
    )
    for copies in (12, 100):
        perfect_sequence = (MONOMER_A + MONOMER_B) * copies
        (outdir / f"cendetecthor_pipeline_perfect_{len(perfect_sequence)}.fa").write_text(
            f">toy:0-{len(perfect_sequence)}\n" + perfect_sequence + "\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    generate(args.outdir)


if __name__ == "__main__":
    main()
