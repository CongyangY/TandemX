"""Audit bounded native-read alignments against a frozen array context.

The local alignment is only a screening result. Its MAPQ cannot establish
genome-wide locus specificity or a physical array copy count.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path


_CIGAR = re.compile(r"(\d+)([MID=X])")
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith(">") or any(x.startswith(">") for x in lines[1:]):
        raise ValueError("expected one FASTA record")
    sequence = "".join(lines[1:]).upper()
    if not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("invalid FASTA sequence")
    return sequence


def read_fastq(path: Path) -> dict[str, str]:
    reads: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="ascii") as handle:
        while header := handle.readline():
            sequence = handle.readline().strip().upper()
            separator = handle.readline()
            quality = handle.readline().strip()
            if not header.startswith("@") or not separator.startswith("+") or len(sequence) != len(quality):
                raise ValueError("invalid selected FASTQ")
            read_id = header[1:].split()[0]
            if read_id in reads:
                raise ValueError("duplicate selected read ID")
            reads[read_id] = sequence
    return reads


def parse_cigar(text: str) -> list[tuple[int, str]]:
    entries = [(int(n), op) for n, op in _CIGAR.findall(text)]
    if not entries or "".join(f"{n}{op}" for n, op in entries) != text:
        raise ValueError("invalid CIGAR")
    return entries


def partition_identity(
    reference: str, query: str, strand: str, query_start: int, query_end: int,
    target_start: int, cigar: str, boundaries: tuple[int, int],
) -> list[dict[str, int | float | None]]:
    """Count exact pairwise matches and alignment columns in three target bins."""
    if strand not in ("+", "-") or not (0 <= query_start <= query_end <= len(query)):
        raise ValueError("invalid query coordinates")
    if not (0 < boundaries[0] < boundaries[1] < len(reference)):
        raise ValueError("invalid target boundaries")
    oriented = query if strand == "+" else query.translate(_RC)[::-1]
    query_position = query_start if strand == "+" else len(query) - query_end
    target_position = target_start
    counts = [[0, 0] for _ in range(3)]

    def bin_index(position: int) -> int:
        return 0 if position < boundaries[0] else 1 if position < boundaries[1] else 2

    for length, operation in parse_cigar(cigar):
        if operation in ("M", "=", "X"):
            for _ in range(length):
                if query_position >= len(oriented) or target_position >= len(reference):
                    raise ValueError("CIGAR extends past sequence")
                bucket = counts[bin_index(target_position)]
                bucket[1] += 1
                bucket[0] += oriented[query_position] == reference[target_position]
                query_position += 1
                target_position += 1
        elif operation == "I":
            if query_position + length > len(oriented):
                raise ValueError("CIGAR extends past query")
            counts[bin_index(min(target_position, len(reference) - 1))][1] += length
            query_position += length
        else:
            if target_position + length > len(reference):
                raise ValueError("CIGAR extends past target")
            for _ in range(length):
                counts[bin_index(target_position)][1] += 1
                target_position += 1
    expected_query_end = query_end if strand == "+" else len(query) - query_start
    if query_position != expected_query_end:
        raise ValueError("CIGAR/query span mismatch")
    return [{"matches": m, "columns": c, "identity": round(m / c, 6) if c else None}
            for m, c in counts]


def audit(paf: Path, context: Path, selected_fastq: Path,
          array_start: int, array_end: int, flank_bp: int = 1000,
          minimum_identity: float = 0.95) -> dict:
    reference = read_fasta(context)
    reads = read_fastq(selected_fastq)
    if not 0 <= minimum_identity <= 1:
        raise ValueError("minimum identity must be between zero and one")
    if not (0 <= array_start - flank_bp < array_start < array_end < array_end + flank_bp <= len(reference)):
        raise ValueError("array or flank outside context")
    paf_rows = 0
    distinct_queries: set[str] = set()
    primary_array_reads: set[str] = set()
    geometric: list[dict] = []
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError("invalid PAF")
            paf_rows += 1
            read_id = fields[0]
            distinct_queries.add(read_id)
            query_length, query_start, query_end = map(int, fields[1:4])
            strand = fields[4]
            target_length, target_start, target_end = map(int, fields[6:9])
            matches, columns, mapq = map(int, fields[9:12])
            tags = {field[:2]: field[5:] for field in fields[12:] if len(field) >= 6 and field[2:5] == ":A:"}
            if target_length != len(reference) or tags.get("tp") != "P":
                continue
            if target_start <= array_start and target_end >= array_end:
                primary_array_reads.add(read_id)
            if target_start > array_start - flank_bp or target_end < array_end + flank_bp:
                continue
            if read_id not in reads or len(reads[read_id]) != query_length:
                raise ValueError(f"missing diagnostic read: {read_id}")
            cigar_tags = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
            if len(cigar_tags) != 1:
                raise ValueError("missing CIGAR for diagnostic alignment")
            partitions = partition_identity(reference, reads[read_id], strand,
                                            query_start, query_end, target_start,
                                            cigar_tags[0], (array_start, array_end))
            if sum(p["matches"] for p in partitions) != matches or sum(p["columns"] for p in partitions) != columns:
                raise ValueError("CIGAR-derived identity disagrees with PAF")
            geometric.append({"read_id": read_id, "query_length": query_length,
                              "strand": strand, "query_start": query_start,
                              "query_end": query_end, "target_start": target_start,
                              "target_end": target_end, "mapq_local": mapq,
                              "identity": round(matches / columns, 6),
                              "left_flank": partitions[0], "array": partitions[1],
                              "right_flank": partitions[2]})
    qualified = sorted({row["read_id"] for row in geometric if row["identity"] >= minimum_identity})
    return {"schema_version": 1, "context_sha256": sha256_file(context),
            "paf_sha256": sha256_file(paf), "diagnostic_reads_sha256": sha256_file(selected_fastq),
            "array_context_start_0based": array_start, "array_context_end_0based": array_end,
            "minimum_natural_flank_bp_each_side": flank_bp,
            "minimum_whole_alignment_identity": minimum_identity,
            "paf_alignment_count": paf_rows, "paf_distinct_query_count": len(distinct_queries),
            "primary_array_overlap_read_count": len(primary_array_reads),
            "geometric_spanner_count": len({row["read_id"] for row in geometric}),
            "qualified_spanner_count": len(qualified), "qualified_read_ids": qualified,
            "geometric_alignments": sorted(geometric, key=lambda row: row["read_id"])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--diagnostic-reads", type=Path, required=True)
    parser.add_argument("--array-start", type=int, required=True)
    parser.add_argument("--array-end", type=int, required=True)
    parser.add_argument("--flank-bp", type=int, default=1000)
    parser.add_argument("--minimum-identity", type=float, default=0.95)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.paf, args.context, args.diagnostic_reads,
                   args.array_start, args.array_end, args.flank_bp,
                   args.minimum_identity)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
