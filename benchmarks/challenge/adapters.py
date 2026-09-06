"""Strict normalization of tool-specific output to 0-based intervals."""

from __future__ import annotations

from pathlib import Path
from collections.abc import Iterator

from .schema import ArrayRecord, iter_table


def iter_tandemx(path: Path) -> Iterator[ArrayRecord]:
    # Candidate TSV does not contain its consensus sequence. Family recovery is
    # evaluated separately using the emitted monomers.fa catalog.
    for r in iter_table(path, {"read_id", "read_start", "read_end", "period_bp"}):
        yield ArrayRecord(r["read_id"], int(r["read_start"]), int(r["read_end"]), int(r["period_bp"]))


def iter_trf(path: Path) -> Iterator[ArrayRecord]:
    read_id = ""
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("@"):
                read_id = line[1:].split()[0]
                continue
            fields = line.split()
            if not read_id or len(fields) < 14:
                raise ValueError(f"Malformed TRF -ngs row {path}:{line_number}")
            yield ArrayRecord(read_id, int(fields[0]) - 1, int(fields[1]), int(fields[2]), fields[13])


def iter_tidehunter(path: Path) -> Iterator[ArrayRecord]:
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 11:
                raise ValueError(f"Expected TideHunter -f 2 (11 fields): {path}:{line_number}")
            yield ArrayRecord(fields[0], int(fields[4]) - 1, int(fields[5]), int(fields[6]), fields[10])


def iter_ultra(path: Path) -> Iterator[ArrayRecord]:
    """ULTRA 1.2.2 TSV: start/end are already zero-based half-open.

    Verified against upstream TabFileWriter.cpp (start + repeatLength). An
    unresolved '*' is represented by N; '.' denotes a missing consensus. Neither
    can be counted as a correct known base in family recovery.
    """
    for r in iter_table(path, {"SeqID", "Start", "End", "Period", "Consensus"}):
        yield ArrayRecord(r["SeqID"], int(r["Start"]), int(r["End"]), int(r["Period"]),
                          "" if r["Consensus"] == "." else r["Consensus"].upper().replace("*", "N"))


def parse_tandemx(path: Path) -> list[ArrayRecord]:
    return list(iter_tandemx(path))


def parse_trf(path: Path) -> list[ArrayRecord]:
    return list(iter_trf(path))


def parse_tidehunter(path: Path) -> list[ArrayRecord]:
    return list(iter_tidehunter(path))


def parse_ultra(path: Path) -> list[ArrayRecord]:
    return list(iter_ultra(path))


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in sequences:
                    raise ValueError(f"Duplicate FASTA identifier in {path}: {name}")
                sequences[name] = ""
            elif not name:
                raise ValueError(f"Sequence before FASTA header: {path}")
            else:
                sequences[name] += line.upper()
    if any(not sequence for sequence in sequences.values()):
        raise ValueError(f"Empty FASTA record: {path}")
    return sequences


def parse_arrays(tool: str, path: Path, min_period: int, max_period: int, min_span: int) -> list[ArrayRecord]:
    return list(iter_arrays(tool, path, min_period, max_period, min_span))


def iter_arrays(tool: str, path: Path, min_period: int, max_period: int, min_span: int) -> Iterator[ArrayRecord]:
    parser = {"tandemx": iter_tandemx, "trf": iter_trf, "tidehunter": iter_tidehunter, "ultra": iter_ultra}[tool]
    for row in parser(path):
        if min_period <= row.period <= max_period and row.end - row.start >= min_span:
            yield row


def build_command(tool: str, executable: str, reads: Path, outdir: Path,
                  min_period: int, max_period: int, min_span: int,
                  threads: int = 1) -> tuple[list[str], Path]:
    if not isinstance(threads, int) or not 1 <= threads <= 64:
        raise ValueError("Threads must be an integer in [1,64]")
    if tool == "tandemx":
        output = outdir / "discover" / "candidate_reads.tsv"
        return [executable, "discover", "--reads", str(reads), "--outdir", str(output.parent),
                "--min-period", str(min_period), "--max-period", str(max_period),
                "--min-repeat-span", str(min_span), "--min-support-reads", "1", "--min-read-length", "1",
                "--kmer-backend", "rust", "--threads", str(threads), "--no-progress"], output
    if tool == "trf":
        return [executable, str(reads), "2", "7", "7", "80", "10", "50", str(max_period), "-ngs", "-h"], outdir / "trf.txt"
    if tool == "tidehunter":
        output = outdir / "tidehunter.tsv"
        return [executable, "-t", str(threads), "-f", "2", "-p", str(min_period), "-P", str(max_period),
                "-m", str(min_period), "-c", "2", "-o", str(output), str(reads)], output
    if tool == "ultra":
        output = outdir / "ultra.tsv"
        return [executable, "-t", str(threads), "-p", str(max_period), "--min_length", str(min_span),
                "--max_consensus", str(max_period), "--tsv", "-o", str(output), str(reads)], output
    raise ValueError(f"Unknown tool: {tool}")
