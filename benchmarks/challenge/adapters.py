"""Strict normalization of tool-specific output to 0-based intervals."""

from __future__ import annotations

from pathlib import Path

from .schema import ArrayRecord, read_table


def parse_tandemx(path: Path) -> list[ArrayRecord]:
    # Candidate TSV does not contain its consensus sequence. Family recovery is
    # evaluated separately using the emitted monomers.fa catalog.
    return [ArrayRecord(r["read_id"], int(r["read_start"]), int(r["read_end"]), int(r["period_bp"]))
            for r in read_table(path, {"read_id", "read_start", "read_end", "period_bp"})]


def parse_trf(path: Path) -> list[ArrayRecord]:
    arrays: list[ArrayRecord] = []
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
            arrays.append(ArrayRecord(read_id, int(fields[0]) - 1, int(fields[1]), int(fields[2]), fields[13]))
    return arrays


def parse_tidehunter(path: Path) -> list[ArrayRecord]:
    arrays: list[ArrayRecord] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 11:
                raise ValueError(f"Expected TideHunter -f 2 (11 fields): {path}:{line_number}")
            arrays.append(ArrayRecord(fields[0], int(fields[4]) - 1, int(fields[5]), int(fields[6]), fields[10]))
    return arrays


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
    parser = {"tandemx": parse_tandemx, "trf": parse_trf, "tidehunter": parse_tidehunter}[tool]
    return [r for r in parser(path) if min_period <= r.period <= max_period and r.end - r.start >= min_span]


def build_command(tool: str, executable: str, reads: Path, outdir: Path,
                  min_period: int, max_period: int, min_span: int) -> tuple[list[str], Path]:
    if tool == "tandemx":
        output = outdir / "discover" / "candidate_reads.tsv"
        return [executable, "discover", "--reads", str(reads), "--outdir", str(output.parent),
                "--min-period", str(min_period), "--max-period", str(max_period),
                "--min-repeat-span", str(min_span), "--min-support-reads", "1", "--min-read-length", "1",
                "--kmer-backend", "rust", "--threads", "1", "--no-progress"], output
    if tool == "trf":
        return [executable, str(reads), "2", "7", "7", "80", "10", "50", str(max_period), "-ngs", "-h"], outdir / "trf.txt"
    if tool == "tidehunter":
        output = outdir / "tidehunter.tsv"
        return [executable, "-t", "1", "-f", "2", "-p", str(min_period), "-P", str(max_period),
                "-m", str(min_period), "-c", "2", "-o", str(output), str(reads)], output
    raise ValueError(f"Unknown tool: {tool}")
