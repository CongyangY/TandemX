"""Command-line interface skeleton for TandemX."""

from __future__ import annotations

import argparse
import logging
import os
import platform as platform_module
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from tandemx import __version__
from tandemx.simulate.toy import ToySimulationConfig, generate_toy_dataset, parse_int_list


COMMANDS = ("discover", "quantify", "locate", "probe", "compare", "visualize", "simulate")


class InputFileError(Exception):
    """Raised when a required input path is missing or not a regular file."""


def _path_value(value: str) -> Path:
    return Path(value)


def _require_existing_files(paths: Iterable[Path]) -> None:
    for path in paths:
        if not path.exists():
            raise InputFileError(f"Input file does not exist: {path}")
        if not path.is_file():
            raise InputFileError(f"Input path is not a file: {path}")


def _yaml_value(value: Any, indent: int = 0) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, Path):
        value = str(value)
    if isinstance(value, (list, tuple)):
        if not value:
            return "[]"
        prefix = " " * indent
        nested = []
        for item in value:
            nested.append(f"{prefix}- {_yaml_value(item, indent + 2)}")
        return "\n".join(nested)
    text = str(value)
    escaped = text.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _write_run_config(
    outdir: Path,
    command: str,
    args: argparse.Namespace,
    status: str = "skeleton_not_implemented",
) -> None:
    config_path = outdir / "run_config.yaml"
    values = {
        key: value
        for key, value in vars(args).items()
        if key != "func" and not key.startswith("_")
    }
    timestamp = datetime.now(timezone.utc).isoformat()
    lines = [
        f"command: {_yaml_value(f'tandemx {command}')}",
        f"subcommand: {_yaml_value(command)}",
        f"version: {_yaml_value(__version__)}",
        f"timestamp_utc: {_yaml_value(timestamp)}",
        f"cwd: {_yaml_value(os.getcwd())}",
        "argv:",
        _yaml_value(getattr(args, "_argv", []), indent=2),
        f"python_version: {_yaml_value(sys.version.split()[0])}",
        f"platform: {_yaml_value(platform_module.platform())}",
        f"status: {_yaml_value(status)}",
        "parameters:",
    ]
    for key in sorted(values):
        value = values[key]
        lines.append(f"  {key}: {_yaml_value(value, indent=4)}")
    config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _configure_log(outdir: Path, command: str) -> logging.Logger:
    logger = logging.getLogger(f"tandemx.{command}")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    handler = logging.FileHandler(outdir / "run.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    return logger


def _prepare_placeholder_run(
    command: str,
    args: argparse.Namespace,
    required_inputs: Iterable[Path],
) -> int:
    _require_existing_files(required_inputs)
    args.outdir.mkdir(parents=True, exist_ok=True)
    _write_run_config(args.outdir, command, args)
    logger = _configure_log(args.outdir, command)
    logger.info("command=tandemx %s", command)
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=skeleton_not_implemented")
    logger.info("Core algorithm is not implemented yet")

    print(
        f"tandemx {command}: not implemented yet. "
        f"Wrote run_config.yaml and run.log to {args.outdir}"
    )
    return 0


def run_discover(args: argparse.Namespace) -> int:
    return _prepare_placeholder_run("discover", args, [args.reads])


def run_quantify(args: argparse.Namespace) -> int:
    return _prepare_placeholder_run(
        "quantify",
        args,
        [args.reads, args.catalogue, args.monomers],
    )


def run_locate(args: argparse.Namespace) -> int:
    return _prepare_placeholder_run(
        "locate",
        args,
        [args.assembly, args.catalogue, args.monomers],
    )


def run_probe(args: argparse.Namespace) -> int:
    required = [args.catalogue, args.monomers, args.copy_number]
    if args.locations is not None:
        required.append(args.locations)
    return _prepare_placeholder_run("probe", args, required)


def run_compare(args: argparse.Namespace) -> int:
    return _prepare_placeholder_run(
        "compare",
        args,
        [args.read_copy_number, args.assembly_density],
    )


def run_visualize(args: argparse.Namespace) -> int:
    required = [args.catalogue, args.copy_number]
    for optional in (args.locations, args.probes, args.comparison):
        if optional is not None:
            required.append(optional)
    return _prepare_placeholder_run("visualize", args, required)


def run_simulate_toy(args: argparse.Namespace) -> int:
    monomer_lengths = parse_int_list(args.monomer_lengths, "--monomer-lengths")
    copies = parse_int_list(args.copies, "--copies")
    config = ToySimulationConfig(
        outdir=args.outdir,
        seed=args.seed,
        num_reads=args.num_reads,
        read_length=args.read_length,
        background_length=args.background_length,
        monomer_lengths=monomer_lengths,
        copies=copies,
        error_rate=args.error_rate,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    generate_toy_dataset(config)
    _write_run_config(args.outdir, "simulate toy", args, status="toy_dataset_generated")
    logger = _configure_log(args.outdir, "simulate.toy")
    logger.info("command=tandemx simulate toy")
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=toy_dataset_generated")
    logger.info("Generated reproducible toy dataset")
    print(f"tandemx simulate toy: wrote toy dataset to {args.outdir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tandemx",
        description="TandemX CLI skeleton for toy-scale tandem repeat analysis workflows.",
    )
    parser.add_argument("--version", action="version", version=f"tandemx {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    discover = subparsers.add_parser(
        "discover",
        help="Discover candidate tandem repeat monomers from reads.",
    )
    discover.add_argument("--reads", required=True, type=_path_value, help="Input HiFi-like reads in FASTA or FASTQ format.")
    discover.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and future discover outputs.")
    discover.add_argument("--min-monomer-len", type=int, default=20, help="Minimum candidate monomer length in bp.")
    discover.add_argument("--max-monomer-len", type=int, default=2000, help="Maximum candidate monomer length in bp.")
    discover.add_argument("--min-support-reads", type=int, default=5, help="Minimum number of reads supporting a candidate family.")
    discover.add_argument("--min-repeat-span", type=int, default=100, help="Minimum repeat-supporting span in a read, in bp.")
    discover.add_argument("--seed", type=int, default=1, help="Random seed for reproducible toy-scale workflows.")
    discover.set_defaults(func=run_discover)

    quantify = subparsers.add_parser(
        "quantify",
        help="Estimate read-based repeat copy number from diagnostic k-mers.",
    )
    quantify.add_argument("--reads", required=True, type=_path_value, help="Input reads used to estimate diagnostic k-mer depth.")
    quantify.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat family catalogue produced by discover.")
    quantify.add_argument("--monomers", required=True, type=_path_value, help="Monomer FASTA produced by discover.")
    quantify.add_argument("--genome-size", required=True, type=int, help="Estimated haploid or target genome size in bp.")
    quantify.add_argument("--k", type=int, default=21, help="Diagnostic k-mer size.")
    quantify.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and future quantify outputs.")
    quantify.set_defaults(func=run_quantify)

    locate = subparsers.add_parser(
        "locate",
        help="Localize repeat family evidence on an assembly.",
    )
    locate.add_argument("--assembly", required=True, type=_path_value, help="Genome assembly FASTA used for repeat localization.")
    locate.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat family catalogue produced by discover.")
    locate.add_argument("--monomers", required=True, type=_path_value, help="Monomer FASTA used as localization query sequences.")
    locate.add_argument("--window-size", type=int, default=100000, help="Window size in bp for assembly repeat density summaries.")
    locate.add_argument("--step-size", type=int, default=10000, help="Step size in bp for sliding-window density summaries.")
    locate.add_argument("--min-identity", type=float, default=0.8, help="Minimum match or alignment identity for future repeat hits.")
    locate.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and future locate outputs.")
    locate.set_defaults(func=run_locate)

    probe = subparsers.add_parser(
        "probe",
        help="Rank candidate FISH probes from repeat families.",
    )
    probe.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat family catalogue produced by discover.")
    probe.add_argument("--monomers", required=True, type=_path_value, help="Monomer FASTA used to derive future probe candidates.")
    probe.add_argument("--copy-number", required=True, type=_path_value, dest="copy_number", help="Read-based copy-number table produced by quantify.")
    probe.add_argument("--locations", type=_path_value, help="Optional repeat density or localization table produced by locate.")
    probe.add_argument("--min-probe-len", type=int, default=80, help="Minimum candidate probe length in bp.")
    probe.add_argument("--max-probe-len", type=int, default=300, help="Maximum candidate probe length in bp.")
    probe.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and future probe outputs.")
    probe.set_defaults(func=run_probe)

    compare = subparsers.add_parser(
        "compare",
        help="Compare read-based and assembly-based repeat abundance.",
    )
    compare.add_argument(
        "--read-copy-number",
        required=True,
        type=_path_value,
        dest="read_copy_number",
        help="Read-based copy-number table produced by quantify.",
    )
    compare.add_argument(
        "--assembly-density",
        required=True,
        type=_path_value,
        dest="assembly_density",
        help="Assembly repeat-density table produced by locate.",
    )
    compare.add_argument("--under-assembly-ratio", type=float, default=2.0, help="Read/assembly ratio threshold for possible under-representation.")
    compare.add_argument("--over-expansion-ratio", type=float, default=0.5, help="Read/assembly ratio threshold for possible over-expansion.")
    compare.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and future compare outputs.")
    compare.set_defaults(func=run_compare)

    visualize = subparsers.add_parser(
        "visualize",
        help="Generate static visual summaries from TandemX outputs.",
    )
    visualize.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat family catalogue produced by discover.")
    visualize.add_argument("--copy-number", required=True, type=_path_value, dest="copy_number", help="Read-based copy-number table produced by quantify.")
    visualize.add_argument("--locations", type=_path_value, help="Optional repeat localization or density table produced by locate.")
    visualize.add_argument("--probes", type=_path_value, help="Optional ranked probe table produced by probe.")
    visualize.add_argument("--comparison", type=_path_value, help="Optional assembly-vs-read comparison table produced by compare.")
    visualize.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and future visualize outputs.")
    visualize.set_defaults(func=run_visualize)

    simulate = subparsers.add_parser(
        "simulate",
        help="Generate simulated toy datasets for TandemX development.",
    )
    simulate_subparsers = simulate.add_subparsers(dest="simulate_command", required=True)
    toy = simulate_subparsers.add_parser(
        "toy",
        help="Generate a reproducible toy tandem-repeat dataset.",
    )
    toy.add_argument("--outdir", required=True, type=_path_value, help="Directory for toy dataset outputs.")
    toy.add_argument("--seed", type=int, default=1, help="Random seed for fully reproducible output.")
    toy.add_argument("--num-reads", type=int, default=40, help="Number of simulated reads to write.")
    toy.add_argument("--read-length", type=int, default=1200, help="Length of each simulated read in bp.")
    toy.add_argument("--background-length", type=int, default=2000, help="Length of random background sequence in bp.")
    toy.add_argument("--monomer-lengths", default="566,350", help="Comma-separated simulated monomer lengths in bp.")
    toy.add_argument("--copies", default="9,7", help="Comma-separated read-truth copy counts for each family.")
    toy.add_argument("--error-rate", type=float, default=0.01, help="Per-base substitution error rate for simulated reads.")
    toy.set_defaults(func=run_simulate_toy)

    return parser


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = list(argv)
    try:
        return args.func(args)
    except (InputFileError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
