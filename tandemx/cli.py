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
from tandemx.discover.mvp import DiscoverConfig, discover_toy_repeats
from tandemx.io.sequences import SequenceFormatError
from tandemx.io.validators import ValidationError, validate_project
from tandemx.locate.mvp import LocateConfig, locate_toy_arrays
from tandemx.probe.mvp import ProbeConfig, rank_toy_probes
from tandemx.quantify.mvp import QuantifyConfig, quantify_toy_copy_number
from tandemx.simulate.toy import ToySimulationConfig, generate_toy_dataset, parse_int_list
from tandemx.visualize.mvp import VisualizeConfig, render_static_plots


COMMANDS = ("discover", "quantify", "locate", "probe", "compare", "visualize", "simulate", "validate")


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

    handler = logging.FileHandler(outdir / "run.log", mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    return logger


def _prepare_deferred_run(
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
    logger.info("Core algorithm is deferred in this MVP")

    print(
        f"tandemx {command}: deferred in this MVP. "
        f"Wrote run_config.yaml and run.log to {args.outdir}"
    )
    return 0


def run_discover(args: argparse.Namespace) -> int:
    _require_existing_files([args.reads])
    args.outdir.mkdir(parents=True, exist_ok=True)
    _write_run_config(args.outdir, "discover", args, status="discover_running")
    logger = _configure_log(args.outdir, "discover")
    logger.info("command=tandemx discover")
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=discover_running")
    logger.info(
        "algorithm_mode=spacing_prefilter kmer_backend=%s min_period=%s max_period=%s",
        args.kmer_backend,
        args.min_monomer_len,
        args.max_monomer_len,
    )
    config = DiscoverConfig(
        reads=args.reads,
        outdir=args.outdir,
        min_monomer_len=args.min_monomer_len,
        max_monomer_len=args.max_monomer_len,
        min_support_reads=args.min_support_reads,
        min_repeat_span=args.min_repeat_span,
        min_read_length=args.min_read_length,
        kmer_size=args.kmer_size,
        top_periods=args.top_periods,
        min_seed_occurrences=args.min_seed_occurrences,
        min_spacing_support=args.min_spacing_support,
        max_pairs_per_kmer=args.max_pairs_per_kmer,
        max_reads=args.max_reads,
        max_read_bases=args.max_read_bases,
        sample_rate=args.sample_rate,
        seed=args.seed,
        progress_every=args.progress_every,
        chunk_size=args.chunk_size,
        kmer_backend=args.kmer_backend,
    )
    try:
        candidates, families = discover_toy_repeats(config, logger=logger)
    except KeyboardInterrupt:
        logger.warning("status=discover_interrupted partial_candidate_output=%s", args.outdir / "candidate_reads.tsv")
        raise
    _write_run_config(args.outdir, "discover", args, status="discover_mvp_completed")
    logger.info("status=discover_mvp_completed")
    logger.info("candidate_count=%s", len(candidates))
    logger.info("family_count=%s", len(families))
    logger.info("Discover MVP supports toy-scale FASTA/FASTQ input, including gzip-compressed files")
    print(
        f"tandemx discover: wrote {len(candidates)} candidates and "
        f"{len(families)} families to {args.outdir}"
    )
    return 0


def run_quantify(args: argparse.Namespace) -> int:
    monomers_path = args.monomers if args.monomers is not None else args.catalogue
    _require_existing_files([args.reads, monomers_path])
    args.outdir.mkdir(parents=True, exist_ok=True)
    estimates = quantify_toy_copy_number(
        QuantifyConfig(
            reads=args.reads,
            monomers=monomers_path,
            genome_size=args.genome_size,
            outdir=args.outdir,
            k=args.k,
            haploid_depth=args.haploid_depth,
        )
    )
    _write_run_config(args.outdir, "quantify", args, status="quantify_mvp_completed")
    logger = _configure_log(args.outdir, "quantify")
    logger.info("command=tandemx quantify")
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=quantify_mvp_completed")
    logger.info("family_count=%s", len(estimates))
    logger.info("Quantify MVP supports toy-scale FASTA/FASTQ read input, including gzip-compressed files")
    print(f"tandemx quantify: wrote copy-number estimates for {len(estimates)} families to {args.outdir}")
    return 0


def run_locate(args: argparse.Namespace) -> int:
    monomers_path = args.monomers if args.monomers is not None else args.catalogue
    required = [args.assembly, monomers_path]
    if args.copy_number is not None:
        required.append(args.copy_number)
    _require_existing_files(required)
    args.outdir.mkdir(parents=True, exist_ok=True)
    density, arrays, comparisons = locate_toy_arrays(
        LocateConfig(
            assembly=args.assembly,
            monomers=monomers_path,
            copy_number=args.copy_number,
            outdir=args.outdir,
            window_size=args.window_size,
            step_size=args.step_size,
            k=args.k,
        )
    )
    _write_run_config(args.outdir, "locate", args, status="locate_mvp_completed")
    logger = _configure_log(args.outdir, "locate")
    logger.info("command=tandemx locate")
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=locate_mvp_completed")
    logger.info("density_window_count=%s", len(density))
    logger.info("array_count=%s", len(arrays))
    logger.info("comparison_count=%s", len(comparisons))
    logger.info("Locate MVP supports toy-scale FASTA assembly input, including gzip-compressed files")
    print(f"tandemx locate: wrote {len(arrays)} arrays and {len(comparisons)} comparisons to {args.outdir}")
    return 0


def run_probe(args: argparse.Namespace) -> int:
    monomers_path = args.monomers if args.monomers is not None else args.catalogue
    arrays_path = args.arrays if args.arrays is not None else args.locations
    if args.assembly is None:
        raise ValueError("--assembly is required for probe MVP")
    if arrays_path is None:
        raise ValueError("--arrays is required for probe MVP")
    _require_existing_files([monomers_path, args.assembly, args.copy_number, arrays_path])
    args.outdir.mkdir(parents=True, exist_ok=True)
    candidates, signals = rank_toy_probes(
        ProbeConfig(
            monomers=monomers_path,
            assembly=args.assembly,
            copy_number=args.copy_number,
            arrays=arrays_path,
            outdir=args.outdir,
            min_len=args.min_len,
            max_len=args.max_len,
        )
    )
    _write_run_config(args.outdir, "probe", args, status="probe_mvp_completed")
    logger = _configure_log(args.outdir, "probe")
    logger.info("command=tandemx probe")
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=probe_mvp_completed")
    logger.info("probe_count=%s", len(candidates))
    logger.info("signal_count=%s", len(signals))
    logger.info("Probe MVP is a toy-scale prioritization heuristic")
    print(f"tandemx probe: wrote {len(candidates)} probes and {len(signals)} predicted signals to {args.outdir}")
    return 0


def run_compare(args: argparse.Namespace) -> int:
    return _prepare_deferred_run(
        "compare",
        args,
        [args.read_copy_number, args.assembly_density],
    )


def run_visualize(args: argparse.Namespace) -> int:
    required = [args.catalogue, args.copy_number]
    for optional in (args.locations, args.probes, args.comparison, args.fish):
        if optional is not None:
            required.append(optional)
    _require_existing_files(required)
    args.outdir.mkdir(parents=True, exist_ok=True)
    outputs = render_static_plots(
        VisualizeConfig(
            copy_number=args.copy_number,
            comparison=args.comparison,
            probes=args.probes,
            fish=args.fish,
            outdir=args.outdir,
        )
    )
    _write_run_config(args.outdir, "visualize", args, status="visualize_mvp_completed")
    logger = _configure_log(args.outdir, "visualize")
    logger.info("command=tandemx visualize")
    logger.info("timestamp_utc=%s", datetime.now(timezone.utc).isoformat())
    logger.info("output_directory=%s", args.outdir)
    logger.info("status=visualize_mvp_completed")
    logger.info("plot_count=%s", len(outputs))
    print(f"tandemx visualize: wrote {len(outputs)} plot files to {args.outdir}")
    return 0


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


def run_validate(args: argparse.Namespace) -> int:
    if not args.project.exists():
        raise InputFileError(f"Project directory does not exist: {args.project}")
    if not args.project.is_dir():
        raise InputFileError(f"Project path is not a directory: {args.project}")
    results = validate_project(args.project)
    total_records = sum(result.record_count for result in results)
    print(
        f"tandemx validate: validated {len(results)} output files "
        f"with {total_records} records under {args.project}"
    )
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
    discover.add_argument("--reads", required=True, type=_path_value, help="Input toy-scale HiFi-like reads in FASTA/FASTQ format, optionally gzip-compressed.")
    discover.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, candidate_reads.tsv, monomers.fa, and families.tsv.")
    discover.add_argument("--min-period", "--min-monomer-len", dest="min_monomer_len", type=int, default=20, help="Minimum candidate repeat period in bp.")
    discover.add_argument("--max-period", "--max-monomer-len", dest="max_monomer_len", type=int, default=2000, help="Maximum candidate repeat period in bp.")
    discover.add_argument("--min-support-reads", type=int, default=5, help="Minimum number of reads supporting a candidate family.")
    discover.add_argument("--min-repeat-span", type=int, default=100, help="Minimum repeat-supporting span in a read, in bp.")
    discover.add_argument("--min-read-length", type=int, default=1, help="Skip reads shorter than this length in bp; set explicitly for real-read pilots.")
    discover.add_argument("--kmer-size", type=int, default=11, help="Canonical seed k-mer size for spacing prefiltering.")
    discover.add_argument("--top-periods", type=int, default=5, help="Maximum spacing peaks refined per read.")
    discover.add_argument("--min-seed-occurrences", type=int, default=2, help="Minimum within-read occurrences required for a seed k-mer.")
    discover.add_argument("--min-spacing-support", type=int, default=2, help="Minimum repeated-seed support required for a spacing peak.")
    discover.add_argument("--max-pairs-per-kmer", type=int, default=100, help="Maximum adjacent position pairs retained per seed k-mer.")
    discover.add_argument("--max-reads", type=int, help="Maximum sampled reads to process for pilot runs.")
    discover.add_argument("--max-read-bases", type=int, help="Maximum cumulative sampled read bases to process for pilot runs.")
    discover.add_argument("--sample-rate", type=float, default=1.0, help="Fraction of input reads sampled reproducibly, in (0, 1].")
    discover.add_argument("--seed", type=int, default=1, help="Random seed for reproducible read sampling.")
    discover.add_argument("--progress-every", type=int, default=1000, help="Log progress after this many processed reads.")
    discover.add_argument("--chunk-size", type=int, default=1000, help="Logical read chunk size reserved for future parallel/checkpoint execution.")
    discover.add_argument("--kmer-backend", choices=("python", "rust"), default="python", help="Read-local seed backend. Rust requires the compiled extension; Python remains the fallback/default.")
    discover.set_defaults(func=run_discover)

    quantify = subparsers.add_parser(
        "quantify",
        help="Estimate read-based repeat copy number from diagnostic k-mers.",
    )
    quantify.add_argument("--reads", required=True, type=_path_value, help="Input toy-scale FASTA/FASTQ reads used to estimate diagnostic k-mer depth, optionally gzip-compressed.")
    quantify.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat catalog FASTA generated by `tandemx discover`. The --catalog spelling is accepted.")
    quantify.add_argument("--monomers", type=_path_value, help="Optional explicit path to the discovered monomer catalog FASTA. If omitted, --catalogue/--catalog is used.")
    quantify.add_argument("--genome-size", required=True, type=int, help="Estimated haploid or target genome size in bp.")
    quantify.add_argument("--k", type=int, default=21, help="Diagnostic k-mer size.")
    quantify.add_argument("--haploid-depth", type=float, help="Optional haploid sequencing depth. If omitted, depth is estimated as total read bases divided by genome size.")
    quantify.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, and copy_number.tsv.")
    quantify.set_defaults(func=run_quantify)

    locate = subparsers.add_parser(
        "locate",
        help="Localize repeat family evidence on an assembly.",
    )
    locate.add_argument("--assembly", required=True, type=_path_value, help="Toy-scale genome assembly FASTA used for repeat localization, optionally gzip-compressed.")
    locate.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat catalog FASTA generated by `tandemx discover`. The --catalog spelling is accepted.")
    locate.add_argument("--monomers", type=_path_value, help="Optional explicit path to the discovered monomer catalog FASTA. If omitted, --catalogue/--catalog is used.")
    locate.add_argument("--copy-number", type=_path_value, dest="copy_number", help="Optional copy_number.tsv produced by quantify for assembly-vs-read comparison.")
    locate.add_argument("--window-size", type=int, default=100000, help="Window size in bp for assembly repeat density summaries.")
    locate.add_argument("--step-size", type=int, default=10000, help="Step size in bp for sliding-window density summaries.")
    locate.add_argument("--k", type=int, default=21, help="K-mer size used to scan assembly for monomer evidence.")
    locate.add_argument("--min-identity", type=float, default=0.8, help="Minimum match or alignment identity for future repeat hits.")
    locate.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, repeat_density.bedgraph, arrays.bed, and assembly_vs_read_cn.tsv.")
    locate.set_defaults(func=run_locate)

    probe = subparsers.add_parser(
        "probe",
        help="Rank candidate FISH probes from repeat families.",
    )
    probe.add_argument("--catalogue", "--catalog", required=True, type=_path_value, help="Repeat catalog FASTA generated by `tandemx discover`. The --catalog spelling is accepted.")
    probe.add_argument("--monomers", type=_path_value, help="Optional explicit path to the discovered monomer catalog FASTA. If omitted, --catalogue/--catalog is used.")
    probe.add_argument("--assembly", type=_path_value, help="Toy-scale assembly FASTA used to estimate probe target and off-target regions.")
    probe.add_argument("--copy-number", required=True, type=_path_value, dest="copy_number", help="Read-based copy-number table produced by quantify.")
    probe.add_argument("--arrays", type=_path_value, help="arrays.bed produced by locate.")
    probe.add_argument("--locations", type=_path_value, help="Backward-compatible alias for --arrays.")
    probe.add_argument("--min-len", "--min-probe-len", type=int, default=80, help="Minimum candidate probe length in bp.")
    probe.add_argument("--max-len", "--max-probe-len", type=int, default=300, help="Maximum candidate probe length in bp.")
    probe.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, probes.fa, probes.rank.tsv, and in_silico_fish.tsv.")
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
    visualize.add_argument("--fish", type=_path_value, help="Optional in_silico_fish.tsv produced by probe.")
    visualize.add_argument("--outdir", required=True, type=_path_value, help="Directory for run_config.yaml, run.log, SVG, and PDF static plots.")
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

    validate = subparsers.add_parser(
        "validate",
        help="Validate TandemX MVP output files under a project directory.",
    )
    validate.add_argument("--project", required=True, type=_path_value, help="Project or output directory to scan for TandemX output files.")
    validate.set_defaults(func=run_validate)

    return parser


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = list(argv)
    try:
        return args.func(args)
    except KeyboardInterrupt:
        parser.exit(130, "error: interrupted by user; partial outputs may be available\n")
    except (InputFileError, SequenceFormatError, ValidationError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
