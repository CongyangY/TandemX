"""Step-level orchestration for TandemX command-line workflows."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence, TextIO

from tandemx.io.validators import ValidationError, validate_project
from tandemx.reporting import write_output_manifest, write_run_report


PIPELINE_STEPS = ("discover", "quantify", "locate", "compare", "probe", "visualize", "validate")
ASSEMBLY_STEPS = {"locate", "compare", "probe", "visualize"}
SUMMARY_FIELDS = (
    "run_id",
    "input_reads",
    "input_assembly",
    "max_reads",
    "max_read_bases",
    "kmer_backend",
    "step",
    "command",
    "start_time",
    "end_time",
    "runtime_seconds",
    "exit_status",
    "output_dir",
    "output_validated",
    "notes",
)


@dataclass(frozen=True)
class PipelineConfig:
    reads: Path
    assembly: Path | None
    genome_size: int | None
    haploid_depth: float | None
    outdir: Path
    max_reads: int | None
    max_read_bases: int | None
    kmer_backend: str
    steps: tuple[str, ...]
    min_period: int
    max_period: int
    top_periods: int
    threads: int
    resume: bool
    force: bool
    profile: bool


@dataclass(frozen=True)
class StepRecord:
    run_id: str
    input_reads: str
    input_assembly: str
    max_reads: int | None
    max_read_bases: int | None
    kmer_backend: str
    step: str
    command: str
    start_time: str
    end_time: str
    runtime_seconds: float
    exit_status: int
    output_dir: str
    output_validated: bool
    notes: str

    def as_serializable_dict(self) -> dict[str, object]:
        row = asdict(self)
        row["runtime_seconds"] = round(self.runtime_seconds, 6)
        return row


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_steps(value: str) -> tuple[str, ...]:
    requested = [step.strip() for step in value.split(",") if step.strip()]
    unknown = sorted(set(requested) - set(PIPELINE_STEPS))
    if unknown:
        raise argparse.ArgumentTypeError(f"Unknown pipeline step(s): {', '.join(unknown)}")
    if not requested:
        raise argparse.ArgumentTypeError("--steps must contain at least one step")
    if len(requested) != len(set(requested)):
        raise argparse.ArgumentTypeError("--steps must not contain duplicates")
    return tuple(step for step in PIPELINE_STEPS if step in requested)


def add_pipeline_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--assembly", type=Path)
    parser.add_argument("--genome-size", type=int)
    parser.add_argument("--haploid-depth", type=float)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--max-reads", type=int)
    parser.add_argument("--max-read-bases", type=int)
    parser.add_argument("--kmer-backend", choices=("python", "rust"), default="python")
    parser.add_argument(
        "--steps",
        type=parse_steps,
        default=parse_steps("discover,quantify,locate,compare,probe,visualize,validate"),
    )
    parser.add_argument("--min-period", type=int, default=20)
    parser.add_argument("--max-period", type=int, default=2000)
    parser.add_argument("--top-periods", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1, help="Recorded for future parallel execution; currently must be 1.")
    parser.add_argument("--resume", action="store_true", help="Skip existing steps only when their expected outputs validate.")
    parser.add_argument("--force", action="store_true", help="Rerun selected steps even when output directories already exist.")
    parser.add_argument("--profile", action="store_true", help="Write a cProfile file for each executed step.")
    parser.add_argument("--validate", action="store_true", dest="validate_run", help="Append the validate step when it is not listed in --steps.")


def config_from_args(args: argparse.Namespace) -> PipelineConfig:
    steps = args.steps
    if args.validate_run and "validate" not in steps:
        steps = tuple(step for step in PIPELINE_STEPS if step in {*steps, "validate"})
    if args.resume and args.force:
        raise ValueError("--resume and --force are mutually exclusive")
    if args.threads != 1:
        raise ValueError("--threads is reserved for future parallel execution and currently must be 1")
    for name in ("max_reads", "max_read_bases", "genome_size"):
        value = getattr(args, name)
        if value is not None and value <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if args.haploid_depth is not None and args.haploid_depth <= 0:
        raise ValueError("--haploid-depth must be positive")
    return PipelineConfig(
        reads=args.reads,
        assembly=args.assembly,
        genome_size=args.genome_size,
        haploid_depth=args.haploid_depth,
        outdir=args.outdir,
        max_reads=args.max_reads,
        max_read_bases=args.max_read_bases,
        kmer_backend=args.kmer_backend,
        steps=steps,
        min_period=args.min_period,
        max_period=args.max_period,
        top_periods=args.top_periods,
        threads=args.threads,
        resume=args.resume,
        force=args.force,
        profile=args.profile,
    )


def build_step_command(config: PipelineConfig, step: str) -> list[str]:
    cli = [sys.executable, "-m", "tandemx.cli"]
    discover_dir = config.outdir / "discover"
    quantify_dir = config.outdir / "quantify"
    locate_dir = config.outdir / "locate"
    compare_dir = config.outdir / "compare"
    probe_dir = config.outdir / "probe"
    visualize_dir = config.outdir / "visualize"
    catalog = discover_dir / "monomers.fa"
    copy_number = quantify_dir / "copy_number.tsv"

    if step == "discover":
        command = [
            *cli,
            "discover",
            "--reads",
            str(config.reads),
            "--outdir",
            str(discover_dir),
            "--kmer-backend",
            config.kmer_backend,
            "--min-period",
            str(config.min_period),
            "--max-period",
            str(config.max_period),
            "--top-periods",
            str(config.top_periods),
        ]
        if config.max_reads is not None:
            command.extend(["--max-reads", str(config.max_reads)])
        if config.max_read_bases is not None:
            command.extend(["--max-read-bases", str(config.max_read_bases)])
        return command
    if step == "quantify":
        if config.genome_size is None:
            raise ValueError("--genome-size is required when quantify is selected")
        command = [
            *cli,
            "quantify",
            "--reads",
            str(config.reads),
            "--catalog",
            str(catalog),
            "--genome-size",
            str(config.genome_size),
            "--outdir",
            str(quantify_dir),
            "--kmer-backend",
            config.kmer_backend,
        ]
        if config.haploid_depth is not None:
            command.extend(["--haploid-depth", str(config.haploid_depth)])
        if config.max_reads is not None:
            command.extend(["--max-reads", str(config.max_reads)])
        if config.max_read_bases is not None:
            command.extend(["--max-read-bases", str(config.max_read_bases)])
        return command
    if step == "locate":
        assert config.assembly is not None
        return [
            *cli,
            "locate",
            "--assembly",
            str(config.assembly),
            "--catalog",
            str(catalog),
            "--copy-number",
            str(copy_number),
            "--outdir",
            str(locate_dir),
        ]
    if step == "compare":
        return [
            *cli,
            "compare",
            "--copy-number",
            str(copy_number),
            "--arrays",
            str(locate_dir / "arrays.bed"),
            "--outdir",
            str(compare_dir),
        ]
    if step == "probe":
        assert config.assembly is not None
        return [
            *cli,
            "probe",
            "--catalog",
            str(catalog),
            "--assembly",
            str(config.assembly),
            "--copy-number",
            str(copy_number),
            "--arrays",
            str(locate_dir / "arrays.bed"),
            "--outdir",
            str(probe_dir),
        ]
    if step == "visualize":
        command = [
            *cli,
            "visualize",
            "--catalog",
            str(catalog),
            "--copy-number",
            str(copy_number),
            "--outdir",
            str(visualize_dir),
        ]
        if config.assembly is not None:
            command.extend(
                [
                    "--comparison",
                    str(compare_dir / "assembly_vs_read_cn.tsv"),
                    "--probes",
                    str(probe_dir / "probes.rank.tsv"),
                    "--fish",
                    str(probe_dir / "in_silico_fish.tsv"),
                ]
            )
        return command
    if step == "validate":
        return [*cli, "validate", "--project", str(config.outdir)]
    raise ValueError(f"Unsupported pipeline step: {step}")


def expected_outputs(config: PipelineConfig, step: str) -> tuple[Path, ...]:
    output_dir = config.outdir / step
    outputs = {
        "discover": (output_dir / "candidate_reads.tsv", output_dir / "monomers.fa", output_dir / "families.tsv"),
        "quantify": (output_dir / "copy_number.tsv",),
        "locate": (output_dir / "repeat_density.bedgraph", output_dir / "arrays.bed", output_dir / "assembly_vs_read_cn.tsv"),
        "compare": (output_dir / "assembly_vs_read_cn.tsv",),
        "probe": (output_dir / "probes.fa", output_dir / "probes.rank.tsv", output_dir / "in_silico_fish.tsv"),
        "visualize": (output_dir / "catalogue_summary.svg", output_dir / "catalogue_summary.pdf"),
        "validate": (),
    }
    return outputs[step]


def step_outputs_validate(config: PipelineConfig, step: str) -> bool:
    outputs = expected_outputs(config, step)
    if not outputs or not all(path.is_file() for path in outputs):
        return False
    if step == "visualize":
        return all(path.stat().st_size > 0 for path in outputs)
    try:
        validate_project(config.outdir / step)
    except ValidationError:
        return False
    return True


def missing_step_inputs(config: PipelineConfig, step: str) -> tuple[Path, ...]:
    if step == "compare":
        required = (config.outdir / "quantify" / "copy_number.tsv", config.outdir / "locate" / "arrays.bed")
        return tuple(path for path in required if not path.is_file())
    return ()


def write_summaries(config: PipelineConfig, records: Sequence[StepRecord]) -> None:
    config.outdir.mkdir(parents=True, exist_ok=True)
    rows = [record.as_serializable_dict() for record in records]
    with (config.outdir / "pipeline_summary.tsv").open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (config.outdir / "pipeline_summary.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def finalize_run_outputs(config: PipelineConfig, records: Sequence[StepRecord]) -> None:
    write_summaries(config, records)
    write_run_report(config, records)
    write_output_manifest(config, records)


def make_record(
    config: PipelineConfig,
    run_id: str,
    step: str,
    *,
    command: Sequence[str] = (),
    start_time: str,
    end_time: str,
    runtime_seconds: float,
    exit_status: int,
    output_validated: bool,
    notes: str,
) -> StepRecord:
    return StepRecord(
        run_id=run_id,
        input_reads=str(config.reads),
        input_assembly=str(config.assembly) if config.assembly else "",
        max_reads=config.max_reads,
        max_read_bases=config.max_read_bases,
        kmer_backend=config.kmer_backend,
        step=step,
        command=shlex.join(command),
        start_time=start_time,
        end_time=end_time,
        runtime_seconds=runtime_seconds,
        exit_status=exit_status,
        output_dir=str(config.outdir / step),
        output_validated=output_validated,
        notes=notes,
    )


def run_command_with_live_logs(
    command: Sequence[str],
    stdout_path: Path,
    stderr_path: Path,
) -> int:
    """Run a command while teeing child output to log files and this terminal."""
    with stdout_path.open("wt", encoding="utf-8") as stdout_log, stderr_path.open(
        "wt", encoding="utf-8"
    ) as stderr_log:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        def forward(pipe: TextIO, log_handle: TextIO, mirror: TextIO) -> None:
            for line in pipe:
                log_handle.write(line)
                log_handle.flush()
                mirror.write(line)
                mirror.flush()

        threads = [
            threading.Thread(target=forward, args=(process.stdout, stdout_log, sys.stdout)),
            threading.Thread(target=forward, args=(process.stderr, stderr_log, sys.stderr)),
        ]
        for thread in threads:
            thread.start()
        returncode = process.wait()
        for thread in threads:
            thread.join()
        return returncode


def run_pipeline(config: PipelineConfig) -> tuple[list[StepRecord], int]:
    if not config.reads.is_file():
        raise ValueError(f"Input reads file does not exist: {config.reads}")
    if config.assembly is not None and not config.assembly.is_file():
        raise ValueError(f"Input assembly file does not exist: {config.assembly}")
    config.outdir.mkdir(parents=True, exist_ok=True)
    logs_dir = config.outdir / "logs"
    profiles_dir = config.outdir / "profiles"
    logs_dir.mkdir(parents=True, exist_ok=True)
    if config.profile:
        profiles_dir.mkdir(parents=True, exist_ok=True)
    run_id = f"tandemx-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{time.time_ns()}"
    records: list[StepRecord] = []
    pipeline_log = config.outdir / "pipeline.log"
    pipeline_log.touch(exist_ok=True)

    for step in config.steps:
        start_time = utc_now()
        if config.assembly is None and step in ASSEMBLY_STEPS:
            record = make_record(
                config,
                run_id,
                step,
                start_time=start_time,
                end_time=utc_now(),
                runtime_seconds=0.0,
                exit_status=0,
                output_validated=False,
                notes="skipped_missing_assembly",
            )
            records.append(record)
            write_summaries(config, records)
            continue

        missing_inputs = missing_step_inputs(config, step)
        if missing_inputs:
            record = make_record(
                config,
                run_id,
                step,
                start_time=start_time,
                end_time=utc_now(),
                runtime_seconds=0.0,
                exit_status=2,
                output_validated=False,
                notes="missing_input:" + ",".join(str(path) for path in missing_inputs),
            )
            records.append(record)
            finalize_run_outputs(config, records)
            return records, 2

        output_dir = config.outdir / step
        if step != "validate" and output_dir.exists() and any(output_dir.iterdir()):
            if config.resume and step_outputs_validate(config, step):
                record = make_record(
                    config,
                    run_id,
                    step,
                    start_time=start_time,
                    end_time=utc_now(),
                    runtime_seconds=0.0,
                    exit_status=0,
                    output_validated=True,
                    notes="skipped_validated_resume",
                )
                records.append(record)
                write_summaries(config, records)
                continue
            if not config.force and not config.resume:
                record = make_record(
                    config,
                    run_id,
                    step,
                    start_time=start_time,
                    end_time=utc_now(),
                    runtime_seconds=0.0,
                    exit_status=2,
                    output_validated=False,
                    notes="output_exists_use_force_or_resume",
                )
                records.append(record)
                finalize_run_outputs(config, records)
                return records, 2

        try:
            command = build_step_command(config, step)
        except ValueError as exc:
            record = make_record(
                config,
                run_id,
                step,
                start_time=start_time,
                end_time=utc_now(),
                runtime_seconds=0.0,
                exit_status=2,
                output_validated=False,
                notes=str(exc),
            )
            records.append(record)
            finalize_run_outputs(config, records)
            return records, 2

        output_dir.mkdir(parents=True, exist_ok=True)

        actual_command = command
        if config.profile:
            actual_command = [
                sys.executable,
                "-m",
                "cProfile",
                "-o",
                str(profiles_dir / f"{step}.prof"),
                *command[1:],
            ]
        print(
            f"tandemx run: starting step {len(records) + 1}/{len(config.steps)}: {step}",
            flush=True,
        )
        started = time.perf_counter()
        returncode = run_command_with_live_logs(
            actual_command,
            logs_dir / f"{step}.stdout.log",
            logs_dir / f"{step}.stderr.log",
        )
        runtime = time.perf_counter() - started
        end_time = utc_now()
        print(
            f"tandemx run: finished step {step} exit_status={returncode} runtime_seconds={runtime:.3f}",
            flush=True,
        )
        validated = returncode == 0 and (
            step == "validate" or step_outputs_validate(config, step)
        )
        notes = f"threads_recorded={config.threads}"
        if config.profile:
            notes += ";cprofile_written"
        record = make_record(
            config,
            run_id,
            step,
            command=actual_command,
            start_time=start_time,
            end_time=end_time,
            runtime_seconds=runtime,
            exit_status=returncode,
            output_validated=validated,
            notes=notes,
        )
        records.append(record)
        with pipeline_log.open("at", encoding="utf-8") as handle:
            handle.write(
                f"{end_time} step={step} exit_status={returncode} "
                f"runtime_seconds={runtime:.6f} validated={str(validated).lower()}\n"
            )
        write_summaries(config, records)
        if returncode != 0:
            finalize_run_outputs(config, records)
            return records, returncode
    finalize_run_outputs(config, records)
    return records, 0


def run_pipeline_cli(args: argparse.Namespace) -> int:
    config = config_from_args(args)
    records, exit_status = run_pipeline(config)
    total_runtime = sum(record.runtime_seconds for record in records)
    print(f"tandemx run: recorded {len(records)} steps in {total_runtime:.3f} seconds")
    print(f"output directory: {config.outdir}")
    for label, path in (
        ("families.tsv", config.outdir / "discover" / "families.tsv"),
        ("monomers.fa", config.outdir / "discover" / "monomers.fa"),
        ("copy_number.tsv", config.outdir / "quantify" / "copy_number.tsv"),
        ("assembly_vs_read_cn.tsv", config.outdir / "compare" / "assembly_vs_read_cn.tsv"),
    ):
        if path.is_file():
            print(f"{label}: {path}")
    print(f"output_manifest.tsv: {config.outdir / 'output_manifest.tsv'}")
    print(f"run_report.md: {config.outdir / 'run_report.md'}")
    print(f"pipeline_summary.tsv: {config.outdir / 'pipeline_summary.tsv'}")
    print(f"validation: tandemx validate --project {shlex.quote(str(config.outdir))}")
    return exit_status
