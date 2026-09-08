"""Frozen multi-seed TideCluster factorial-assembly execution."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

from benchmarks.challenge.schema import digest_file, write_table
from benchmarks.scripts.profile_stage_resources import profile
from benchmarks.scripts.run_tidecluster_docker_reference import docker_metadata, parse_gnu_time
from benchmarks.tidecluster.factorial_evaluate import evaluate_factorial


SUMMARY_FIELDS = [
    "seed", "setting", "status", "tidehunter_wall_seconds", "tidehunter_maximum_rss_kb",
    "clustering_wall_seconds", "clustering_maximum_rss_kb", "truth_array_count",
    "predicted_array_count", "matched_array_count", "array_recall", "array_precision",
    "base_union_recall", "base_union_precision", "matched_boundary_mae_bp",
    "matched_period_mae_bp", "truth_family_count", "recovered_family_count",
    "cyclic_monomer_recall", "distinct_consensus_count", "homologous_consensus_fraction",
    "operational_family_count", "warning",
]


def command_output(command: list[str]) -> str:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return (completed.stdout or completed.stderr).strip()


def validate_datasets(config: dict[str, object]) -> list[dict[str, object]]:
    datasets = config.get("datasets")
    if not isinstance(datasets, list) or not datasets:
        raise ValueError("Config requires non-empty datasets")
    seen: set[int] = set()
    validated: list[dict[str, object]] = []
    for row in datasets:
        seed = int(row["seed"])
        if seed in seen:
            raise ValueError(f"Duplicate dataset seed: {seed}")
        seen.add(seed)
        genome_dir = Path(row["genome_dir"]).resolve()
        manifest = json.loads((genome_dir / "manifest.json").read_text(encoding="utf-8"))
        expected = {
            "genome.fa": row["genome_sha256"],
            "catalogue.fa": row["catalogue_sha256"],
            "truth_copy_number.tsv": row["truth_sha256"],
        }
        if (
            manifest.get("generator") != "streamed_factorial_genome_v1"
            or int(manifest.get("seed", -1)) != seed
            or int(manifest.get("genome_bp", -1)) != 10_000_000
        ):
            raise ValueError(f"Unexpected factorial manifest for seed {seed}")
        for name, expected_sha256 in expected.items():
            path = genome_dir / name
            if manifest["files"].get(name) != expected_sha256 or digest_file(path) != expected_sha256:
                raise ValueError(f"Frozen dataset file differs for seed {seed}: {name}")
        validated.append({"seed": seed, "genome_dir": genome_dir, "manifest": manifest})
    return validated


def build_stage_specs(
    config: dict[str, object],
    datasets: list[dict[str, object]],
    outdir: Path,
    docker: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    settings = config.get("settings")
    if not isinstance(settings, dict) or not settings:
        raise ValueError("Config requires non-empty settings")
    cpus = int(config["cpus"])
    minimum_length = int(config["minimum_length"])
    if not 1 <= cpus <= 64 or minimum_length <= 0:
        raise ValueError("Invalid CPU or minimum-length setting")
    stages: list[dict[str, object]] = []
    runs: list[dict[str, object]] = []
    for dataset in datasets:
        seed = int(dataset["seed"])
        assembly = Path(dataset["genome_dir"]) / "genome.fa"
        for setting_name, setting in settings.items():
            if not setting_name.replace("_", "").isalnum():
                raise ValueError(f"Unsafe setting name: {setting_name}")
            tidehunter_arguments = setting.get("tidehunter_arguments")
            if not isinstance(tidehunter_arguments, str) or not tidehunter_arguments:
                raise ValueError(f"Missing TideHunter arguments: {setting_name}")
            run_dir = outdir / f"seed{seed}" / setting_name
            run_dir.mkdir(parents=True)
            volume_input = f"{assembly}:/input/genome.fa:ro"
            volume_output = f"{run_dir}:/output"
            base = [
                docker, "run", "--rm", "--platform", "linux/amd64",
                "-v", volume_input, "-v", volume_output, "-w", "/output",
                "--entrypoint", "/opt/conda/envs/tidecluster/bin/time",
                config["image"], "-v",
            ]
            prefix = f"s{seed}_{setting_name}"
            stages.extend(
                [
                    {
                        "name": f"{prefix}_tidehunter",
                        "command": [
                            *base, "-o", "/output/tidehunter.gnu_time.txt",
                            "TideCluster.py", "tidehunter", "-c", str(cpus),
                            "-T", tidehunter_arguments, "-pr", "/output/tc",
                            "-f", "/input/genome.fa",
                        ],
                        "cwd": str(run_dir),
                        "scratch_dir": str(run_dir),
                    },
                    {
                        "name": f"{prefix}_clustering",
                        "command": [
                            *base, "-o", "/output/clustering.gnu_time.txt",
                            "TideCluster.py", "clustering", "-c", str(cpus),
                            "-m", str(minimum_length), "-pr", "/output/tc",
                            "-f", "/input/genome.fa",
                        ],
                        "cwd": str(run_dir),
                        "scratch_dir": str(run_dir),
                    },
                ]
            )
            runs.append(
                {
                    "seed": seed,
                    "setting": setting_name,
                    "run_dir": run_dir,
                    "genome_dir": Path(dataset["genome_dir"]),
                }
            )
    return stages, runs


def run_validation(
    config_path: Path,
    outdir: Path,
    docker: str = "docker",
    interval: float = 0.2,
) -> list[dict[str, object]]:
    if outdir.exists():
        raise FileExistsError(f"Choose a new output directory: {outdir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    datasets = validate_datasets(config)
    metadata = docker_metadata(docker, config["image"])
    if metadata["image_id"] != config["image_id"]:
        raise ValueError("Docker image ID differs from frozen config")
    project_root = Path(__file__).resolve().parents[2]
    project_commit = command_output(["git", "-C", str(project_root), "rev-parse", "HEAD"])
    project_status = command_output(
        ["git", "-C", str(project_root), "status", "--short", "--untracked-files=no"]
    )
    if project_status:
        raise ValueError("Tracked worktree must be clean before validation execution")

    outdir.mkdir(parents=True)
    stages, runs = build_stage_specs(config, datasets, outdir, docker)
    stage_manifest = outdir / "stage_manifest.json"
    stage_manifest.write_text(json.dumps({"stages": stages}, indent=2) + "\n", encoding="utf-8")
    snapshot = outdir / "source_snapshot"
    snapshot.mkdir()
    source_paths = (
        config_path,
        Path(__file__),
        project_root / "benchmarks/tidecluster/factorial_evaluate.py",
        project_root / "benchmarks/tidecluster/normalize.py",
        project_root / "benchmarks/scripts/profile_stage_resources.py",
        project_root / "benchmarks/scripts/run_tidecluster_docker_reference.py",
    )
    for source in source_paths:
        shutil.copyfile(source, snapshot / source.name)
    environment = {
        "schema_version": 1,
        "complete": False,
        "experiment_id": config["experiment_id"],
        "project_commit": project_commit,
        "project_status": project_status,
        "config_sha256": digest_file(config_path),
        "tool": metadata,
        "datasets": [
            {
                "seed": row["seed"],
                "genome_dir": str(row["genome_dir"]),
                "manifest_sha256": digest_file(Path(row["genome_dir"]) / "manifest.json"),
            }
            for row in datasets
        ],
        "boundary": config["boundary"],
    }
    environment_path = outdir / "environment.json"
    environment_path.write_text(json.dumps(environment, indent=2) + "\n", encoding="utf-8")
    profile_exit = profile(stage_manifest, outdir / "profile", interval)
    profile_receipt = json.loads((outdir / "profile/receipt.json").read_text(encoding="utf-8"))
    if profile_exit != 0 or profile_receipt.get("complete") is not True:
        (outdir / "run_receipt.json").write_text(
            json.dumps(
                {
                    "complete": False,
                    "fate": "external_process_failure",
                    "profile": profile_receipt,
                    "warning": "failed_runs_are_missing_measurements_not_zero_accuracy",
                },
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
        raise RuntimeError("TideCluster validation stage failed; retained outputs are authoritative")

    rows: list[dict[str, object]] = []
    try:
        for run in runs:
            run_dir = Path(run["run_dir"])
            metrics = evaluate_factorial(
                run_dir / "tc_tidehunter.gff3",
                run_dir / "tc_clustering.gff3_1.gff3",
                run_dir / "tc_clustering.gff3",
                run_dir / "tc_consensus/consensus_sequences_all.fasta",
                Path(run["genome_dir"]),
                run_dir / "evaluation",
            )
            resources = {
                name: parse_gnu_time(run_dir / f"{name}.gnu_time.txt")
                for name in ("tidehunter", "clustering")
            }
            rows.append(
                {
                    "seed": run["seed"],
                    "setting": run["setting"],
                    "status": "ok",
                    "tidehunter_wall_seconds": resources["tidehunter"]["wall_seconds"],
                    "tidehunter_maximum_rss_kb": resources["tidehunter"]["maximum_rss_kb"],
                    "clustering_wall_seconds": resources["clustering"]["wall_seconds"],
                    "clustering_maximum_rss_kb": resources["clustering"]["maximum_rss_kb"],
                    **{name: metrics.get(name) for name in SUMMARY_FIELDS[7:]},
                }
            )
    except Exception as error:
        (outdir / "run_receipt.json").write_text(
            json.dumps(
                {
                    "complete": False,
                    "fate": "normalization_or_evaluation_failure",
                    "error": str(error),
                    "completed_summary_rows": rows,
                },
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
        raise
    write_table(outdir / "summary.tsv", rows, SUMMARY_FIELDS)
    for row in datasets:
        genome_dir = Path(row["genome_dir"])
        for name, expected in (
            ("genome.fa", row["manifest"]["files"]["genome.fa"]),
            ("catalogue.fa", row["manifest"]["files"]["catalogue.fa"]),
            ("truth_copy_number.tsv", row["manifest"]["files"]["truth_copy_number.tsv"]),
        ):
            if digest_file(genome_dir / name) != expected:
                raise ValueError("A frozen validation input changed during execution")
    environment["complete"] = True
    environment_path.write_text(json.dumps(environment, indent=2) + "\n", encoding="utf-8")
    receipt = {
        "schema_version": 1,
        "complete": True,
        "run_count": len(rows),
        "successful_run_count": sum(row["status"] == "ok" for row in rows),
        "summary_sha256": digest_file(outdir / "summary.tsv"),
        "warning": (
            "same_process_simulated_validation_genomes_not_independent_plants;"
            "technical_repetitions_not_biological_replicates"
        ),
    }
    (outdir / "run_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return rows
