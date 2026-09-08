"""Static and Git-provenance audit of an unmodified unitFinder checkout."""
from __future__ import annotations

import ast
import json
from pathlib import Path
import re
import subprocess

from benchmarks.challenge.schema import digest_file, write_table


def command_output(command: list[str], cwd: Path) -> str:
    completed = subprocess.run(command, cwd=cwd, check=True, capture_output=True, text=True)
    return completed.stdout.strip()


def scan_python_sources(source_dir: Path, relative_paths: list[str]) -> dict[str, object]:
    syntax_failures: list[dict[str, object]] = []
    imports: set[str] = set()
    os_system_calls = 0
    for relative in relative_paths:
        if not relative.endswith(".py"):
            continue
        path = source_dir / relative
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=relative)
        except SyntaxError as error:
            syntax_failures.append(
                {
                    "file": relative,
                    "line": error.lineno,
                    "error": error.msg,
                }
            )
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
            elif (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "os"
                and node.func.attr == "system"
            ):
                os_system_calls += 1
    return {
        "python_file_count": sum(relative.endswith(".py") for relative in relative_paths),
        "syntax_failure_count": len(syntax_failures),
        "syntax_failures": syntax_failures,
        "import_roots": sorted(imports),
        "os_system_call_count": os_system_calls,
    }


def audit_source(source_dir: Path, config_path: Path, outdir: Path) -> dict[str, object]:
    if outdir.exists():
        raise FileExistsError(f"Choose a new output directory: {outdir}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    commit = command_output(["git", "rev-parse", "HEAD"], source_dir)
    status = command_output(["git", "status", "--short"], source_dir)
    if commit != config["expected_commit"] or status:
        raise ValueError("unitFinder checkout does not match the frozen clean commit")
    raw_paths = subprocess.run(
        ["git", "ls-files", "-z"], cwd=source_dir, check=True, capture_output=True
    ).stdout
    relative_paths = sorted(raw_paths.decode("utf-8").rstrip("\0").split("\0"))
    if len(relative_paths) != int(config["expected_tracked_file_count"]):
        raise ValueError("unitFinder tracked-file count differs from frozen config")
    file_rows = [
        {
            "file": relative,
            "bytes": (source_dir / relative).stat().st_size,
            "sha256": digest_file(source_dir / relative),
        }
        for relative in relative_paths
    ]
    scan = scan_python_sources(source_dir, relative_paths)
    combined_text = "\n".join(
        (source_dir / relative).read_text(encoding="utf-8", errors="replace")
        for relative in relative_paths
        if relative.endswith((".py", ".md"))
    )
    command_evidence = {
        command: bool(re.search(rf"(?<![A-Za-z0-9_.-]){re.escape(command)}(?![A-Za-z0-9_.-])", combined_text))
        for command in config["external_commands"]
    }
    if not all(command_evidence.values()):
        raise ValueError("A declared external command lacks source or README evidence")
    entrypoint = source_dir / config["entrypoint"]
    if not entrypoint.is_file():
        raise ValueError("Frozen unitFinder entrypoint is missing")
    third_party_imports = sorted(
        set(scan["import_roots"]) & {"Bio", "matplotlib", "numpy", "sklearn"}
    )
    receipt = {
        "schema_version": 1,
        "complete": True,
        "scope": "source_dependency_and_interface_audit_no_accuracy_execution",
        "official_url": config["official_url"],
        "commit": commit,
        "checkout_status": status,
        "tracked_file_count": len(relative_paths),
        "tracked_python_file_count": scan["python_file_count"],
        "packaged_python2_pyc_count": sum(relative.endswith(".pyc") for relative in relative_paths),
        "python_syntax_failure_count": scan["syntax_failure_count"],
        "python_syntax_failures": scan["syntax_failures"],
        "os_system_call_count": scan["os_system_call_count"],
        "third_party_python_imports": third_party_imports,
        "external_command_evidence": command_evidence,
        "has_dependency_manifest": any(
            Path(relative).name in {"environment.yml", "requirements.txt", "pyproject.toml", "setup.py"}
            for relative in relative_paths
        ),
        "has_automated_test_tree": any(Path(relative).parts[0].lower() in {"test", "tests"} for relative in relative_paths),
        "readme_manual_steps": {
            "shell_concatenation": "cat ../Chr*/Chr*.Type.fa" in combined_text,
            "script_modification": "Modify the ./unify/change.stat.py" in combined_text,
            "clustal_omega": "clustalo -i" in combined_text,
        },
        "config_sha256": digest_file(config_path),
        "entrypoint_sha256": digest_file(entrypoint),
        "boundary": config["boundary"],
    }
    outdir.mkdir(parents=True)
    write_table(outdir / "tracked_files.tsv", file_rows, ["file", "bytes", "sha256"])
    (outdir / "source_audit.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt
