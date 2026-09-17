"""Prepare isolated pinned competitor workdirs and immutable run manifests."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONTROL = ROOT / "benchmarks/assembly_hor_comparator_v1/control/control.chr1.fasta"
EXPECTED = "60e65dd68af2af8774a872ea25656225714a01415d840b4aec59cfb09a9d1079"
OUT = ROOT / "benchmarks/assembly_hor_comparator_v1/run_manifests"
CEN = ROOT / "benchmarks/competitor_envs/cendetecthor"
TIDE = ROOT / "benchmarks/competitor_envs/tidecluster"
CEN_SOURCE = Path("/private/tmp/tandemx_cendetecthor_20260916/pipeline")
TIDE_SOURCE = Path("/private/tmp/tandemx_tidecluster_20260916")


def write_manifest(path: Path, name: str, command: list[str], cwd: Path, scratch: Path) -> None:
    path.write_text(json.dumps({"stages": [{"name": name, "command": command,
                                          "cwd": str(cwd), "scratch_dir": str(scratch)}]}, indent=2) + "\n")


def main() -> None:
    if hashlib.sha256(CONTROL.read_bytes()).hexdigest() != EXPECTED:
        raise ValueError("Frozen assembly input changed")
    if not (CEN_SOURCE / "Snakefile").is_file() or not (TIDE_SOURCE / "TideCluster.py").is_file():
        raise FileNotFoundError("Pinned competitor source checkout missing")
    OUT.mkdir(parents=True, exist_ok=True)
    cen_work = CEN / "work/assembly_hor_v1"
    tide_work = TIDE / "work/assembly_hor_v1"
    for work in (cen_work, tide_work):
        if work.exists():
            raise FileExistsError(f"Refusing to overwrite a competitor run: {work}")
        (work / "fasta").mkdir(parents=True)
        shutil.copyfile(CONTROL, work / "fasta/control.chr1.fasta")
    (cen_work / "config").mkdir()
    shutil.copyfile(CEN_SOURCE / "config/config.yaml", cen_work / "config/config.yaml")
    cen_path = ":".join((str(CEN / ".venv/bin"), str(CEN / ".conda-r/bin"),
                         str(CEN / ".conda-tools/bin"),
                         "/private/tmp/tandemx_cendetecthor_20260916/stringdecomposer/bin",
                         os.environ.get("PATH", "/usr/bin:/bin")))
    cen_command = ["/usr/bin/env", f"PATH={cen_path}", f"MPLCONFIGDIR={CEN / 'work/mpl'}",
                   "TMPDIR=/private/tmp", "PYTHONHASHSEED=0", "OPENBLAS_NUM_THREADS=1",
                   str(CEN / ".venv/bin/snakemake"), "--snakefile", str(CEN_SOURCE / "Snakefile"),
                   "--cores", "1", "--shared-fs-usage", "input-output", "persistence",
                   "software-deployment", "software-deployment-cache", "sources",
                   "storage-local-copies", "--config", f"workdir={cen_work}/", "consFile=no",
                   "CHRcons=control.chr1.fasta", "HORdet_threads=1", "windowSize=1000"]
    write_manifest(OUT / "cendetecthor.json", "cendetecthor_full_pipeline", cen_command,
                   cen_work, cen_work)
    tide_path = ":".join((str(TIDE / ".conda-clustering/bin"),
                          str(TIDE / ".local-bin"), os.environ.get("PATH", "/usr/bin:/bin")))
    tide_base = ["/usr/bin/env", f"PATH={tide_path}", "OPENBLAS_NUM_THREADS=1",
                 str(TIDE / ".venv/bin/python"), str(TIDE_SOURCE / "TideCluster.py")]
    tide_file = str(tide_work / "fasta/control.chr1.fasta")
    write_manifest(OUT / "tidehunter.json", "tidehunter",
                   tide_base + ["tidehunter", "-c", "1", "-f", tide_file,
                                "-pr", "control", "-T", "-p 100 -P 2000 -c 2 -e 0.25",
                                "--max_memory", "1"], tide_work, tide_work)
    write_manifest(OUT / "tidecluster.json", "tidecluster_clustering",
                   tide_base + ["clustering", "-c", "1", "-f", tide_file,
                                "-pr", "control", "-m", "100"], tide_work, tide_work)


if __name__ == "__main__":
    main()
