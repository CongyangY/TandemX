# AGENTS.md

## Project identity

This repository develops TandemX, an open-source computational genomics tool for large plant genomes.

TandemX is designed for tandem repeat and satellite repeat analysis in species such as wheat, rye, barley, oat, maize, and other large repetitive plant genomes.

The scientific goal is to build a read-first, assembly-aware platform that can:

1. Discover tandem/satellite repeat monomers directly from HiFi raw reads.
2. Estimate read-corrected repeat copy number using diagnostic k-mer depth.
3. Localize tandem repeat arrays on genome assemblies when assemblies are provided.
4. Compare read-based and assembly-based copy number to detect possible array collapse or under-assembly.
5. Prioritize FISH probe candidates with predicted signal regions and specificity scores.
6. Provide publication-ready visualization and reproducible benchmark results.

## Development principles

1. Do not build a simple wrapper around existing tools.
2. Do not generate fake outputs just to make examples look complete.
3. Do not silently ignore errors.
4. Do not write monolithic scripts.
5. Do not load whole large genomes or all raw reads into memory unless explicitly justified.
6. Core algorithms should be designed for streaming, chunking, and parallel execution.
7. Every public command must have documented input files, output files, and field definitions.
8. Every new function should have unit tests when practical.
9. Every command should have at least one toy example.
10. Every major module should expose uncertainty or confidence labels when the biological interpretation is uncertain.
11. Avoid overstating what the software can infer from raw reads or assemblies.

## Scientific boundaries

TandemX should not claim that HiFi reads can fully reconstruct megabase-scale satellite arrays by default.

Preferred wording:

- estimate copy number
- infer array organization
- identify candidate monomers
- anchor repeat arrays when sufficient evidence exists
- detect possible assembly under-representation

Avoid wording:

- completely resolve all tandem arrays
- precisely locate every repeat copy
- fully assemble all satellite arrays from raw reads

## Repository structure expectations

The repository should be organized as follows:

TandemX/
├── README.md
├── LICENSE
├── pyproject.toml
├── environment.yml
├── AGENTS.md
├── docs/
│ ├── design.md
│ ├── algorithms.md
│ ├── file_formats.md
│ ├── benchmark_plan.md
│ └── publication_plan.md
├── tandemx/
│ ├── **init**.py
│ ├── cli.py
│ ├── discover/
│ ├── quantify/
│ ├── locate/
│ ├── probe/
│ ├── compare/
│ ├── visualize/
│ ├── io/
│ └── utils/
├── tests/
│ ├── unit/
│ ├── integration/
│ └── data/
├── benchmarks/
│ ├── simulated/
│ ├── real_genomes/
│ ├── scripts/
│ └── results/
├── examples/
│ ├── toy/
│ ├── rye/
│ └── wheat/
├── workflows/
│ ├── Snakefile
│ └── nextflow/
└── paper/
├── outline.md
├── figures/
└── tables/

## Planned CLI

The target command-line interface is:

tandemx discover
tandemx quantify
tandemx locate
tandemx probe
tandemx compare
tandemx visualize

Each command must:

1. Validate input files.
2. Create the output directory if needed.
3. Write a run_config.yaml file.
4. Write a log file.
5. Produce documented outputs.
6. Fail with clear error messages when inputs are invalid.

## MVP scope

The first MVP should only support toy-scale data.

The MVP should include:

1. A toy HiFi-like read simulator.
2. Detection of simple tandem repeat monomers from toy reads.
3. Diagnostic k-mer based copy-number estimation on toy reads.
4. Assembly repeat density estimation on toy assembly.
5. Basic FISH probe scoring.
6. One end-to-end tutorial.

The MVP should not include:

1. Full optimization for 16 Gb genomes.
2. Complex higher-order repeat inference.
3. Interactive web dashboard.
4. Real wheat or rye analysis.
5. Rust/C++ reimplementation unless the Python prototype is stable.

## Development environment

All development and testing should be performed inside the dedicated conda environment `tandemx-dev`.

Use:

```bash
conda env create -f environment.yml
conda activate tandemx-dev
pip install -e .
pytest
```

Rules:

1. Do not install packages into the base conda environment.
2. Do not modify global shell configuration files such as `~/.bashrc`, `~/.zshrc`, or `~/.condarc`.
3. Do not install system-level packages unless explicitly requested.
4. Do not add heavy dependencies without explaining why they are necessary.
5. Keep `environment.yml` minimal during the MVP stage.
6. If a new dependency is required, update `environment.yml`, README.md, and the relevant documentation.
7. All tests should be run with Python >=3.10, preferably Python 3.11.
8. The conda environment isolates software dependencies, but it does not replace git commits, backups, or careful file management.

## Testing requirements

Use pytest for Python tests.

Every module should include:

1. Unit tests for small functions.
2. Integration tests for command-line execution.
3. Toy data smaller than 1 MB.
4. Tests for empty input, invalid input, low-complexity sequences, reverse-complement sequences, and reproducibility.

Do not modify unrelated modules when adding tests.

## Output file rules

All TSV outputs must:

1. Have a header line.
2. Use stable column names.
3. Be documented in docs/file_formats.md.
4. Avoid ambiguous field names.
5. Include confidence or warning columns when relevant.

FASTA headers should be structured and machine-readable.

BED files should be 0-based and documented.

## Visualization goals

TandemX visualization should eventually include:

1. Repeat catalogue view.
2. Genome landscape view.
3. Array inspector.
4. Read evidence view.
5. Assembly-vs-read copy-number view.
6. Sample comparison view.
7. FISH probe design view.
8. Publication-ready PDF/SVG output.

SVG text should remain editable when possible.

## Benchmark goals

Benchmark should compare TandemX with TRF, TideHunter, TRASH, and RepeatExplorer2/TAREAN when applicable.

Metrics should include:

1. Monomer length accuracy.
2. Repeat family recall.
3. False positive rate.
4. Runtime.
5. Peak memory.
6. Read-based copy-number error.
7. Assembly-collapse detection accuracy.
8. Probe specificity.
9. FISH validation success rate when experimental data are available.

## Coding style

1. Prefer clear, maintainable code over clever code.
2. Use type hints for Python functions when practical.
3. Keep functions small and testable.
4. Avoid global state.
5. Avoid hidden dependencies.
6. Use pathlib for file paths.
7. Use logging instead of print for pipeline messages.
8. Do not introduce heavy dependencies without explaining why.

## Documentation requirements

When changing public behavior, update:

1. README.md
2. docs/file_formats.md
3. docs/algorithms.md
4. relevant example commands

When adding a new command, include:

1. CLI help text.
2. Input description.
3. Output description.
4. Minimal example.
5. Test command.

## Current development priority

The immediate goal is not to analyze real large genomes.

The immediate goal is to create a clean, testable, publication-oriented MVP that can run end-to-end on toy data and establish the architecture for future large-genome optimization.

## Codex session memory

These repository-level instructions are the persistent project memory for Codex.
At the start of every new Codex conversation in this repository:

1. Read and follow this `AGENTS.md` before making changes.
2. Read `docs/current_status.md` completely before planning or changing code. Verify its Git and test claims with `git status`, `git log`, and the relevant checks because the document is a handoff snapshot.
3. Preserve the scientific boundaries above and avoid overclaiming biological inference.
4. Treat generated outputs, benchmark scratch directories, and local environment caches as non-source artifacts unless explicitly requested.
5. Use the `tandemx-dev` conda environment for development checks and run `pytest` before committing when practical.
6. Commit meaningful code, tests, and documentation together so GitHub remains the durable project record across compressed or restarted conversations.
7. Do not modify global shell configuration or install packages outside the project environment without explicit approval.
