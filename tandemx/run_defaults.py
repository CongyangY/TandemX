"""Strict advanced configuration and recorded defaults for ``tandemx run``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import yaml


# These options only tune the existing pipeline surface.  Inputs and output
# locations stay on the command line so an advanced file cannot silently run
# against a different sample or overwrite a different result directory.
ADVANCED_OPTION_TYPES: dict[str, type | tuple[type, ...]] = {
    "genome_size": int,
    "haploid_depth": float,
    "max_reads": int,
    "max_read_bases": int,
    "kmer_backend": str,
    "steps": (str, list),
    "min_period": int,
    "max_period": int,
    "top_periods": int,
    "threads": int,
    "discovery_method": str,
    "clustering_method": str,
    "family_audit": str,
    "cluster_identity": float,
    "single_copy_kmers": str,
    "read_error_rate": float,
    "disable_quality_correction": bool,
    "profile": bool,
}


class _StrictConfigLoader(yaml.SafeLoader):
    """Safe YAML loader that refuses duplicate or non-string option keys."""

    def construct_mapping(self, node: yaml.MappingNode, deep: bool = False) -> dict[str, object]:
        mapping: dict[str, object] = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, str):
                raise ValueError("--config keys must be strings")
            if key in mapping:
                raise ValueError(f"Duplicate --config key: {key}")
            mapping[key] = self.construct_object(value_node, deep=deep)
        return mapping


def _type_name(expected: type | tuple[type, ...]) -> str:
    if expected == (str, list):
        return "string or list"
    return expected.__name__ if isinstance(expected, type) else " or ".join(item.__name__ for item in expected)


def _valid_type(value: object, expected: type | tuple[type, ...]) -> bool:
    if expected is float:
        return type(value) in {int, float}
    if expected is int:
        return type(value) is int
    if expected is bool:
        return type(value) is bool
    return isinstance(value, expected)


def _explicit_options(argv: list[str]) -> set[str]:
    explicit: set[str] = set()
    for token in argv:
        if token.startswith("--"):
            explicit.add(token[2:].split("=", 1)[0].replace("-", "_"))
        elif token == "-o":
            explicit.add("outdir")
    return explicit


def apply_advanced_config(args: argparse.Namespace) -> None:
    """Apply a YAML object while keeping explicit command-line values first."""
    config_path = getattr(args, "config", None)
    if config_path is None:
        return
    try:
        loaded = yaml.load(Path(config_path).read_text(encoding="utf-8"), Loader=_StrictConfigLoader)
    except OSError as exc:
        raise ValueError(f"Cannot read --config file: {config_path}") from exc
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML in --config file: {config_path}") from exc
    if not isinstance(loaded, dict):
        raise ValueError("--config must contain one mapping of advanced run options")
    if not all(isinstance(key, str) for key in loaded):
        raise ValueError("--config keys must be strings")
    unknown = sorted(set(loaded) - set(ADVANCED_OPTION_TYPES))
    if unknown:
        raise ValueError("Unknown --config key(s): " + ", ".join(unknown))
    explicit = _explicit_options(getattr(args, "_argv", []))
    for key, value in loaded.items():
        expected = ADVANCED_OPTION_TYPES[key]
        if not _valid_type(value, expected):
            raise ValueError(f"--config key '{key}' must have type {_type_name(expected)}")
        if key in explicit:
            continue
        if key == "steps":
            value = ",".join(value) if isinstance(value, list) and all(isinstance(x, str) for x in value) else value
            if not isinstance(value, str):
                raise ValueError("--config key 'steps' must be a comma-separated string or a list of strings")
        if key == "single_copy_kmers":
            value = Path(value)
        setattr(args, key, value)


def write_automatic_defaults(outdir: Path, *, genome_size: int | None, genome_size_source: str,
                             threads: int, kmer_backend: str, config_path: Path | None) -> None:
    """Record fixed automatic choices without suggesting they are optimized."""
    payload = {
        "schema_version": 1,
        "diagnostic_k": 21,
        "diagnostic_k_selection": "fixed_existing_default_not_an_optimality_claim",
        "threads": threads,
        "kmer_backend": kmer_backend,
        "genome_size_bp": genome_size,
        "genome_size_source": genome_size_source,
        "warnings": (
            ["assembly_total_length_is_a_provisional_normalization_proxy; it can underestimate a target or haploid genome and does not establish ploidy"]
            if genome_size_source == "assembly_total_length_provisional" else []
        ),
        "advanced_config": str(config_path) if config_path is not None else None,
    }
    (outdir / "automatic_defaults.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
