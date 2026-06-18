#!/usr/bin/env python3
"""Run a step-level TandemX pipeline benchmark."""

from __future__ import annotations

import argparse

from tandemx.pipeline import add_pipeline_arguments, config_from_args, run_pipeline


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_pipeline_arguments(parser)
    args = parser.parse_args()
    config = config_from_args(args)
    _, exit_status = run_pipeline(config)
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
