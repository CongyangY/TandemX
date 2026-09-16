# Pipeline interruption cleanup, 2026-09-16

Scope: production `tandemx run` orchestration only; no scientific estimator changed.

The command runner now starts each step in a separate POSIX process group. If
the parent is interrupted while waiting for the child or forwarding its output,
it sends SIGTERM to the group, waits up to five seconds, and escalates to
SIGKILL if needed. The pipeline records exit status 130, marks the step
unvalidated, writes the summary/report/manifest, and exits 130. Partial stage
files remain visible for diagnosis; stage-level resume still requires validated
outputs and a matching fingerprint. This does not add intra-stage checkpoints.

Verification: `conda run -n tandemx-dev pytest -q
tests/unit/test_pipeline_process_cleanup.py tests/unit/test_pipeline_summary.py`
returned 5 passed in 11.30 seconds. The focused checks cover normal stdout and
stderr forwarding, a real SIGINT delivered to a wrapper with a live child,
absence of that child after the wrapper exits, and the interruption receipt.
They do not establish cleanup behavior after a hard power failure or SIGKILL
of the parent.
