# Quantify TSV partial-write failure and bounded repair

On 2026-09-16, a local fault injection against the prior
`tandemx.quantify.mvp.write_copy_number` used a patched `Path.write_text` that
wrote the first 12 bytes, then raised `OSError(5, "injected EIO after partial
write")`. The call propagated `OSError 5`, but the public `copy_number.tsv`
already existed and contained exactly `b'family_id\tmo'`. A failed write could
therefore leave a truncated final-named TSV in a fresh output directory.
This was a software fault injection on the internal disk, not an observed T7
write or a biological result.

The repair writes the unchanged TSV bytes into a uniquely named temporary file
in the destination directory, flushes and `fsync`s it, and uses `os.replace`
only after the write succeeds. `finally` removes an uncommitted temporary file
after write, sync, replace, or `KeyboardInterrupt` failures. No estimator,
threshold, normalization, input selection, or TSV field was changed.

Focused fault tests inject partial `write` EIO, `fsync` EIO/ENOSPC,
`os.replace` EIO, and `KeyboardInterrupt`. They verify nonzero propagation,
no new final output, and no residual temporary file; a replace failure leaves
an older completed target byte-identical. A successful call writes the same
header and leaves no temporary file. The existing quantify unit tests also
pass: `conda run -n tandemx-dev pytest -q
tests/unit/test_quantify_output_atomic.py tests/unit/test_quantify_mvp.py`
gave **32 passed**.
The complete repository suite then passed **940 tests** in `tandemx-dev`.

This is an atomic single-file output fix, not stage-level checkpointing. A
`SIGKILL` or an unhandled terminating signal can still leave an orphaned
temporary file because Python cleanup cannot run; the final `copy_number.tsv`
is not replaced by that orphan. A
failed rerun in a reused output directory may still leave a prior completed
`copy_number.tsv`; callers must inspect process exit status and run configuration
or use a fresh output directory. Recovery/fingerprint gates remain separate
work. No T7 file was modified.
