# Minimal offline report example

Run inside `tandemx-dev` after installing the repository:

```bash
tandemx simulate toy --outdir /tmp/tandemx-simple-input
tandemx run --reads /tmp/tandemx-simple-input/reads.fa \
  --assembly /tmp/tandemx-simple-input/assembly.fa -o /tmp/tandemx-simple-report
tandemx validate --project /tmp/tandemx-simple-report
```

Open `/tmp/tandemx-simple-report/report.html` locally. The assembly-length
normalization warning is expected because no independent genome size was given.
Repeat the same `run` command to exercise fingerprint-validated automatic resume.
Add `--force` to rerun, or use `--no-resume` to refuse existing outputs.

Reads-only example:

```bash
tandemx run --reads /tmp/tandemx-simple-input/reads.fa \
  -o /tmp/tandemx-simple-reads-only
```

The catalogue/report are available; absolute abundance is unavailable until a
genome size is supplied. Missing assembly/abundance is not plotted as zero.

Recovery is a separately enrolled research PoC, not an automatic patch stage.
Its interface is documented in `docs/recovery_and_reporting_plan.md` and
`docs/file_formats.md`. No original assembly is modified.
