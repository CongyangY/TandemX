# Figure 3: operating range across plant long-read datasets

This is the manuscript candidate for Figure 3. It replaces the rejected
editorial draft in `figure3_operating_range/` and is intentionally limited to
the manuscript's operating-range claim.

Panels A--B use the ten complete public HiFi libraries that passed source-byte,
checksum, gzip and full FASTQ checks. Panel C shows two TandemX outputs in
closely matched 11.4--12.7 Mb nested inputs from wheat, barley, rye and oat;
the metrics are given separate aligned scales. Panel D connects only nested
inputs from the same library and includes the 1.17-Gb Morex run. Neither panel
is a cross-species accuracy comparison, a biological-replication analysis or a
runtime/memory ranking.

Regenerate with:

```bash
MPLCONFIGDIR=/private/tmp/tandemx_mpl conda run -n tandemx-dev python build_figure3_v2.py
```

Outputs are editable SVG/PDF plus a 350-dpi PNG. `panel_source.tsv` records
the authoritative source table and panel-level evidence boundary.
