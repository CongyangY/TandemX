# Original known-repeat records

`records/` contains four official ENA EMBL/FASTA pairs and their original
retrieval receipt. `archive_manifest.json` records byte/SHA parity of each
copy from T7. Accession versions, materials and source sequences are checked
by `benchmarks/scripts/curate_known_monomers.py` and its source-fixture tests.
See `docs/known_monomer_sources.md` for individual evidence limits, exclusions,
source papers and the repeat-query extraction command. These are not a complete
repeat truth catalogue and do not independently establish presence in the
materials being benchmarked.
