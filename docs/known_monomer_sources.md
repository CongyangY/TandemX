# Source-backed repeat queries

These records provide small, post hoc biological checks. They are not discovery
inputs, a complete plant repeat catalogue, or a denominator for genome-wide
family recall. A source-cultivar sequence need not occur unchanged in another
cultivar. A miss cannot be relabelled a tool false negative without evidence of
presence in the tested material. All original records and retrieval hashes are
retained under `paper/evidence/known_repeat_sources`.

| Accession | Deposited material | Included query | Evidence boundary |
| --- | --- | --- | --- |
| AF078922.1 | Maize Seneca 60, chromosome 9 | CentC, 156 bp | One historical variant, not a Mo17 truth catalogue |
| X16095.1 | Barley Donetsky A. | HvT01, bases 1–118 of 208 bp | Explicit monomer feature; not a Morex donor |
| X55642.1 | Rice Cigalon | Deposited 358-bp repeat record | Native unit boundaries unvalidated; record lists unpublished manuscript |
| AF058902.1 | Indica rice IR-BB21 | Excluded | 639-bp multiunit RCS2/CentO clone, not a 155-bp monomer |

CentC's original [sequence study](https://doi.org/10.1073/pnas.95.22.13073),
the original [RCS2 clone study](https://doi.org/10.1073/pnas.95.14.8135), and
the [Morex gap study](https://doi.org/10.1111/pbi.13816) supply biological context.
The [X55642 record](https://www.ebi.ac.uk/ena/browser/view/X55642.1) supports
only the deposited query and associated metadata.

```bash
python -m benchmarks.scripts.curate_known_monomers \
  --records-dir paper/evidence/known_repeat_sources/records \
  --outdir /path/to/new/known_repeat_queries
pytest -q tests/unit/test_known_monomer_sources.py
```

The curator requires each complete EMBL/FASTA pair, exact accession/version,
source species/material, published sequence MD5, and retrieval SHA-256. EMBL
and FASTA bases must agree; the 118-bp extraction must match the source feature.
It writes `known_monomers.fa`, `source_units.tsv` and `curation_receipt.json` to
a new directory. Failed source validation produces no curated bank. No motif
matching result or completed accuracy assessment is implied by curation alone.

Actual consensus outputs from multiple tools can be evaluated with the same
independent cyclic edit-distance criterion:

```bash
python -m benchmarks.scripts.evaluate_known_repeats_across_tools \
  --known /path/to/known_monomers.fa \
  --query-id HvT01_X16095_1_1_118 \
  --catalog tandemx=/path/to/tandemx/monomers.fa \
  --consensus-table trf=/path/to/trf/normalized_arrays.tsv \
  --consensus-table tidehunter=/path/to/tidehunter/normalized_arrays.tsv \
  --material Morex \
  --evidence-boundary "HvT01 donor differs from Morex" \
  --outdir /path/to/new/evaluation
```

The output is source-query recovery only. The command does not establish that
the historical sequence is present unchanged in the tested material, and it
does not interpret unmatched predictions as false positives.
