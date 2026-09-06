# Nipponbare AGIS1.0 reference identity and evidence limits

Verified2026-09-06 through the original paper and NCBI's versioned assembly API.
GenBank **GCA_034140825.1**, assemblyAGIS1.0, cultivarNipponbare/isolateAGIS-1.0,
belongs toPRJNA953663. The database expects12 chromosomes and385,710,679 bp.
The paired RefSeq **GCF_034140825.1 differs**: NCBI explicitly records added
chromosomesMT, Pltd andB1. We acquire the original GCA version and keep organellar
and other reference additions separate. They are not interchangeable denominators.

The assembly BioSampleSAMN36344332 records50-day leaf material collected
2023-06-16. The included HiFi runSRR25241090 belongs toSAMN36368305, which records
15-day leaf material collected2022-03-07. Cultivar, isolate and study agree;
the same donor plant/extraction is **unresolved**. Metadata should not be
rewritten to imply exact specimen matching.

The original [Molecular Plant paper](https://doi.org/10.1016/j.molp.2023.08.003)
(first page, accessible in the [author-uploaded primary paper](https://www.researchgate.net/publication/372967140_A_complete_assembly_of_the_rice_Nipponbare_reference_genome))
states that the chr9 terminal45S rDNA array was filled using estimated copy
number and must be treated as model sequence. Exclude that array from exact
assembly-copy truth. The study's CENH3 ChIP, CentO annotations and published FISH
comparison are useful orthogonal context, but do not supply a complete real
repeat-family truth set. Supplementary coordinates still need independent curation.

The BioProject free-text HiFi yield(1043 Gb) conflicts with the included ENA
HiFi-labelled WGS run's verified32,966,159,623 bp. Its current ten-row ENA snapshot
has one PacBio entry; neither the free text nor the assembly API's220 coverage
field is used as measured HiFi/nuclear depth. Full-file QC observed1,785,885 reads.
The precise protocol/whole-study yield discrepancy remains unresolved.

```bash
python -m benchmarks.scripts.fetch_ncbi_reference \
  --accession GCA_034140825.1 --bioproject PRJNA953663 \
  --source-paper https://doi.org/10.1016/j.molp.2023.08.003 \
  --outdir /path/to/references/Nipponbare_AGIS1_GCA034140825v1 --max-download-gb 0.3
pytest -q tests/unit/test_fetch_ncbi_reference.py
```

Official source records:
[assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_034140825.1/),
[raw-read BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN36368305),
[assembly BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN36344332),
[BioProject](https://www.ncbi.nlm.nih.gov/bioproject/953663).
The downloader archives the exact API response, reports and official checksums.
