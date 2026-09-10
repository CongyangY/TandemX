# TandemX manuscript v5 supplementary tables

## Supplementary Table S1. Public plant HiFi libraries used to establish the operating range of TandemX

All rows passed complete-file checksum and FASTQ validation. The table describes
input libraries and does not treat species, pooled samples or technical libraries
as biological accuracy replicates. The complete machine-readable table, including
file identities and validation receipts, is available in
`paper/tables/input_cohort.tsv`.

| Species | Material | Run accession | Reads | Total bases (Gb) | Median read length (kb) | Read N50 (kb) | GC (%) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| *Zea mays* | Mo17 | SRR15447419 | 407,670 | 5.625 | 13.455 | 13.968 | 45.922 |
| *Arabidopsis thaliana* | Col-0N | ERR6210723 | 933,904 | 14.647 | 15.514 | 15.663 | 36.762 |
| *Oryza sativa* | Nipponbare | SRR25241090 | 1,785,885 | 32.966 | 17.489 | 19.159 | 43.502 |
| *Hordeum vulgare* | Morex | ERR4659246 | 896,701 | 19.525 | 21.499 | 21.875 | 44.278 |
| *Secale cereale* | Lo7 | ERR15194059 | 4,576,975 | 77.092 | 16.757 | 16.896 | 45.575 |
| *Triticum aestivum* | Chinese Spring | SRR28200549 | 1,500,000 | 24.954 | 15.985 | 16.579 | 45.957 |
| *Avena sativa* | Victoria | ERR10422581 | 398,850 | 7.346 | 18.137 | 18.295 | 44.172 |
| *Arabidopsis thaliana* | Col-0R | ERR8666127 | 982,810 | 17.747 | 16.635 | 18.027 | 37.392 |
| *Arabidopsis thaliana* | Ey15-2R | ERR8666125 | 837,586 | 18.637 | 21.506 | 22.587 | 37.146 |
| *Glycine soja* | YSD56 | SRR28726931 | 2,617,227 | 44.193 | 16.224 | 17.474 | 34.749 |

## Supplementary Table S2. Unified abundance comparison on newly generated planted-read data

Values are mean absolute relative error across three independently simulated
datasets per condition. Every positive-truth family enters the mean, including
families not recovered by a method. The competitive-mapping baseline uses the
TandemX-discovered catalogue and is therefore an abundance baseline rather than
an independent discovery method.

| Condition | TandemX | SRF k=151 | SRF k=101 | Competitive mapping |
| --- | ---: | ---: | ---: | ---: |
| Clean | 0.008361 | 0.045623 | 0.049840 | 0.000359 |
| 1% substitutions | 0.193880 | 0.057697 | 0.052819 | 0.000261 |
| 0.1% insertion + 0.1% deletion | 0.047374 | 0.055464 | 0.049672 | 0.000251 |
| 2% unit divergence | 0.349075 | 1.000000 | 0.162351 | 0.000127 |
| 10% positive reads | 0.003704 | 0.035118 | 0.047998 | 0.000394 |
| Shared-fragment background | 0.013280 | 0.058794 | 0.046549 | 0.000437 |

The complete native-catalogue, unmatched-abundance and negative-read attribution
endpoints are retained in `paper/evidence/srf_formal_unified_v1/summary.tsv` and
`paper/evidence/srf_formal_unified_v1/condition_summary.tsv`.

## Supplementary Table S3. Discovery yield and abundance-stratified recovery across nested read depths

Values are median [range] across three independently generated 10-Mb validation
genomes. Recall was calculated against all 55 planted families. Low-, medium-
and high-abundance tiers contained 18, 18 and 19 families with 20, 80 and at
least 200 planted copies, respectively.

| Depth | Raw candidates | Unique candidates | Operational families | Overall recall | Low recall | Medium recall | High recall |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5x | 105 [99--128] | 36 [33--44] | 34 [32--44] | 0.618 [0.582--0.800] | 0.500 [0.333--0.667] | 0.722 [0.500--0.778] | 0.842 [0.684--0.947] |
| 1x | 214 [211--245] | 48 [46--55] | 46 [44--51] | 0.836 [0.800--0.927] | 0.611 [0.556--0.889] | 0.889 [0.889--0.889] | 1.000 [0.947--1.000] |
| 2x | 435 [429--463] | 57 [55--62] | 53 [52--54] | 0.964 [0.945--0.982] | 0.889 [0.833--0.944] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] |
| 5x | 1,086 [1,031--1,088] | 70 [70--73] | 55 [55--57] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] |
| 10x | 2,121 [2,082--2,128] | 79 [78--84] | 56 [55--57] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] |
| 20x | 4,240 [4,192--4,261] | 101 [101--107] | 57 [57--58] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] |
| 30x | 6,343 [6,302--6,460] | 119 [117--125] | 58 [57--58] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] | 1.000 [1.000--1.000] |

## Supplementary Table S4. Marginal family discovery and catalogue stability between adjacent read depths

Values are median [range] across three genomes. `New operational families`
counts higher-depth representatives without a one-to-one match at cyclic global
edit identity at least 0.90. `New planted families` counts truth families first
recovered at the higher depth. A seed passed when new-family fraction was below
0.05 and family-level Jaccard was above 0.95.

| Transition | New operational families | Operational families per added 1x | New planted families | Planted families per added 1x | New-family fraction | Family Jaccard | Passing seeds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5x--1x | 12 [7--12] | 24.000 [14.000--24.000] | 12 [7--12] | 24.000 [14.000--24.000] | 0.261 [0.137--0.273] | 0.739 [0.727--0.863] | 0/3 |
| 1x--2x | 7 [3--8] | 7.000 [3.000--8.000] | 7 [3--8] | 7.000 [3.000--8.000] | 0.132 [0.056--0.154] | 0.868 [0.846--0.944] | 0/3 |
| 2x--5x | 3 [1--4] | 1.000 [0.333--1.333] | 2 [1--3] | 0.667 [0.333--1.000] | 0.055 [0.018--0.070] | 0.945 [0.930--0.982] | 1/3 |
| 5x--10x | 0 [0--1] | 0.000 [0.000--0.200] | 0 [0--0] | 0.000 [0.000--0.000] | 0.000 [0.000--0.018] | 1.000 [0.982--1.000] | 3/3 |
| 10x--20x | 1 [1--2] | 0.100 [0.100--0.200] | 0 [0--0] | 0.000 [0.000--0.000] | 0.018 [0.017--0.035] | 0.982 [0.965--0.983] | 3/3 |
| 20x--30x | 0 [0--1] | 0.000 [0.000--0.100] | 0 [0--0] | 0.000 [0.000--0.000] | 0.000 [0.000--0.017] | 1.000 [0.983--1.000] | 3/3 |

The earliest depth with two consecutive all-seed passing transitions was 20x.
Complete per-seed values and family-level recovery states are provided in
`paper/evidence/discovery_saturation_validation_v1/`.
