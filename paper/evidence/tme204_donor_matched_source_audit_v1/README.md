# TME204 donor-matched source audit v1

This is a negative source-eligibility result, not a TandemX accuracy result.
The cassava TME204 study is scientifically attractive because the authors state
that 121x PacBio CLR and 42x HiFi reads were generated from the same DNA sample,
and that the CLR-Falcon/Falcon-Unzip assembly was less haplotype-resolved than
the HiFi-hifiasm assembly.  The sample relationship therefore passes the donor
match screen at the study level.

The public data inventory does not yet pass the assembly-availability gate.
The GigaDB 102193 API returned 86 files.  It includes chromosome, haplotig and
soft-masked forms of the two final HiFi+Hi-C haplotypes, but no filename
containing `CLR`, `Falcon` or `Unzip`, and no file described as the older
CLR-Falcon/Falcon-Unzip assembly.  The linked Mendeley v2 record exposes
annotation and supplementary-data folders rather than that old assembly.  The
paper reports quality statistics for the old assembly but does not provide a
public sequence accession for it.

The exact ENA runs are available and share BioSample `SAMEA8245905`, but the
combined public downloads would be about 56.14 GB compressed.  Downloading the
reads cannot repair the missing historical assembly, and rebuilding it would
change the validation question from a retrospective old-versus-new comparison
to a new assembler experiment.  No raw data were therefore downloaded and no
benchmark was preregistered or run.

Current fate: `metadata_blocked_no_public_old_assembly`.  TME204 can be
reconsidered if the published CLR-Falcon/Falcon-Unzip FASTA and its exact
sequence provenance become available.  Its high heterozygosity also requires a
separate haplotype-aware scoring design; the two HiFi haplotypes must not be
summed or treated as interchangeable truth without a frozen ploidy model.

