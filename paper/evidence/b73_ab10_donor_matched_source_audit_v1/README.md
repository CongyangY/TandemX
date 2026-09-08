# B73-Ab10 donor-matched source audit v1

This is a retained source/reference failure, not a TandemX accuracy result.
B73-Ab10 is unusually relevant because the v2 study explicitly states that its
PacBio HiFi library used the same high-molecular-weight genomic DNA as the CLR-
based v1 assembly. Both assemblies and the raw sequencing deposit are public.

Two independent gates nevertheless fail. First, the archival records use
different BioSamples (`SAMEA6236504` for v1 and `SAMN48109500` for the v2 HiFi
deposit), and the newer record reports collection year 2021. This does not
cancel the paper's direct same-DNA statement, but the identifier discrepancy
must be resolved before calling the archival metadata donor-identical. Second,
and decisively, the newer assembly is not high-quality collapse truth: the
authors report N-gaps predominantly within tandem-repeat arrays, slightly
smaller assembled knob180/TR-1 totals than v1, and only about 28% assembly of
the approximately 30.67-Mb Ab10 knob. A newer assembly that is explicitly
incomplete at the target arrays cannot label v1 collapse.

The public BioProject `PRJNA1254310` identifies four B73-Ab10 Sequel II runs as
HiFi/CCS-mode sequencing, but the deposited original objects are split
`subreads` BAM files and ENA exposes subread FASTQ. The four ENA conversions sum
to 855,681,258,665 bases, 57,504,646 records and 257,275,342,048 compressed
bytes. They are not a ready CCS-read input for TandemX; reconstructing CCS would
require the complete per-movie subread sets and a separate validated conversion.

Current fate: `reference_ineligible_new_assembly_incomplete`, with an additional
`input_blocked_no_public_ccs_product` warning. No large file was downloaded, no
benchmark was preregistered, and no accuracy or resource result was generated.

