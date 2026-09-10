# Competitor gap audit — 2026-09-10

## Scope and conclusion

This is a read-only audit of three supplied reviewer-style suggestion notes and
the current repository/T7 evidence. It checks claims against original papers
and official implementations; it does not repeat the notes' scores, publication
advice, or unverified citations. No software was installed, no sequence was
downloaded, and no benchmark was run.

The central positioning is supportable only at the **family-resolved audit**
endpoint: TandemX links a discovered operational family, its read-derived
abundance estimate, assembly-localized representation, and recorded uncertainty.
It must not claim that read-first satellite discovery, generic read--assembly
k-mer completeness, or local extra-long-tandem-repeat assembly assessment is
new. SRF is the one direct read-first satellite workflow that needs a bounded
held-out comparison. Merqury/KAT, TandemTools, and AniAnn's are important
conceptual/feature comparators, but do not share a single fair de novo,
family-level endpoint with TandemX.

## Verified capabilities and comparable endpoints

| Method | Verified capability | Fair endpoint against TandemX | Endpoint that is not interchangeable |
| --- | --- | --- | --- |
| [SRF](https://github.com/lh3/srf) | Builds abundant satellite motifs from high-occurrence k-mers in short reads, accurate long reads, or high-quality contigs. Its documented workflow maps elongated motifs back to source sequences, resolves competing motif mappings heuristically, and reports mapped-bp abundance. It may emit HORs rather than minimal monomers. | Read-first satellite motif/array recovery, motif-to-truth sequence correspondence, read interval/base recovery, native mapped-bp abundance error, complete workflow resources. | A minimal TandemX monomer missed by SRF is not automatically an SRF failure if SRF recovered a valid containing HOR; a family/HOR-aware correspondence rule is required. SRF's mapped-bp fraction is not a calibrated physical copy-number truth. |
| [Merqury](https://github.com/marbl/merqury) | Reference-free assembly QV, k-mer completeness and spectra analyses from assembly plus high-accuracy read k-mers; optionally hap-mer/phasing analyses. The primary paper explicitly frames it as assembly-level assessment. | Whole-assembly k-mer completeness/QV or spectra, using the same accurate read set and assembly. | It does not discover repeat families from reads or assign a missing k-mer signal to a de novo repeat family. Its global spectra cannot be scored as TandemX family recall or family-specific deficit accuracy. |
| [KAT](https://github.com/EarlhamInst/KAT) | Counts, compares and visualizes k-mer spectra and GC composition for reads and assemblies; its paper demonstrates pairwise read--assembly composition/QC. | Whole-assembly spectra/assembly composition check on identical k-mer inputs. | It is a k-mer QC toolkit, not a de novo satellite-family caller or a locus/family deficit estimator. A spectra-cn signal is not a family label. |
| [TandemTools](https://github.com/ablab/TandemTools) | TandemMapper maps long reads to assembled extra-long tandem repeats; TandemQUAST assesses and can polish ETR assemblies. | Assembly-locus quality evidence after a target ETR interval and matching long reads are supplied. | It is not a read-first de novo satellite/family discovery or whole-sample family-abundance method. It should not be forced into the synthetic family-discovery scorecard. |
| [AniAnn's](https://github.com/marbl/anianns) | Assembly/FASTA annotation of tandem/satellite arrays using alignment-free ANI/modimizers; optional class assignment requires a supplied satellite k-mer database. | Array boundaries and optional supplied-library classification on a fully resolved assembly. | It does not provide a read-first de novo abundance estimate. Its required class database makes its classification mode unlike TandemX de novo family construction. |

Primary literature confirms the endpoint distinctions: SRF reconstructs satellite
units/HORs from accurate reads or assemblies ([Zhang *et al.*, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760446/));
Merqury evaluates reference-free assembly quality/completeness from read and
assembly k-mers ([Rhie *et al.*, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7488777/));
KAT provides read/assembly k-mer QC ([Mapleson *et al.*, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5408915/));
TandemTools evaluates assembled ETRs with long reads ([Mikheenko *et al.*, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7355294/));
and AniAnn's annotates arrays in FASTA assemblies ([McCartney *et al.*, 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12873891/)).

## Existing evidence found

The supplied notes repeatedly state that SRF lacks a direct benchmark. That is
now incomplete rather than wholly correct:

* `docs/srf_workflow.md`, `benchmarks/scripts/run_srf_pilot.py`, and
  `paper/evidence/srf_pilot/` document the upstream seven-stage workflow:
  KMC, `kmc_dump`, SRF, `enlong`, minimap2, `paf2bed`, and `bed2abun`.
  Native files, commands, hashes, status and stage resources are retained.
* `/Volumes/T7/Codex/TandemX/results/srf_workflow_pilot_v1_20260906/`
  has four successful 500,000-base planted-truth runs (clean and 0.1% indel,
  each at `ci20` and `ci100`). Its `validation.json` says explicitly that it
  has no timing repetitions and no calibrated copy-number claim. All four
  recovered 70/70 scored arrays; these are pilot observations, not a general
  advantage claim.
* `/Volumes/T7/Codex/TandemX/results/srf_workflow_development_v1_20260906/`
  records 32 development invocations, including 11 native failures caused by
  empty eligible-k-mer dumps. The guarded rerun preserves this distinction
  (`no_eligible_kmers`) instead of fabricating an empty SRF catalogue.
* `docs/comparator_matrix.md` already recognizes SRF as the essential
  read-first comparator and explicitly records that k=151/ci20/ci100 is a
  development condition, not sufficient tuning. It also states that TandemTools
  is applicable only once a target array exists. No local/T7 evidence was found
  for a Merqury, KAT, TandemTools, or AniAnn's execution, which is appropriate
  to their non-overlapping primary endpoints but leaves literature/feature
  positioning to add to the manuscript.

Thus the real gap is **not** basic SRF installation or a zero-versus-nonzero
toy run. It is a prespecified, held-out, family/HOR-aware SRF comparison with a
native abundance-error endpoint and an honest statement of where methods have
different output units.

## Minimal fair SRF pilot (design only)

1. **Hold-out inputs.** Use new, deterministic simulator seeds and families not
   used by the existing SRF development suite. Retain the same observable
   accurate-read FASTA for both tools, plus the truth only for downstream
   scoring. Include one clean and one modest-error condition, and one
   low-abundance condition; do not add real data or change TandemX rules.
2. **Development before hold-out.** The current runner is fixed at `k=151`, and
   the existing 500,000-base toy input does not supply an empirical genome
   coverage from which a new cutoff can be frozen. Calibrate a finite SRF
   k/count selection set first, including the paper-motivated k=101 condition
   only after the runner supports it. Freeze the selected condition and all
   correspondence rules before generating or opening the one-time hold-out.
   Record native assertion failures and guarded empty-count outcomes separately.
3. **Endpoint normalization.** Preserve raw SRF FASTA/BED/abundance output.
   Score (a) any overlapping tandem interval/base recovery, (b) motif or
   HOR-to-truth correspondence with a declared cyclic-sequence criterion, and
   (c) per-native-motif mapped-bp fraction versus simulated repeat bp. Report
   the correspondence table rather than collapsing SRF HORs into TandemX
   monomers. TandemX is scored at its existing family and interval outputs;
   neither method receives truth-derived motifs as input.
4. **Resources and failure accounting.** Include all seven SRF stages and
   TandemX's complete equivalent workflow. Preserve per-stage and end-to-end
   wall time, direct-child peak RSS, retained output bytes, commands, hashes,
   exit status, and `no_eligible_kmers`/native-failure states. Repeat only after
   the pilot establishes machine variability; repeats are technical, not
   biological replicates.

### Bounded resources

The pilot is limited to the existing toy-scale envelope: at most 500,000 input
bases per dataset, one thread, KMC's observed minimum `-m2` allocation, and no
downloads or new dependencies. Four held-out datasets x two frozen SRF
conditions is at most eight SRF workflows; any paired TandemX invocation uses
the same FASTA and does not introduce a new algorithm or threshold. Stop and
report resource/infrastructure failure rather than expanding memory, changing
the benchmark after inspection, or treating a failed executable as zero
accuracy.

## Manuscript-facing implications

State the following narrowly: TandemX provides a family-resolved audit view
that combines read-first family discovery, read-derived family abundance and
assembly representation with recorded confidence. Cite SRF for read-first
satellite reconstruction and mapped-bp abundance; cite Merqury/KAT for global
k-mer assembly evaluation; cite TandemTools for locus-level ETR assessment; and
cite AniAnn's for assembly FASTA array annotation. Do not claim that any of
these tools is directly outperformed until the corresponding shared endpoint
has been run and reported.
