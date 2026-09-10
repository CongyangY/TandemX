# Read-occupancy input audit (2026-09-10)

## Scope and conclusion

This is a read-only audit of existing inputs and code. It does not introduce an
estimator, select a threshold, inspect a new holdout result, or rerun a
benchmark. A valid family read-occupancy quantity is

```
sum over reads of union_bp(accepted family intervals on that read) / sum read length
```

with each base assigned at most once within each family on a read. Distinct
non-overlapping intervals from different families may each contribute to their
respective family totals; a policy may instead exclude multi-family reads, but
that is a conservative mapping rule rather than a mathematical requirement.
The numerator is not the number of reads with a hit. The latter is a support
rate and is length- and threshold-dependent.

The repository has one directly reusable implementation of this numerator:
`benchmarks/scripts/evaluate_macadamia_ont_occupancy.py:summarize_paf`.
It streams query-grouped PAF, forms a per-read interval union, and uses the
additional conservative policy of excluding a read when accepted primary
mappings name more than one family; it records accepted query bp, library
fraction, and source hashes. It is an orthogonal,
direction-only mapping audit, not a general TandemX abundance baseline.

## Existing interfaces and what they can provide

| Material | Join/fields | Occupancy readiness | Boundary |
| --- | --- | --- | --- |
| `discover/candidate_reads.tsv` | `read_id`, `candidate_id`, `read_start`, `read_end`, `repeat_span_bp` | Intervals exist, but no `family_id`. A family mapping must be supplied through the candidate/monomer-membership and family tables. Intervals need unioning per `(read_id, family_id)`. | Candidate calls are not competitive alignments and can overlap; `repeat_span_bp` rows must not be summed directly. |
| `candidate_monomers.fa`, `monomer_membership.tsv`, `families.tsv` | candidate/monomer/family provenance | Can construct the missing candidate-to-family crosswalk when its identifiers are complete and one-to-one. | Must reject unmapped or multiply assigned candidates rather than silently choose a family. |
| Challenge `truth_arrays.tsv` | `read_id`, `start`, `end`, `family_id` | Direct oracle for sampled per-read family union bp. `benchmarks/challenge/simulate.py` writes it. | Truth is evaluation-only and must never enter a TandemX command. |
| Abundance read-condition `manifest.json` | sampled repeat bp and actual source/observed coverage | Supplies a second, independent simulation oracle. `benchmarks/abundance/run_stream_quantify.py` explicitly confines it to scoring. | It is not a read-interval table, so it cannot validate overlapping predicted intervals by itself. |
| ONT/HiFi competitive PAF | query coordinates, target family, primary tag | `summarize_paf` already implements the correct union numerator and multi-family exclusion. | Requires PAF sorted/grouped by query; primary status and supplied mapping filters determine the measured occupancy. |

The current public quantifier (`tandemx/quantify/mvp.py`) is diagnostic-k-mer
depth normalized by a genome-size/depth denominator. It has no whole-read
alignment occupancy estimator. `tandemx/quantify/multik.py` and
`joint_multik.py` provide experimental multi-k diagnostic-k-mer extrapolation
and joint read moments; these are neither whole-read alignments nor
multi-representative mapping models. Their documented uncertainty excludes
calibration, background specificity and biological sampling bias.

## Reusable code and known error sources

* Reuse `union_length` and the query-grouped `summarize_paf` flow only after
  retaining its safeguards: primary-filtering, block/identity validation,
  the chosen multi-family policy, PAF query grouping, and file hashes.
* `candidate_reads.tsv` is a discovery product. Its coordinate interval may be
  duplicated, nested, harmonic, or overlap another candidate. It therefore
  needs a family crosswalk and per-read union; a hit-rate denominator is wrong.
* Mapping occupancy has repeat-template composition, divergence, read-error,
  mappability, primary-alignment and multi-family-exclusion bias. Multiplying a
  library fraction by a genome size yields an abundance-like bp scale only under
  explicit calibration; it is not physical copy-number truth.
* The existing Macadamia audit applies one frozen global HiFi mapping-efficiency
  factor and labels the output direction-only. It does not establish a
  family-specific correction, a new estimator, or an interchangeability claim
  across sequencing platforms.

## Prior consumption and remaining holdout space

The following seeds/results are already consumed and must not be treated as new
holdout material: abundance genomes `s4101` onward in the existing baseline
and development archives; classifier development seeds 5601--5603; the first
classifier heldout 5701--5703, which became development after its failed gate;
and depth-gated validation seeds 5801--5803 (consumed after the earlier
"reserved" historical status); and factorial-scale seeds 6301--6303 and
6401--6403. Seeds 7401--7403 are rejected for reuse. Existing heldout result
directories were deliberately not opened for this audit.

The specifically retained unconsumed option is factorial-scale seeds 3201--3203.
A new simulation holdout must still preregister its configuration, write a new
generation receipt before output inspection, and keep `truth_arrays.tsv` and
sampled-repeat metadata unavailable to the command being evaluated. The 3201--
3203 option is not authorization to alter frozen estimator rules.

## Small T7 pilot material

A small already-consumed, hash-bound input suitable only for plumbing/parity
checks is:

* `/Volumes/T7/Codex/TandemX/results/abundance_baseline_v1_20260906/genomes/s4101/`
  with `genome.fa` SHA-256
  `26a4179bda66336256ead3bbaa1fe981b1b396a43c35ed53ee549e55961855fd`,
  `catalogue.fa` SHA-256
  `582d27c3c9fb67482e7c67a0feaf1d72b860915bf7e993d92feb7af1683963b4`, and
  `truth_copy_number.tsv` SHA-256
  `64232ca7c69757b2c6522a9a67424e414fa0424013487051cc10dd51b16d90da`, as
  declared by its `manifest.json`.

This is a known-catalogue simulation input, not a new holdout and it does not
by itself contain read-level intervals. The T7 real pilot
`results/ERR6210723_prefix1000_baseline/discover/candidate_reads.tsv` has read
intervals but no family column, so it is unsuitable for occupancy calculation
until a validated candidate-to-family crosswalk and the exact input-read
length denominator are enrolled. Neither should be used to claim biological
abundance or to reopen frozen benchmark conclusions.
