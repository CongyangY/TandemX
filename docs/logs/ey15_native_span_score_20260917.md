# Ey15-2 native-read span diagnostic, 2026-09-17

The edit archive and score protocol were committed as `e4e7760` before
this scoring run. The source package was committed as `39b1f2e`; all nine
edited assembly sequences and seven original native HiFi reads retain their
frozen hashes. The scorer independently re-verifies every edit, source
context, original FASTQ, and context PAF before scoring. A clean score replay
matched all three output files byte for byte.

The source-guided 1,024-bp flank rule found unique near-exact left and right
hits in 7/7 original read molecules. Their read-array spans were 3,244,
3,240, 3,246, 3,252, 3,244, 3,241, and 3,248 bp (median 3,244 bp).
Each flank-derived boundary agreed within 10 bp with a CIGAR projection from
the read's independent spanning context alignment. The intact assembly
interval is also 3,244 bp. Using the predeclared 5% rule gives a 163-bp
discordance threshold.

For each edited assembly, the scorer measures its array span directly from
the edited FASTA using the same two natural flanks and only then checks the
measured span against the independently verified edit ledger. A hostile-review
check caught an initial implementation that read the assembly span from the
ledger; that implementation was replaced before committing the score. The
technical classifications did not change, and the final code/output replay
is byte-identical.

Among the nine edited **development cases**, the length rule classified the
eight exact injected deletions as discordant (8/8) and the intact control as
supported (1/1), with no abstentions. The smallest engineered deletion is
811 bp, roughly five times the threshold. These nine cases reuse one array
and the same seven reads; their classifications are correlated. The result
therefore only demonstrates that this source-guided span comparison detects
these large artificial deletions. It does not establish sensitivity near the
threshold, performance on unselected arrays, or performance on independent
donors.

The frozen M2 route evaluated **0/9** cases: the catalogue's 420-bp period
exceeds its 300-bp cap. Monomer boundaries and physical missing-copy truth
are unverified. The assembly and reads are attributed to the same published
sample 9994, but exact extraction pairing is not separately established;
biological accuracy remains **NOT EVALUATED**. This is no evidence of M2
superiority or a rescued Gate B. The 9-case truth is the injected assembly
bp delta only.

Evidence: `benchmarks/controlled_collapse/ey15_native_span_evidence_20260917/`
(`read_trims.json`, `per_case.jsonl`, `summary.json`). The primary case
denominator is 8 injected deletions plus one intact control, not seven
reads times nine assemblies. No TandemX production command was altered.
