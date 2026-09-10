# Ey15-2 targeted recovery: bounded negative proof of concept

The candidate was generated from the historical assembly and original HiFi
reads only. Newer sequence was used after `candidate_lock.json` was complete.
Previously known family-level newer-assembly abundance values mean this was
recovery-input separation, not a wholly unobserved validation set.

- Selected families: TXF000002 and TXF000154.
- Six historical localized intervals; 20 of 48 eligible flank windows.
- One 5,681-bp unpolished span candidate with 26 dual-anchor reads.
- Old/candidate/newer proxy span lengths: 5,679 / 5,681 / 5,679 bp.
- Frozen full-catalogue localized repeat bases: 679 / 644 / 679 bp.
- Candidate-excluded read support: 97 / 97 / 97 of 120 reads.
- Validated repeat recovery successes: **0**. Stop algorithm expansion.

The 35-bp difference is a localization result on an unpolished read, not a
physical deletion estimate. Family abundance is not assigned to this locus.
No original assembly was patched and no biological HOR order is inferred.

`recovery/recovery_assessment.json` carries the final decision; the original
candidate lock/table remains unchanged. A `partially_resolved` candidate is
not a validated missing-repeat repair. All technical failure/continuation
fates are retained separately. `manifest.json` binds every selected artifact
to its full T7 source. Large recruited-read evidence and full-read PAF remain
on T7 and are hash-bound by the original candidate/mapping receipts.

Reproduction entry points are `tandemx.recovery.poc` and the three benchmark
scripts `extract_recovery_validation_reads.py`, `validate_recovery_proxy.py`
and `remap_recovery_validation.py`. The source configurations preserve the
long-template preparation, short-template cost continuation and resource-only
streaming continuation. See `docs/recovery_and_reporting_plan.md` for the
predefined biological rules, input restrictions and stop decision.
