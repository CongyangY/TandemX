# Real 30x recovery inventory

The original machine-readable interruption inventory is
`recovery_manifest_20260913.json`. The direct-database and local-file re-audit is
`recovery_manifest_v2_20260913.json`; see `cohort_reaudit_20260913.md` for the
evidence boundary. Endpoint byte counts are download sizes, while sampled 30x
FASTA/FASTQ requires roughly the listed target bases plus staging.

At the 2026-09-13 audit, K30076's original
`block_01_attempt_003/block_receipt.json` was 0 B and could not be reused.
Nonzero chunk files were only recovery evidence until the later independently
validated block-1 recovery described below.

## K30076 terminal recovery update — 2026-09-14

The block-1 recovery subsequently passed its independent gate, but block 2 did
not. Although attempt 001 downloaded 34/34 chunks, independent validation found
a stable SHA-256 mismatch and gzip corruption in chunk 3. Attempt 002 failed
with zero rows from direct `vdb-dump`; attempt 003 failed with a T7 `Errno 5`
write error followed by disappearance of the mount. Block 3 was not started,
and the K30076 common 30x input and comparator remain blocked. See
`k30076_block2_terminal_failure_20260914.md` and
`k30076_block2_terminal_failure_20260914.json`. The current status is a
technical acquisition/storage failure and must not be interpreted as a
comparator or biological result.

Recovery order: retain the validated YSD56 receipts; repair K30076; reconcile the
remaining V14167 TandemX/TRF cells after its three TideHunter partitions closed;
then enroll S245 and HN51. Direct no-proxy access is mandatory by default. Proxy
use needs a receipted <=500 MB A/B test, stable >=2x acceleration and the 50 GB
project cap.

## Read-only recovery integrity gate

`benchmarks/scripts/audit_real_30x_recovery.py` audits copied or mounted
explicit material blocks without writing into them. It never recursively scans
an entire mounted results tree. It distinguishes a complete,
hash-bearing receipt (`valid`) from `recoverable_invalid`, `missing`,
`corrupt`, and a non-mounted/inaccessible root (`blocked`). A nonzero chunk set
beside a zero-byte receipt is always `recoverable_invalid`, never valid; it can
only become reusable after a later, independently regenerated aggregate
receipt verifies record/base totals and hashes.

When a receipt declares `required_files` or `required_artifacts`, the gate also
requires each declared file to be nonzero and checks any declared SHA-256. It
does not invent tool-specific required outputs for an older receipt schema; a
receipt lacking the portable completion, positive record/base-total, and hash
fields remains invalid rather than being silently accepted.

Reproduce the current mount check, writing its report to the internal repo:

```bash
python benchmarks/scripts/audit_real_30x_recovery.py \
  --root /Volumes/T7/Codex/TandemX/results \
  --block peanut_k30076_30x_comparison_v1_20260912/blocks/block_01_attempt_003 \
  --out tmp/real_30x_recovery_audit.json
```
