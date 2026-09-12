# Real 30x recovery inventory

The machine-readable inventory is `recovery_manifest_20260913.json`. All endpoint sizes are estimated downloads; sampled 30x FASTQ requires roughly the listed target bases plus staging.

K30076 correction: `block_01_attempt_003/block_receipt.json` is 0 B, so the block is invalid and cannot be reused. Nonzero chunk files are only recovery evidence until a regenerated receipt passes aggregate record, base-count and hash validation.

Recovery order: validate YSD56 receipts; repair K30076; finish V14167; enroll S245 then HN51. Direct no-proxy access is mandatory by default. Proxy use needs a receipted <=500 MB A/B test, stable >=2x acceleration and the 50 GB project cap.

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
