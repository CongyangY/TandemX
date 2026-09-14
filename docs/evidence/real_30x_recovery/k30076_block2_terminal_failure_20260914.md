# K30076 block-2 terminal recovery state — 2026-09-14

This record supersedes the earlier `planned_not_started` description for
K30076 block 2. It was written from internal-disk evidence after the T7 mount
disappeared. No T7 path was accessed while preparing this record.

## Attempt 001: download completed, independent integrity gate failed

The frozen block-2 runner downloaded 34/34 chunks and initially emitted an
aggregate receipt for 845,098 reads and 15,027,572,584 bases
(10.1181243578x). That aggregate is **not acceptable**. Independent full-file
validation stopped at chunk 3, spots 4,615,107--4,640,106. Its file size still
matched the receipt at 370,177,772 bytes, but the receipt SHA-256
`2e3fb14389a5674ec1a62c32db13e0476e5aa62aa1c8501bd86968af10734867`
had changed to the same independently observed value on two reads,
`06bdc011167b5ddbf301da9fca458a79a4264a2ec86b74194cf6fe653cce4e56`.
Gzip decompression also failed with invalid distance/literal-length-distance
codes. The first two chunks passed both SHA-256 and gzip checks. Later chunks
were not promoted after the first mismatch.

## Attempt 002: direct service failure

Only the invalid chunk was requested in a new isolated directory. Direct
`vdb-dump` with proxy variables cleared returned zero of 25,000 required rows
and wrote `Failed to call external services.` to stderr. The retained partial
gzip was 20 bytes and the failed receipt was 947 bytes. Attempt 001 was not
overwritten.

## Attempt 003: T7 I/O failure and disconnect

Only the same invalid chunk was requested in another isolated directory. The
user-authorized proxy exception applied only to that chunk; no proxy address or
credential was written to logs or receipts. During streaming gzip output, the
write failed with `OSError: [Errno 5] Input/output error`. The `/Volumes/T7`
mount then disappeared before a pending failure receipt could be written.
Session 89647 exited 1, no downloader remained, and no automatic mount or retry
was attempted. At the recorded post-failure snapshot the Samsung PSSD T7 still
enumerated at the USB layer, but the T7 filesystem was not mounted. The current
readback state of attempts 001 and 002 is therefore not verified.

## Terminal boundary

K30076 block 2 is `invalid_no_acceptable_aggregate`. Block 3 remains
`planned_not_started`; no block-3 payload or frozen comparator was started.
The common 30x source and comparator remain blocked until all three frozen
blocks are mounted, independently hash-valid, semantically valid, and combined
under the byte-identical input rule. The frozen comparator status is
`blocked_not_run`. These failures are technical acquisition and storage
failures, not biological or comparator results. No production algorithm or
frozen parameter was changed.

Internal raw evidence retained under `tmp/`:

- `tmp/k30076_block1_recovery_validation_20260914.json`, SHA-256
  `014d314061fd445c3cfd74038edfbfe2cb7c6030acd5999a1cf677884f1704a6`.
- `tmp/k30076_block2_integrity_failure_20260914.json`, SHA-256
  `1682809c43e3c86d88f9ab188279a995eb8b844338dfbb85380baceb432231b2`.
- `tmp/k30076_block2_attempt002_failure_20260914.json`, SHA-256
  `f0652ed972df0ad11dfc27c939b91d88031be0ce27d60356caf02f3a66b7af43`.
- `tmp/k30076_block2_attempt003_t7_io_failure_20260914.json`, SHA-256
  `8b3821bc71ae09330b3ed23acb0555187805dd907161d366ba0cd24d8780236e`.
