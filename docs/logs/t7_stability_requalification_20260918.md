# T7 bounded stability requalification, 2026-09-18

The user's revised gate distinguishes storage integrity from interface speed.
USB2 alone is an I/O throughput limitation, not a scientific-accuracy NO-GO.
The qualifying commands and machine-readable inputs/results are in
`docs/evidence/t7_stability_20260918/`; the probe implementation is
`benchmarks/scripts/audit_t7_stability.py`.

At 09:13 Asia/Shanghai, `/dev/disk5s2` was mounted on `/Volumes/T7` as exFAT,
with approximately 1.9 TiB free. `system_profiler SPUSBDataType` placed the
Samsung PSSD T7 under a USB2.0 Hub; `ioreg -p IOUSB -l -w0` reported the T7
`USBSpeed` property as 3 (USB high speed). This records the current path, not
the device's intrinsic capability or scientific input integrity.

Eight predeclared, actually relevant bounded source files from Col-CEN,
Ey15-2, Mo17 and rice were read sequentially. Their total on-disk size was
5,257,735,817 bytes. Every full-file SHA-256 matched its prior independent
receipt; every independently specified size matched. The Col-CEN read subset
had a prior hash but no prior byte-size field, so its 970,235,145-byte size is
an observation, not a historical size comparison. All seven gzip files were
decompressed completely (11,879,760,176 expanded bytes in total) with no CRC,
trailer or truncation error. The Ey15 assembly is plain FASTA and was fully
hashed. The sequential full-file reads also served as non-destructive
continuous read tests; the four read subsets ranged from 970 MB to 1.15 GB.
No BAM or database is part of this next bounded source input manifest.

After read checks passed, a unique 3 GiB temporary file was written in the
existing T7 results directory in 4 MiB blocks of new random bytes. Its file
descriptor completed `fsync` and macOS `F_FULLFSYNC` (88.735 s including
write/sync). The file was closed, reopened with macOS `F_NOCACHE`, and read
sequentially (72.866 s); written and read SHA-256 both equalled
`fd4edc9e45cd1870ca34fda1c44c9a23cef407f1f661a29e6e2ac5395a2ba5cc`.
The path reported the same mounted device before and after, and the disposable
file was deleted. A subsequent mount check still showed `/dev/disk5s2` at
`/Volumes/T7`; no probe file remained. The 20-minute system-log predicate for
`I/O error`, `disk5` and `unexpected unmount` returned no matching failure
event. That predicate is supplementary and not an exhaustive hardware log
audit.

**Current gate:** the tested T7 connection is usable for bounded scientific
accuracy, M1/M2 development and controlled-collapse work whose actual input
files have passed checks. Each newly enrolled large input still requires its
own size/hash/gzip verification. The probe does not prove permanent future
reliability or validate every old file on the volume. If a future I/O error,
unmount, hash mismatch or corruption appears, new large T7 writes stop and
valid inputs are moved to stable internal/NVMe storage when feasible. Formal
runtime/scalability comparisons must use a fixed reliable storage environment
without this USB2 I/O bottleneck; no performance conclusion is drawn from the
probe timings.

After this gate, Macadamia source verification and assembly-only interval
selection resumed as a separately recorded development experiment.
