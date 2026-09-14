# T7 read-only volume verification — 2026-09-13

The external `T7` volume was reconnected as `/dev/disk4s2`, an ExFAT volume
mounted at `/Volumes/T7`. No process had a current working directory on the
volume when verification started.

The following read-only command was run before resuming TandemX writes:

```text
diskutil verifyVolume /Volumes/T7
```

Observed terminal result:

```text
Started file system verification on disk4s2 (T7)
Verifying file system
Volume was successfully unmounted
Performing fsck_exfat -n -x /dev/rdisk4s2
Checking volume
Checking main boot region
Checking system files
Volume name is T7
Checking upper case translation table
Checking file system hierarchy
Checking active bitmap
Rechecking main boot region
Rechecking alternate boot region
The volume T7 appears to be OK
File system check exit code is 0
Restoring the original state found as mounted
Finished file system verification on disk4s2 (T7)
```

This verifies the ExFAT filesystem structure at that time. It does not restore
missing or zero-byte files and does not replace per-file checksums, atomic
receipts, or source-to-output provenance checks. Previously invalidated
V14167/K30076 cells therefore remain invalid until their own input and receipt
gates pass.
