# Historical prioritization source and native-output re-audit

This is a read-only re-audit of the formal historical-assembly SRF/ordinary-
mapping prioritization endpoint after the T7 remount. It does not run a native
comparator and does not change the production algorithm.

The three enrolled HiFi FASTQ files were read in full on 2026-09-13. Their
current SHA-256 values match the completed source/QC receipts. The two frozen
TandemX catalogues and primary family-metric tables match their compact archive
manifests. All six comparator executable hashes also match the formal simulation
environment. No external source file needs to be downloaded again.

Formal historical output is still absent for both SRF k settings in both
species and for the Ey15 ordinary-mapping arm. Macadamia retains a hash-valid
HiFi-to-43-family PAF and occupancy table, but those files belong to the frozen
orthogonal direction-only audit. That run used eligible-only 250-kb templates,
`-N 5` and additional filtering; it lacks a complete command/resource/source
receipt and differs from the formal simulation mapping rule. It must not be
relabeled as the missing formal historical arm.

The current documents do not constitute a runnable frozen historical protocol.
The formal SRF protocol explicitly excludes historical prioritization, accepts
one uncompressed simulation FASTA, and uses a 180-second stage cap. It does not
define historical FASTQ conversion/order, Macadamia two-file handling, real-
library normalization, SRF family/HOR correspondence, zero/unmatched states, or
which of the two incompatible competitive-mapping rules controls this endpoint.
Executing now would be post-hoc method construction.

The exact verified compressed-read payload is 35,547,861,768 bytes. Including
the two frozen catalogues, a separate staging copy would require at least
35,548,536,680 bytes before any KMC temporary database, dumped k-mers, SRF
catalogue/mapping output, logs or receipts. The internal filesystem had
39,066,935,296 bytes free at audit, leaving only 3,518,398,616 bytes after that
minimum copy. KMC temporary/output demand has no frozen bound. The main task
subsequently reported `fsck_exfat` exit 0 and a healthy T7 remount, restoring a
writable large-output route with sufficient reported free space. Therefore the
minimal external acquisition is zero bytes and storage is no longer the active
blocker; the missing historical-specific execution contract remains the reason
no formal comparator was launched.

`source_reaudit.json` records every current hash, the native-output inventory,
the exact protocol gaps and the minimum unblocking actions. Any later endpoint
must be described as agreement with a retrospective newer-assembly proxy, not
physical-copy or population accuracy.
