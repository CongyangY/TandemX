# MorexV3 reference acquisition and file QC

The original IPK MorexV3 pseudomolecule FASTA was acquired from data
publication DOI `10.5447/ipk/2021/3` under CC BY 4.0. The completed file is
4,296,032,540 bytes and its SHA-256 is
`54c98a04d13ff97350f5f3a5bfa45ac395ad640df8bb1f7598eca4e7edb437c1`,
identical to the checksum published on the archived source page.

Streaming FASTA QC found eight unique records (`chr1H`–`chr7H` and `chrUn`),
4,225,605,719 sequence bases, 1,353,994 `N` bases and no other ambiguity codes.
The IPK server ignored byte-range requests; the successful transfer first
verified the retained 2.126-GB prefix against the new full response and then
downloaded only the remaining suffix. The final receipt records
`downloaded_complete_after_verified_prefix`.

`archive_manifest.json` re-hashes the compact source page, plan and completion
receipt and independently re-hashes the external FASTA. The 4.30-GB FASTA is
kept at the configurable data root and is intentionally not copied into Git.
The raw-read and reference projects share Morex study context, but this does not
prove an identical DNA donor or exact satellite-array copy truth.
