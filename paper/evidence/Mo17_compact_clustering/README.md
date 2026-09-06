# Mo17 85,663-candidate compact-index replay

Source `9d6ace8`, one fixed-order execution per variant, same rounded candidate
exports from the 1.129-Gb run. All family and membership values agree exactly:
28,586 families, output SHA-256
`f0d1e88a78345707bcca3c40c948b7050a474e6b8201966543e1693b731b7039`.

| Index | Clustering seconds | Child seconds | Peak RSS MiB |
| --- | ---: | ---: | ---: |
| Python dictionary postings | 172.4202 | 174.3350 | 361.5000 |
| Python packed uint64 postings | 199.7638 | 202.6802 | 247.8594 |

The measured peak decreases 31.44%, while clustering time increases 15.86%.
Both changes are retained. Child RSS includes imports, candidate loading and
output encoding. Other jobs were active; these are diagnostic measurements,
not publication timing rankings. Serialized candidate scores have exported
precision. This is narrower than full live-pipeline parity and says nothing
about biological family truth or superiority over external software.

`archive_manifest.json` verifies the copied compact records against the T7 run.
Large result payloads and frozen source/native snapshots remain at the recorded
T7 paths. The later native-index implementation was not used in this replay.
