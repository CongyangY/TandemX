#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-/tmp/tandemx-tidehunter-import}"
TIDEHUNTER="${2:-TideHunter}"

if ! command -v "${TIDEHUNTER}" >/dev/null 2>&1; then
  echo "TideHunter executable not found: ${TIDEHUNTER}" >&2
  exit 2
fi

tandemx simulate toy --outdir "${OUTDIR}/simulated" --seed 17 --num-reads 120
"${TIDEHUNTER}" \
  -t 1 -f 2 -p 20 -P 500 -m 20 -c 2 \
  -o "${OUTDIR}/tidehunter.tsv" \
  "${OUTDIR}/simulated/reads.fa"
tandemx import tidehunter \
  --input "${OUTDIR}/tidehunter.tsv" \
  --reads "${OUTDIR}/simulated/reads.fa" \
  --min-support-reads 2 \
  --outdir "${OUTDIR}/imported"
tandemx validate --project "${OUTDIR}/imported"
tandemx quantify \
  --reads "${OUTDIR}/simulated/reads.fa" \
  --catalog "${OUTDIR}/imported/monomers.fa" \
  --genome-size 7744 \
  --outdir "${OUTDIR}/quantify"

echo "TideHunter import workflow complete: ${OUTDIR}"
