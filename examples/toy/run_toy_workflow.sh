#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-examples/toy/results}"
TANDEMX_CMD="${TANDEMX_CMD:-tandemx}"
read -r -a TANDEMX <<< "${TANDEMX_CMD}"

SIMULATED="${OUTDIR}/simulated"
DISCOVER="${OUTDIR}/discover"
QUANTIFY="${OUTDIR}/quantify"
LOCATE="${OUTDIR}/locate"
COMPARE="${OUTDIR}/compare"
PROBE="${OUTDIR}/probe"
VISUALIZE="${OUTDIR}/visualize"

mkdir -p "${OUTDIR}"

"${TANDEMX[@]}" simulate toy \
  --outdir "${SIMULATED}" \
  --seed 42 \
  --num-reads 120 \
  --read-length 1200 \
  --background-length 2000 \
  --monomer-lengths 566,350 \
  --copies 9,7 \
  --error-rate 0.005

"${TANDEMX[@]}" discover \
  --reads "${SIMULATED}/reads.fa" \
  --outdir "${DISCOVER}" \
  --min-monomer-len 300 \
  --max-monomer-len 700 \
  --min-support-reads 3 \
  --min-repeat-span 600

# For the default toy simulator parameters, the read-generating source length is:
# background_length + 566*9 + 350*7 + two 100 bp spacers = 7744 bp.
GENOME_SIZE=7744
HAPLOID_DEPTH=$(python - <<'PY'
print(f"{(120 * 1200) / 7744:.6f}")
PY
)

"${TANDEMX[@]}" quantify \
  --reads "${SIMULATED}/reads.fa" \
  --catalog "${DISCOVER}/monomers.fa" \
  --genome-size "${GENOME_SIZE}" \
  --haploid-depth "${HAPLOID_DEPTH}" \
  --outdir "${QUANTIFY}"

"${TANDEMX[@]}" locate \
  --assembly "${SIMULATED}/assembly.fa" \
  --catalog "${DISCOVER}/monomers.fa" \
  --copy-number "${QUANTIFY}/copy_number.tsv" \
  --window-size 500 \
  --step-size 250 \
  --outdir "${LOCATE}"

"${TANDEMX[@]}" compare \
  --copy-number "${QUANTIFY}/copy_number.tsv" \
  --arrays "${LOCATE}/arrays.bed" \
  --outdir "${COMPARE}"

"${TANDEMX[@]}" probe \
  --catalog "${DISCOVER}/monomers.fa" \
  --assembly "${SIMULATED}/assembly.fa" \
  --copy-number "${QUANTIFY}/copy_number.tsv" \
  --arrays "${LOCATE}/arrays.bed" \
  --outdir "${PROBE}"

"${TANDEMX[@]}" visualize \
  --catalog "${DISCOVER}/monomers.fa" \
  --copy-number "${QUANTIFY}/copy_number.tsv" \
  --comparison "${COMPARE}/assembly_vs_read_cn.tsv" \
  --probes "${PROBE}/probes.rank.tsv" \
  --fish "${PROBE}/in_silico_fish.tsv" \
  --outdir "${VISUALIZE}"

"${TANDEMX[@]}" validate --project "${OUTDIR}"

echo "Toy workflow complete: ${OUTDIR}"
