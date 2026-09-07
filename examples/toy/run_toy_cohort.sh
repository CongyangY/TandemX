#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-examples/toy/cohort-results}"
TANDEMX_CMD="${TANDEMX_CMD:-tandemx}"
read -r -a TANDEMX <<< "${TANDEMX_CMD}"

run_sample() {
  local sample_id="$1"
  local copies="$2"
  local genome_size="$3"
  local sample_dir="${OUTDIR}/${sample_id}"
  local haploid_depth
  haploid_depth=$(python -c "print(f'{(120 * 1200) / ${genome_size}:.6f}')")

  "${TANDEMX[@]}" simulate toy \
    --outdir "${sample_dir}/simulated" \
    --seed 42 \
    --num-reads 120 \
    --read-length 1200 \
    --background-length 2000 \
    --monomer-lengths 566,350 \
    --copies "${copies}" \
    --error-rate 0.005

  "${TANDEMX[@]}" discover \
    --reads "${sample_dir}/simulated/reads.fa" \
    --outdir "${sample_dir}/discover" \
    --min-monomer-len 300 \
    --max-monomer-len 700 \
    --min-support-reads 3 \
    --min-repeat-span 600

  "${TANDEMX[@]}" quantify \
    --reads "${sample_dir}/simulated/reads.fa" \
    --catalog "${sample_dir}/discover/monomers.fa" \
    --genome-size "${genome_size}" \
    --haploid-depth "${haploid_depth}" \
    --outdir "${sample_dir}/quantify"

  "${TANDEMX[@]}" locate \
    --assembly "${sample_dir}/simulated/assembly.fa" \
    --catalog "${sample_dir}/discover/monomers.fa" \
    --copy-number "${sample_dir}/quantify/copy_number.tsv" \
    --window-size 500 \
    --step-size 250 \
    --outdir "${sample_dir}/locate"

  "${TANDEMX[@]}" compare \
    --copy-number "${sample_dir}/quantify/copy_number.tsv" \
    --arrays "${sample_dir}/locate/arrays.bed" \
    --outdir "${sample_dir}/compare"
}

mkdir -p "${OUTDIR}"

# The same seed fixes the two monomer sequences. Different copy counts create
# two abundance profiles. Genome sizes include the 2-kb background and spacers.
run_sample sample_a 9,7 7744
run_sample sample_b 5,11 8880

MANIFEST="${OUTDIR}/samples.tsv"
cat > "${MANIFEST}" <<EOF
sample_id	monomers	copy_number	comparison
sample_a	sample_a/discover/monomers.fa	sample_a/quantify/copy_number.tsv	sample_a/compare/assembly_vs_read_cn.tsv
sample_b	sample_b/discover/monomers.fa	sample_b/quantify/copy_number.tsv	sample_b/compare/assembly_vs_read_cn.tsv
EOF

"${TANDEMX[@]}" cohort \
  --manifest "${MANIFEST}" \
  --outdir "${OUTDIR}/cohort"

"${TANDEMX[@]}" validate --project "${OUTDIR}/cohort"

echo "Toy cohort workflow complete: ${OUTDIR}/cohort"
