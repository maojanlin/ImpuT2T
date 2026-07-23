#!/usr/bin/env bash
# Extract per-donor chromosome FASTAs from Index-zone AGC archives into database/chr_split/.
#
# Usage (from repository root):
#   database/extract_chr_pangenome.sh [donor_list.tsv] [chr ...]
#
# If no chromosomes are listed after the donor list, extract all of
# chr1–22,X,Y,M (skipping any chromosome whose AGC is missing).
#
# Default donor list: snakemake/subsample_lists/sample_id.POP-10.log
# Requires: agc (bioconda), archives in database/chr_agc/pangenome_<chr>.agc

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
AGC_DIR="${REPO_ROOT}/database/chr_agc"
OUT_ROOT="${REPO_ROOT}/database/chr_split"
DONOR_LIST="${1:-${REPO_ROOT}/snakemake/subsample_lists/sample_id.POP-10.log}"
if [[ $# -gt 0 ]]; then
  shift
fi

if [[ $# -gt 0 ]]; then
  CHROMS=("$@")
else
  CHROMS=()
  for i in $(seq 1 22); do CHROMS+=("chr${i}"); done
  CHROMS+=("chrX" "chrY" "chrM")
fi

if [[ ! -f "$DONOR_LIST" ]]; then
  echo "Donor list not found: $DONOR_LIST" >&2
  exit 1
fi

mapfile -t DONORS < <(grep -v '^[[:space:]]*#' "$DONOR_LIST" | grep -v '^[[:space:]]*$' | sed 's/\r$//')

for chr in "${CHROMS[@]}"; do
  agc_file="${AGC_DIR}/pangenome_${chr}.agc"
  if [[ ! -f "$agc_file" ]]; then
    echo "Skip ${chr}: missing ${agc_file}" >&2
    continue
  fi
  out_dir="${OUT_ROOT}/${chr}"
  mkdir -p "$out_dir"
  echo "Extracting ${chr} ..."
  for donor in "${DONORS[@]}"; do
    out_fa="${out_dir}/${donor}.fasta"
    if agc getset "$agc_file" "$donor" > "$out_fa" 2>/dev/null; then
      :
    else
      rm -f "$out_fa"
    fi
  done
done

echo "Done. FASTAs under ${OUT_ROOT}/"
