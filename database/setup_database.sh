#!/usr/bin/env bash
# Download and extract ImpuT2T reference data into database/.
#
# Usage (from repository root):
#   database/setup_database.sh              # full panel (all chromosomes + donors)
#   database/setup_database.sh --test       # smole-test: download AGCs; extract only POP-10 chr20+chr22
#   database/setup_database.sh --help
#
# Requires: wget (or curl), agc (bioconda: conda install -c bioconda agc)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DB="${REPO_ROOT}/database"
HQ="${DB}/high_quality_set"
CHR_AGC="${DB}/chr_agc"
CHR_SPLIT="${DB}/chr_split"
EXTRACT_SCRIPT="${DB}/extract_chr_pangenome.sh"

ZENODO_AGC_URL="https://zenodo.org/records/14854401/files/human472.agc?download=1"
HPRC2_AGC_BASE="https://genome-idx.s3.amazonaws.com/agc/hprc2"

MODE="full"   # full | test
DONOR_LIST="${REPO_ROOT}/snakemake/subsample_lists/sample_id.log"
# Chromosomes to download into chr_agc/ (always full panel by default / --test).
CHROMS=()
for i in $(seq 1 22); do CHROMS+=("chr${i}"); done
CHROMS+=("chrX" "chrY" "chrM")
# Chromosomes to extract into chr_split/ (empty = same as CHROMS).
EXTRACT_CHROMS=()

usage() {
  cat <<'EOF'
Download ImpuT2T reference databases into database/.

Usage:
  database/setup_database.sh [options]

Options:
  (default) Full extraction: all chromosomes (chr1–22,X,Y,M) + sample_id.log donors.
  --test    Smoke-test data: download all chromosome AGCs, extract
            POP-10 donors for chr20 and chr22 (matches config.POP-10.yaml).
  --donor-list PATH
            Override donor list for chr_split extraction.
  -h, --help
            Show this help.

Examples:
  database/setup_database.sh
  database/setup_database.sh --test
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --test|--smoke)
      MODE="test"
      DONOR_LIST="${REPO_ROOT}/snakemake/subsample_lists/sample_id.POP-10.log"
      # Still download the full AGC set; only extract smoke-test chromosomes.
      EXTRACT_CHROMS=(chr20 chr22)
      shift
      ;;
    --full)
      # Kept for compatibility; same as default.
      MODE="full"
      DONOR_LIST="${REPO_ROOT}/snakemake/subsample_lists/sample_id.log"
      EXTRACT_CHROMS=()
      shift
      ;;
    --donor-list)
      DONOR_LIST="$2"
      shift 2
      ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

download() {
  local url="$1"
  local out="$2"
  if [[ -f "$out" ]]; then
    echo "Exists, skip: $out"
    return 0
  fi
  echo "Downloading: $out"
  mkdir -p "$(dirname "$out")"
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "$out" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --retry 3 -o "$out" "$url"
  else
    echo "Need wget or curl on PATH" >&2
    exit 1
  fi
}

need_agc() {
  if ! command -v agc >/dev/null 2>&1; then
    echo "agc not found. Install with: conda install -c bioconda agc" >&2
    exit 1
  fi
}

echo "=== ImpuT2T database setup (${MODE}) ==="
echo "Repo: ${REPO_ROOT}"
mkdir -p "$HQ" "$CHR_AGC" "$CHR_SPLIT"

need_agc

# --- 1. human472 → high_quality_set (RagTag small_pan) ---
echo ""
echo "=== [1/3] Small pangenome (Zenodo human472.agc) ==="
download "$ZENODO_AGC_URL" "${DB}/human472.agc"

SMALL_PAN_IDS=(
  100000_CHM13.pri
  110000_GRCh38.pri
  120001_CN1.pat
  120002_CN1.mat
  120003_YAO.pat
  120004_YAO.mat
  120005_KSA001.pat
  120006_KSA001.mat
  200001_HG002.pat
  200002_HG002.mat
)
for id in "${SMALL_PAN_IDS[@]}"; do
  out="${HQ}/${id}.fa"
  if [[ -f "$out" ]]; then
    echo "Exists, skip: $out"
    continue
  fi
  echo "Extracting small_pan: $id"
  agc getset "${DB}/human472.agc" "$id" > "$out"
done

# --- 2. Per-chromosome AGC archives (Index zone) ---
echo ""
echo "=== [2/3] Chromosome AGC archives (Index zone) ==="
for chr in "${CHROMS[@]}"; do
  download \
    "${HPRC2_AGC_BASE}/pangenome_${chr}.agc" \
    "${CHR_AGC}/pangenome_${chr}.agc"
done

# --- 3. Extract donors into chr_split ---
echo ""
echo "=== [3/3] Extract donors into chr_split ==="
if [[ ! -x "$EXTRACT_SCRIPT" ]]; then
  chmod +x "$EXTRACT_SCRIPT" || true
fi

if [[ ${#EXTRACT_CHROMS[@]} -eq 0 ]]; then
  EXTRACT_CHROMS=("${CHROMS[@]}")
fi
"$EXTRACT_SCRIPT" "$DONOR_LIST" "${EXTRACT_CHROMS[@]}"

echo ""
echo "=== Done ==="
echo "  high_quality_set: ${HQ}"
echo "  chr_agc:          ${CHR_AGC} ($(printf '%s ' "${CHROMS[@]}"))"
echo "  chr_split:        ${CHR_SPLIT} (extracted: $(printf '%s ' "${EXTRACT_CHROMS[@]}"))"
echo "  donor list:       ${DONOR_LIST}"
if [[ "$MODE" == "test" ]]; then
  echo ""
  echo "Next: cd snakemake && snakemake --cores 48 --configfile config.POP-10.yaml"
else
  echo ""
  echo "Next: edit snakemake/config.yaml (queries, sample), then"
  echo "      cd snakemake && snakemake --cores 48 --configfile config.yaml"
fi
