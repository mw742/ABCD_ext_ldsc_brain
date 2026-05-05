#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
export LC_ALL=C LANG=C

die(){ echo "[err] $*" >&2; exit 2; }
diag(){ echo "[diag] $*" >&2; }

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT_DIR="${SCRIPT_DIR}"
PATHS_ENV="${ROOT_DIR}/_paths.env"
[[ -s "$PATHS_ENV" ]] || die "Missing _paths.env at $PATHS_ENV"

# shellcheck source=/dev/null
source "$PATHS_ENV"

[[ -s "${PHENO_IMG_REGENIE}" ]] || die "Missing PHENO_IMG_REGENIE TSV: ${PHENO_IMG_REGENIE}"

OUT_LIST="${ROOT_DIR}/IMG_COLS_ALL.txt"

diag "Extracting phenotype columns from header…"
# 1) Get header, split into columns skipping FID/IID
mapfile -t COLS < <(head -n1 "${PHENO_IMG_REGENIE}" | tr -d '\r' | awk -F'\t' '{for (i=3;i<=NF;i++) print $i}')

# 2) Keep only columns ending with _baseline|_mean|_slope
mapfile -t CANDIDATES < <(printf "%s\n" "${COLS[@]}" | grep -E '_(baseline|mean|slope)$' || true)
[[ "${#CANDIDATES[@]}" -gt 0 ]] || die "No *_baseline/mean/slope columns detected."

# 3) Group by IDP base (remove the suffix) and keep stable first 70 bases
declare -A seen_base
bases_ordered=()
while IFS= read -r c; do
  base="${c%_baseline}"; base="${base%_mean}"; base="${base%_slope}"
  if [[ -z "${seen_base[$base]:-}" ]]; then
    seen_base[$base]=1
    bases_ordered+=("$base")
  fi
done < <(printf "%s\n" "${CANDIDATES[@]}")

if (( ${#bases_ordered[@]} < 70 )); then
  diag "[warn] Found only ${#bases_ordered[@]} IDP bases (<70). Will use all of them."
fi

# Cap to first 70 bases (or fewer if fewer exist)
want=$(( ${#bases_ordered[@]} < 70 ? ${#bases_ordered[@]} : 70 ))
declare -A keep_base
for ((i=0;i<want;i++)); do keep_base["${bases_ordered[$i]}"]=1; done

# 4) Emit exactly baseline, mean, slope for each kept base if present (preserve input order)
> "$OUT_LIST"
for base in "${bases_ordered[@]}"; do
  [[ -n "${keep_base[$base]:-}" ]] || continue
  for suf in baseline mean slope; do
    cand="${base}_${suf}"
    if printf "%s\n" "${CANDIDATES[@]}" | grep -Fxq "$cand"; then
      echo "$cand" >> "$OUT_LIST"
    else
      diag "[warn] Missing ${cand} in header; skipping."
    fi
  done
done

N=$(wc -l < "$OUT_LIST")
diag "Wrote ${N} phenotype columns to ${OUT_LIST}"

# 5) Quick TSV integrity check (constant #fields)
diag "Sanity checking TSV field counts…"
bad=$(awk -F'\t' 'NR==1{n=NF; next} NF!=n{print NR}' "${PHENO_IMG_REGENIE}" | head -n 1 || true)
if [[ -n "$bad" ]]; then
  diag "[warn] TSV has at least on

