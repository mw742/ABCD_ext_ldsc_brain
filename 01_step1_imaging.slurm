#!/usr/bin/env bash
#SBATCH --job-name=regenie_s1_img
#SBATCH --output=regenie_s1_img_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=shared
# Must point to the folder containing _paths.env
#SBATCH --chdir=/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe

# ----------------------- SAFETY -----------------------
set -euo pipefail
IFS=$'\n\t'
export LC_ALL=C LANG=C
trap 'echo "[fatal] line $LINENO failed" >&2' ERR
diag(){ echo "[diag] $*" >&2; }
die(){  echo "[err]  $*" >&2; exit 2; }

# ----------------------- MODULES -----------------------
module purge || true
module load plink   || true
module load regenie || true

# ----------------------- _paths.env (use PWD/SLURM_SUBMIT_DIR) -----------------------
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-}"
CANDIDATES=( "${PWD}/_paths.env" )
[[ -n "$SUBMIT_DIR" ]] && CANDIDATES+=("${SUBMIT_DIR}/_paths.env")
PATHS_ENV=""
for c in "${CANDIDATES[@]}"; do
  if [[ -f "$c" ]]; then PATHS_ENV="$c"; break; fi
done
[[ -n "$PATHS_ENV" ]] || { echo "[err]  _paths.env not found. Tried: ${CANDIDATES[*]}"; exit 2; }
diag "Using PATHS_ENV=$PATHS_ENV"
# shellcheck source=/dev/null
source "$PATHS_ENV"

# Expect in _paths.env: OUTDIR, WORK, BFILE, PHENO_IMG_REGENIE, IMG_COLS_ALL
[[ -n "${OUTDIR:-}" && -n "${WORK:-}" && -n "${BFILE:-}" ]] || die "OUTDIR/WORK/BFILE not set in _paths.env"
[[ -s "${PHENO_IMG_REGENIE:-}" ]] || die "Missing pheno file: ${PHENO_IMG_REGENIE:-<unset>}"
[[ -s "${IMG_COLS_ALL:-}"      ]] || die "Missing column list: ${IMG_COLS_ALL:-<unset>}"
[[ -s "${BFILE}.bed" && -s "${BFILE}.bim" && -s "${BFILE}.fam" ]] || die "BFILE trio not found at ${BFILE}.{bed,bim,fam}"

PHENO_FILE="${PHENO_IMG_REGENIE}"
ALL_COLS="${IMG_COLS_ALL}"

# ----------------------- LIMIT TO FIRST N PHENOS -----------------------
MAX_IMG="${MAX_IMG:-ALL}"  # set to an integer to cap; use ALL for no cap
TMP_COLS="${SLURM_TMPDIR:-/tmp}/img_cols_${SLURM_JOB_ID:-$$}.lst"

if [[ "$MAX_IMG" == "ALL" ]]; then
  cp "$ALL_COLS" "$TMP_COLS"
else
  head -n "$MAX_IMG" "$ALL_COLS" > "$TMP_COLS"
fi

N_COLS=$(wc -l < "$TMP_COLS")
[[ "$N_COLS" -gt 0 ]] || die "No phenotype columns after selection (MAX_IMG=${MAX_IMG})"
PHENO_LIST_COMMA="$(paste -sd, "$TMP_COLS")"
[[ -n "$PHENO_LIST_COMMA" ]] || die "Empty pheno list after selection."


# ----------------------- OUTPUTS -----------------------
S1_OUT="${OUTDIR}/regenie_step1_imaging"
mkdir -p "$S1_OUT"
LM_PREFIX="${S1_OUT}/tmp_chunks"

# ----------------------- BUILD/REUSE LD-PRUNED TRAINING SET -----------------------
TRAIN_SNPS="${TRAIN_SNPS:-${S1_OUT}/ldtrain.prune.in}"
if [[ ! -s "$TRAIN_SNPS" ]]; then
  diag "Building LD-pruned training SNPs from $BFILE (this is done once)…"
  # First pass: fairly standard pruning (~0.2 LD), plus basic QC
  plink \
    --bfile "$BFILE" \
    --maf 0.01 \
    --geno 0.01 \
    --hwe 1e-6 midp \
    --indep-pairwise 1000 50 0.2 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --out "${S1_OUT}/ldtrain"

  # Extract pruned variants
  plink \
    --bfile "$BFILE" \
    --extract "${S1_OUT}/ldtrain.prune.in" \
    --make-just-bim \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --out "${S1_OUT}/ldtrain_keep"

  awk '{print $2}' "${S1_OUT}/ldtrain_keep.bim" > "$TRAIN_SNPS"

  N_KEEP=$(wc -l < "$TRAIN_SNPS" || echo 0)
  diag "Initial pruned count: ${N_KEEP}"
  if (( N_KEEP > 1000000 )); then
    diag "Still >1,000,000 SNPs; tightening pruning to r^2=0.1…"
    plink \
      --bfile "$BFILE" \
      --maf 0.01 \
      --geno 0.01 \
      --hwe 1e-6 midp \
      --indep-pairwise 1000 50 0.1 \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --memory 60000 \
      --out "${S1_OUT}/ldtrain_tight"

    plink \
      --bfile "$BFILE" \
      --extract "${S1_OUT}/ldtrain_tight.prune.in" \
      --make-just-bim \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --memory 60000 \
      --out "${S1_OUT}/ldtrain_keep_tight"

    awk '{print $2}' "${S1_OUT}/ldtrain_keep_tight.bim" > "$TRAIN_SNPS"
    N_KEEP=$(wc -l < "$TRAIN_SNPS" || echo 0)
    diag "Tightened pruned count: ${N_KEEP}"
  fi
fi

[[ -s "$TRAIN_SNPS" ]] || die "Failed to build pruned training SNP list."
N_KEEP=$(wc -l < "$TRAIN_SNPS")
diag "Using TRAIN_SNPS=$TRAIN_SNPS (n=${N_KEEP})"
(( N_KEEP <= 1000000 )) || die "Training SNPs still >1,000,000 (n=${N_KEEP}); adjust pruning."

# ----------------------- REGENIE STEP 1 -----------------------
diag "Running Step-1 for ${N_COLS} imaging phenotypes"
diag "PHENO_FILE = $PHENO_FILE"
diag "PHENOS     = $PHENO_LIST_COMMA"
diag "BFILE      = $BFILE"
diag "OUT        = $S1_OUT"

regenie \
  --step 1 \
  --bed "$BFILE" \
  --extract "$TRAIN_SNPS" \
  --phenoFile "$PHENO_FILE" \
  --phenoColList "$PHENO_LIST_COMMA" \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "$LM_PREFIX" \
  --gz \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out "${S1_OUT}/s1"

# ----------------------- ENSURE s1_pred.list EXISTS -----------------------
if compgen -G "${S1_OUT}/s1_pred.list*" > /dev/null; then
  diag "Found Step-1 pred list(s) in ${S1_OUT}"
else
  if ls "${S1_OUT}"/s1_*_pred.list* >/dev/null 2>&1; then
    ls "${S1_OUT}"/s1_*_pred.list* > "${S1_OUT}/s1_pred.list"
    diag "Wrote ${S1_OUT}/s1_pred.list from chunked preds"
  else
    die "No Step-1 prediction files found in ${S1_OUT}; verify Step-1 completed."
  fi
fi

diag "Step-1 imaging done. Outputs in ${S1_OUT}"

